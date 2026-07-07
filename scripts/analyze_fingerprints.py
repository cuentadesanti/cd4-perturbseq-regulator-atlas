#!/usr/bin/env python3
"""Side analysis — mapa de programas regulatorios a partir de huellas transcriptómicas.

Pasa del ranking ("quién es fuerte") al mapa de programas ("qué tipo de efecto produce
cada perturbación"). Usa la matriz cacheada layers/log_fc: filas = perturbación×condición,
columnas = genes medidos, valores = log fold-change.

Análisis:
  1. Matriz de huellas por regulador (peak condition, KD on-target significativo).
  2. Similitud coseno entre reguladores + clustering jerárquico.
  3. PCA/SVD: varianza explicada y test "¿PC1 es sólo magnitud?" (corr con n_downstream).
  4. Validación: ¿miembros del mismo complejo (SAGA / Mediador / TCR) tienen huellas
     más parecidas que gene-sets aleatorios? (null por permutación).
  5. Global vs context-specific: para reguladores en >1 condición, coseno de su huella
     entre condiciones.

    python scripts/analyze_fingerprints.py

Salidas: docs/tables/fingerprint_*.csv, docs/figures/fingerprint_*.png
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parent.parent
CACHE = Path("/Users/cuentadesanti/code/hackaton/data/cache")
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
RNG = np.random.RandomState(0)

# Complejos conocidos (para validación). Se intersectan con reguladores presentes.
COMPLEXES = {
    "SAGA": ["TADA1", "TADA2A", "TADA2B", "TADA3", "SUPT20H", "SUPT7L", "TAF5L",
             "TAF6L", "SGF29", "ATXN7", "ATXN7L3", "USP22", "ENY2", "KAT2A",
             "KAT2B", "SUPT3H"],
    "Mediator": ["MED1", "MED4", "MED6", "MED7", "MED8", "MED9", "MED10", "MED11",
                 "MED12", "MED13", "MED14", "MED15", "MED16", "MED17", "MED18",
                 "MED19", "MED20", "MED21", "MED22", "MED23", "MED24", "MED25",
                 "MED26", "MED27", "MED28", "MED29", "MED30", "MED31", "CCNC",
                 "CDK8", "CDK19"],
    "TCR": ["ZAP70", "LCK", "LAT", "CD3D", "CD3E", "CD3G", "CD247", "FYN", "ITK",
            "PLCG1", "LCP2", "VAV1", "PIK3CD", "PRKCQ", "CARD11", "BCL10", "MALT1"],
    "SWI/SNF": ["SMARCA4", "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCD2",
                "SMARCE1", "ARID1A", "ARID1B", "ARID2", "PBRM1", "DPF2", "BCL7A"],
}

N_TOP = 200          # reguladores para similitud/PCA (por regpower/n_downstream)
MIN_DOWNSTREAM = 20  # descartar huellas casi vacías
N_PERM = 5000        # permutaciones para el null de complejos


def cosine_sim(M):
    n = np.linalg.norm(M, axis=1, keepdims=True)
    n[n == 0] = 1.0
    U = M / n
    return U @ U.T


def load():
    log_fc = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv")
    var = pd.read_csv(CACHE / "fingerprint_var.csv")
    obs["row"] = np.arange(len(obs))
    return log_fc, obs, var


def main():
    log_fc, obs, var = load()
    genes = var["gene_name"].values
    print(f"log_fc {log_fc.shape} · obs {obs.shape} · genes {len(genes)}")

    # --- filtrar a KD on-target significativo ---
    sig = obs[obs["ontarget_significant"].astype(str).str.lower().eq("true")].copy()
    print(f"filas con KD significativo: {len(sig)}  · reguladores únicos: "
          f"{sig['target_contrast_gene_name'].nunique()}")

    # --- huella por regulador en su PEAK condition ---
    peak = (sig.sort_values("n_downstream", ascending=False)
               .drop_duplicates("target_contrast_gene_name"))
    peak = peak[peak["n_downstream"] >= MIN_DOWNSTREAM]
    # top-N por efecto UNION miembros de complejos conocidos (para validación/PCA)
    complex_members = {g for members in COMPLEXES.values() for g in members}
    top = peak.sort_values("n_downstream", ascending=False).head(N_TOP)
    extra = peak[peak["target_contrast_gene_name"].isin(complex_members)
                 & ~peak["target_contrast_gene_name"].isin(top["target_contrast_gene_name"])]
    peak = pd.concat([top, extra]).sort_values("n_downstream", ascending=False)
    rows = peak["row"].values
    M = np.asarray(log_fc[rows, :], dtype=np.float32)   # (N_TOP, n_genes)
    reg = peak["target_contrast_gene_name"].values
    cond = peak["culture_condition"].values
    ndown = peak["n_downstream"].values.astype(float)
    print(f"matriz de huellas: {M.shape} · reguladores {len(reg)}")

    # descartar columnas (genes) sin variación entre estos reguladores
    keep = M.std(axis=0) > 0
    Mv = M[:, keep]
    genes_v = genes[keep]
    print(f"genes con varianza: {Mv.shape[1]}")

    # --- (2) similitud + clustering ---
    S = cosine_sim(Mv)
    D = 1.0 - S
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    clusters = fcluster(Z, t=0.85, criterion="distance")
    sim_df = pd.DataFrame(S, index=reg, columns=reg)
    sim_df.to_csv(TAB / "fingerprint_similarity.csv")

    reg_df = pd.DataFrame({
        "regulator": reg, "peak_condition": cond, "n_downstream": ndown.astype(int),
        "cluster": clusters,
    })

    # --- (3) PCA ---
    pca = PCA(n_components=min(20, Mv.shape[0]), random_state=0)
    scores = pca.fit_transform(Mv - Mv.mean(0))
    evr = pca.explained_variance_ratio_
    for k in range(min(5, scores.shape[1])):
        reg_df[f"PC{k+1}"] = scores[:, k]
    # test: ¿PC1 es sólo magnitud del efecto?
    from scipy.stats import spearmanr, pearsonr
    fp_norm = np.linalg.norm(Mv, axis=1)
    pc1_ndown_s = spearmanr(np.abs(scores[:, 0]), ndown).correlation
    pc1_norm_s = spearmanr(np.abs(scores[:, 0]), fp_norm).correlation
    print(f"varianza explicada PC1..5: {np.round(evr[:5], 3)}")
    print(f"|PC1| vs n_downstream (spearman): {pc1_ndown_s:.3f}")
    print(f"|PC1| vs ||huella|| (spearman):   {pc1_norm_s:.3f}")
    reg_df.to_csv(TAB / "fingerprint_regulators.csv", index=False)

    # top loadings de los primeros componentes (genes que definen cada eje)
    load_recs = []
    for k in range(min(6, pca.components_.shape[0])):
        comp = pca.components_[k]
        for sign, order in [("+", np.argsort(comp)[::-1][:15]),
                            ("-", np.argsort(comp)[:15])]:
            for gi in order:
                load_recs.append({"PC": k + 1, "sign": sign,
                                  "gene": genes_v[gi], "loading": round(float(comp[gi]), 4)})
    pd.DataFrame(load_recs).to_csv(TAB / "fingerprint_pc_loadings.csv", index=False)

    # --- (4) validación por complejos ---
    reg_set = set(reg)
    idx_of = {g: i for i, g in enumerate(reg)}
    comp_recs = []
    for name, members in COMPLEXES.items():
        present = [g for g in members if g in reg_set]
        if len(present) < 2:
            comp_recs.append({"complex": name, "n_present": len(present),
                              "mean_intra_cosine": None, "null_mean": None,
                              "z": None, "p_perm": None, "members": ";".join(present)})
            continue
        ii = [idx_of[g] for g in present]
        sub = S[np.ix_(ii, ii)]
        triu = sub[np.triu_indices(len(ii), k=1)]
        obs_mean = float(triu.mean())
        # null: gene-sets aleatorios del mismo tamaño
        null = np.empty(N_PERM)
        allidx = np.arange(len(reg))
        for p in range(N_PERM):
            jj = RNG.choice(allidx, size=len(ii), replace=False)
            s2 = S[np.ix_(jj, jj)]
            null[p] = s2[np.triu_indices(len(jj), k=1)].mean()
        z = (obs_mean - null.mean()) / (null.std() + 1e-9)
        p_perm = float((np.sum(null >= obs_mean) + 1) / (N_PERM + 1))
        comp_recs.append({"complex": name, "n_present": len(present),
                          "mean_intra_cosine": round(obs_mean, 4),
                          "null_mean": round(float(null.mean()), 4),
                          "z": round(float(z), 2), "p_perm": round(p_perm, 5),
                          "members": ";".join(present)})
        print(f"  [{name}] n={len(present)} intra={obs_mean:.3f} null={null.mean():.3f} "
              f"z={z:.2f} p={p_perm:.4f}")
    pd.DataFrame(comp_recs).to_csv(TAB / "fingerprint_complex_validation.csv", index=False)

    # --- (5) global vs context-specific ---
    # para reguladores presentes (KD sig) en >=2 condiciones: coseno entre condiciones
    ctx_recs = []
    for g, grp in sig.groupby("target_contrast_gene_name"):
        grp = grp[grp["n_downstream"] >= MIN_DOWNSTREAM]
        if grp["culture_condition"].nunique() < 2:
            continue
        rws = grp["row"].values
        Fg = np.asarray(log_fc[rws, :], dtype=np.float32)
        Sg = cosine_sim(Fg)
        conds = grp["culture_condition"].values
        pairs = Sg[np.triu_indices(len(conds), k=1)]
        ctx_recs.append({"regulator": g,
                         "conditions": ";".join(conds),
                         "n_conditions": len(conds),
                         "mean_cross_condition_cosine": round(float(pairs.mean()), 4),
                         "min_cross_condition_cosine": round(float(pairs.min()), 4),
                         "max_n_downstream": int(grp["n_downstream"].max())})
    ctx_df = pd.DataFrame(ctx_recs).sort_values("mean_cross_condition_cosine")
    ctx_df.to_csv(TAB / "fingerprint_context_specificity.csv", index=False)
    print(f"reguladores en >=2 condiciones: {len(ctx_df)}")

    # ---------- figuras ----------
    _fig_dendrogram(Z, reg, cond, clusters)
    _fig_pca(scores, evr, reg, clusters, ndown, COMPLEXES, idx_of, reg_set)
    _fig_pc1_test(scores, ndown, pc1_ndown_s)
    _fig_context(ctx_df)

    summary = {
        "n_regulators": int(len(reg)),
        "n_genes_variable": int(Mv.shape[1]),
        "evr_pc1_5": [round(float(x), 4) for x in evr[:5]],
        "pc1_abs_vs_ndownstream_spearman": round(float(pc1_ndown_s), 3),
        "complex_validation": comp_recs,
        "n_context_regulators": int(len(ctx_df)),
        "most_context_specific": ctx_df.head(10)["regulator"].tolist(),
        "most_global": ctx_df.tail(10)["regulator"].tolist(),
    }
    (TAB / "fingerprint_summary.json").write_text(json.dumps(summary, indent=2))
    print("OK · tablas y figuras escritas")


def _fig_dendrogram(Z, reg, cond, clusters):
    fig, ax = plt.subplots(figsize=(10, 14))
    labels = [f"{r}  [{c}]" for r, c in zip(reg, cond)]
    dendrogram(Z, labels=labels, orientation="left", ax=ax, leaf_font_size=5,
               color_threshold=0.85)
    ax.set_title("Clustering jerárquico de huellas transcriptómicas (coseno, avg linkage)")
    ax.set_xlabel("distancia coseno")
    fig.tight_layout()
    fig.savefig(FIG / "fingerprint_dendrogram.png", dpi=140)
    plt.close(fig)


def _fig_pca(scores, evr, reg, clusters, ndown, complexes, idx_of, reg_set):
    fig, ax = plt.subplots(figsize=(11, 9))
    ax.scatter(scores[:, 0], scores[:, 1], c=clusters, cmap="tab20", s=30,
               alpha=0.6, edgecolor="none")
    # resaltar complejos conocidos
    colors = {"SAGA": "#d62728", "Mediator": "#1f77b4", "TCR": "#2ca02c",
              "SWI/SNF": "#9467bd"}
    for name, members in complexes.items():
        pts = [idx_of[g] for g in members if g in reg_set]
        if not pts:
            continue
        ax.scatter(scores[pts, 0], scores[pts, 1], s=90, facecolor="none",
                   edgecolor=colors.get(name, "k"), linewidth=1.8, label=name)
        for gi in pts:
            ax.annotate(reg[gi], (scores[gi, 0], scores[gi, 1]), fontsize=6,
                        color=colors.get(name, "k"))
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}%)")
    ax.set_title("PCA del espacio de perturbaciones (color = cluster; anillos = complejo)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "fingerprint_pca.png", dpi=140)
    plt.close(fig)


def _fig_pc1_test(scores, ndown, rho):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(ndown, np.abs(scores[:, 0]), s=25, alpha=0.6)
    ax.set_xlabel("n_downstream (magnitud del efecto)")
    ax.set_ylabel("|PC1 score|")
    ax.set_title(f"¿PC1 es sólo magnitud?  spearman(|PC1|, n_downstream) = {rho:.2f}")
    fig.tight_layout()
    fig.savefig(FIG / "fingerprint_pc1_vs_magnitude.png", dpi=140)
    plt.close(fig)


def _fig_context(ctx_df):
    if len(ctx_df) == 0:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(ctx_df["mean_cross_condition_cosine"], bins=30, color="#4c78a8")
    ax.set_xlabel("coseno medio entre condiciones (huella)")
    ax.set_ylabel("nº reguladores")
    ax.set_title("Global (coseno alto) vs context-specific (coseno bajo)")
    fig.tight_layout()
    fig.savefig(FIG / "fingerprint_context_specificity.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
