#!/usr/bin/env python3
"""Side analysis — transcriptional fingerprints / perturbation program similarity.

Pasa de "qué reguladores son fuertes" a "qué programas transcriptómicos produce cada
perturbación y qué reguladores se parecen entre sí". NO reemplaza el ranking core, NO
agrega modelos pesados, NO descarga el h5ad completo.

La huella de una perturbación = su vector de efectos sobre los ~10k genes medidos
(layers/zscore o layers/log_fc). Reguladores con huellas parecidas actúan sobre el
mismo programa. Sobre un panel balanceado corremos PCA, similitud y clustering.

Panel BALANCEADO (no puro top-EB, que daría puro cromatina/SAGA/Mediator):
    top global + top context-specific + promovidos/demotados por la auditoría de repro.

    python scripts/analyze_fingerprints.py --n 200 --matrix zscore --top-genes 2000

Fuente de datos:
    --matrix log_fc  → cache local data/cache/log_fc.f32.npy (instantáneo)
    --matrix zscore  → slice remoto de layers/zscore para las filas del panel (~pocos min),
                       cacheado en data/cache/panel_<matrix>_<hash>.npy

Salidas: docs/tables/fingerprint_*.csv|json · docs/figures/2{0..3}_fingerprint_*.png ·
         docs/FINGERPRINT_ANALYSIS.md (se escribe aparte).
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parent.parent
CACHE = Path("/Users/cuentadesanti/code/hackaton/data/cache")
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
S3_URL = "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"
RNG = np.random.RandomState(0)

COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]
COND_COLOR = {"Rest": "#7f7f7f", "Stim8hr": "#ff7f0e", "Stim48hr": "#d62728"}
CLASS_COLOR = {"global": "#1f77b4", "condition-specific": "#2ca02c"}
SOURCE_COLOR = {"global": "#1f77b4", "context-specific": "#2ca02c",
                "promoted": "#9467bd", "demoted": "#8c564b"}

# Complejos conocidos — validación de que el espacio recupera estructura real.
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
}


# --------------------------------------------------------------------------- IO
def _read(name, base=TAB):
    p = base / name
    return pd.read_csv(p) if p.exists() else None


def load_obs():
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv")
    obs = obs.reset_index(drop=True)
    obs["row"] = np.arange(len(obs))  # posición == fila en cache/h5ad
    var = pd.read_csv(CACHE / "fingerprint_var.csv")
    return obs, var


def cosine_sim(M):
    n = np.linalg.norm(M, axis=1, keepdims=True)
    n[n == 0] = 1.0
    U = M / n
    return U @ U.T


# --------------------------------------------------------- selección del panel
def build_panel(obs, n_cap):
    """Panel balanceado: global + context-specific + promovidos + demotados."""
    classes = _read("regulator_classes.csv")
    hub = _read("hub_ranking_bayes.csv")
    repro = _read("hub_ranking_bayes_reproducibility_aware.csv")
    audit = _read("reproducibility_audit.csv")
    if classes is None or hub is None:
        raise SystemExit("Faltan regulator_classes.csv / hub_ranking_bayes.csv (corre `make audit`).")

    hub = hub.rename(columns={"target_contrast_gene_name": "gene"})
    regpower = hub.set_index("gene")["regpower_eb_mean"]
    cls = classes.set_index("gene")["regulator_class"]

    def by_class(name, k):
        genes = cls[cls == name].index
        s = regpower.reindex(genes).dropna().sort_values(ascending=False)
        return list(s.head(k).index)

    picks = {}
    for g in by_class("global", 75):
        picks.setdefault(g, "global")
    for g in by_class("condition-specific", 75):
        picks.setdefault(g, "context-specific")

    # promovidos / demotados: primero las etiquetas explícitas de la auditoría,
    # luego los mayores rank_shift del ranking reproducibility-aware
    if audit is not None:
        for g in audit.loc[audit.status == "promovido", "gene"]:
            picks.setdefault(g, "promoted")
        for g in audit.loc[audit.status == "demotado", "gene"]:
            picks.setdefault(g, "demoted")
    if repro is not None and "rank_shift" in repro.columns:
        rs = repro.set_index("gene")["rank_shift"]
        for g in rs.sort_values(ascending=False).head(25).index:
            picks.setdefault(g, "promoted")
        for g in rs.sort_values().head(25).index:
            picks.setdefault(g, "demoted")

    # (gene, condición) — condición pico = max n_downstream entre KD significativos en obs
    sig = obs[obs["ontarget_significant"].astype(str).str.lower().eq("true")]
    peak = (sig.sort_values("n_downstream", ascending=False)
               .drop_duplicates("target_contrast_gene_name")
               .set_index("target_contrast_gene_name"))
    rows = []
    for g, source in picks.items():
        if g not in peak.index:
            continue
        r = peak.loc[g]
        rows.append({"gene": g, "condition": r["culture_condition"], "row": int(r["row"]),
                     "source": source, "regulator_class": cls.get(g, None),
                     "regpower_eb_mean": float(regpower.get(g, np.nan)),
                     "n_downstream": int(r["n_downstream"])})
    panel = pd.DataFrame(rows).drop_duplicates("gene")
    # cap: prioriza cobertura balanceada por fuente
    if len(panel) > n_cap:
        panel = (panel.sort_values(["source", "regpower_eb_mean"], ascending=[True, False])
                      .groupby("source", group_keys=False)
                      .head(int(np.ceil(n_cap / panel["source"].nunique())))
                      .head(n_cap))
    panel = panel.reset_index(drop=True)
    print(f"panel: {len(panel)} reguladores · fuentes {panel['source'].value_counts().to_dict()}")
    return panel


# ------------------------------------------------------------ lectura de datos
def read_matrix(panel, matrix):
    """Devuelve (M float32 [n_panel × n_genes]) para la métrica pedida."""
    rows = panel["row"].values
    if matrix == "log_fc" and (CACHE / "log_fc.f32.npy").exists():
        full = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
        return np.asarray(full[rows, :], dtype=np.float32)

    # slice remoto de layers/<matrix>, cacheado por hash del set de filas
    key = hashlib.md5((matrix + ",".join(map(str, sorted(rows)))).encode()).hexdigest()[:10]
    cpath = CACHE / f"panel_{matrix}_{key}.npy"
    if cpath.exists():
        print(f"  panel {matrix} desde cache {cpath.name}")
        return np.load(cpath)
    import h5py, fsspec  # noqa
    print(f"  leyendo {len(rows)} filas de layers/{matrix} (slice remoto)…")
    f = fsspec.open(S3_URL, anon=True, default_cache_type="readahead").open()
    h5 = h5py.File(f, "r")
    lay = h5["layers"][matrix]
    order = np.argsort(rows)
    M = np.empty((len(rows), lay.shape[1]), dtype=np.float32)
    M[order] = lay[np.sort(rows), :].astype(np.float32)
    np.save(cpath, M)
    return M


# --------------------------------------------------------------------- análisis
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200, help="tope de reguladores en el panel")
    ap.add_argument("--matrix", choices=["zscore", "log_fc"], default="zscore")
    ap.add_argument("--top-genes", type=int, default=2000, help="genes por varianza en el panel")
    ap.add_argument("--k-neighbors", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    obs, var = load_obs()
    genes_all = var["gene_name"].values
    panel = build_panel(obs, args.n)
    M = read_matrix(panel, args.matrix)
    M = np.nan_to_num(M, nan=0.0, posinf=0.0, neginf=0.0)
    print(f"matriz cruda: {M.shape} ({args.matrix})")

    # genes (columnas): top por varianza en el panel
    v = M.var(axis=0)
    keep = np.argsort(v)[::-1][:args.top_genes]
    keep = keep[v[keep] > 0]
    Mg = M[:, keep]
    genes_v = genes_all[keep]
    # estandarizar columnas (comparables entre genes) antes de PCA/similitud
    mu, sd = Mg.mean(0), Mg.std(0)
    sd[sd == 0] = 1.0
    X = (Mg - mu) / sd
    print(f"matriz estandarizada: {X.shape} (top-{args.top_genes} genes por varianza)")

    reg = panel["gene"].values
    cond = panel["condition"].values
    source = panel["source"].values
    rclass = panel["regulator_class"].fillna("NA").values
    ndown = panel["n_downstream"].values.astype(float)
    idx_of = {g: i for i, g in enumerate(reg)}

    # ---- PCA / SVD ----
    from sklearn.decomposition import PCA
    pca = PCA(n_components=min(20, X.shape[0] - 1), random_state=0)
    scores = pca.fit_transform(X)
    evr = pca.explained_variance_ratio_
    from scipy.stats import spearmanr
    pc1_ndown = spearmanr(np.abs(scores[:, 0]), ndown).correlation
    print(f"varianza explicada PC1..5: {np.round(evr[:5], 3)}")
    print(f"|PC1| vs n_downstream (spearman): {pc1_ndown:.3f}  (bajo = PC1 no es magnitud)")

    # ---- similitud + clustering ----
    S = cosine_sim(X)
    D = 1.0 - S
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    n_clusters = min(12, max(2, len(reg) // 15))
    clusters = fcluster(Z, t=n_clusters, criterion="maxclust")

    # ---- nearest neighbors ----
    nn_recs = []
    for i, g in enumerate(reg):
        order = np.argsort(S[i])[::-1]
        order = [j for j in order if j != i][:args.k_neighbors]
        for rank, j in enumerate(order, 1):
            nn_recs.append({"gene": g, "condition": cond[i], "neighbor": reg[j],
                            "neighbor_condition": cond[j], "similarity": round(float(S[i, j]), 4),
                            "rank": rank})
    nn_df = pd.DataFrame(nn_recs)

    # ---- similarity edges (para red): top-k por nodo, deduplicado ----
    edge_set = {}
    for r in nn_recs:
        a, b = sorted([r["gene"], r["neighbor"]])
        edge_set[(a, b)] = max(edge_set.get((a, b), 0), r["similarity"])
    edges_df = pd.DataFrame([{"source_gene": a, "target_gene": b, "similarity": s}
                             for (a, b), s in edge_set.items()])
    edges_df = edges_df.sort_values("similarity", ascending=False)

    # ---- validación por complejos (permutación) ----
    reg_set = set(reg)
    comp_recs = []
    for name, members in COMPLEXES.items():
        present = [g for g in members if g in reg_set]
        if len(present) < 2:
            comp_recs.append({"complex": name, "n_present": len(present),
                              "mean_intra_cosine": None, "null_mean": None, "z": None,
                              "p_perm": None, "members": ";".join(present)})
            continue
        ii = [idx_of[g] for g in present]
        obs_mean = float(S[np.ix_(ii, ii)][np.triu_indices(len(ii), 1)].mean())
        null = np.empty(5000)
        allidx = np.arange(len(reg))
        for p in range(5000):
            jj = RNG.choice(allidx, size=len(ii), replace=False)
            null[p] = S[np.ix_(jj, jj)][np.triu_indices(len(jj), 1)].mean()
        z = (obs_mean - null.mean()) / (null.std() + 1e-9)
        p_perm = float((np.sum(null >= obs_mean) + 1) / 5001)
        comp_recs.append({"complex": name, "n_present": len(present),
                          "mean_intra_cosine": round(obs_mean, 4),
                          "null_mean": round(float(null.mean()), 4), "z": round(float(z), 2),
                          "p_perm": round(p_perm, 5), "members": ";".join(present)})
        print(f"  [{name}] n={len(present)} intra={obs_mean:.3f} null={null.mean():.3f} "
              f"z={z:.2f} p={p_perm:.4f}")

    # ---- tablas ----
    pca_df = panel.copy()
    pca_df["cluster"] = clusters
    for k in range(min(6, scores.shape[1])):
        pca_df[f"PC{k+1}"] = np.round(scores[:, k], 4)
    pca_df.to_csv(TAB / "fingerprint_pca_scores.csv", index=False)
    panel.to_csv(TAB / "fingerprint_panel.csv", index=False)
    nn_df.to_csv(TAB / "fingerprint_neighbors.csv", index=False)
    edges_df.to_csv(TAB / "fingerprint_similarity_edges.csv", index=False)

    # clusters: resumen por cluster (miembros, fuente dominante)
    cl_recs = []
    for c in sorted(set(clusters)):
        mem = reg[clusters == c]
        srcs = pd.Series(source[clusters == c]).value_counts().to_dict()
        cl_recs.append({"cluster": int(c), "size": int(len(mem)),
                        "dominant_source": max(srcs, key=srcs.get),
                        "members": ";".join(mem[:25])})
    pd.DataFrame(cl_recs).to_csv(TAB / "fingerprint_clusters.csv", index=False)

    pd.DataFrame(comp_recs).to_csv(TAB / "fingerprint_complex_validation.csv", index=False)

    load_recs = []
    for k in range(min(6, pca.components_.shape[0])):
        comp = pca.components_[k]
        for sign, order in [("+", np.argsort(comp)[::-1][:15]), ("-", np.argsort(comp)[:15])]:
            for gi in order:
                load_recs.append({"PC": k + 1, "sign": sign, "gene": genes_v[gi],
                                  "loading": round(float(comp[gi]), 4)})
    pd.DataFrame(load_recs).to_csv(TAB / "fingerprint_pc_loadings.csv", index=False)

    # ---- figuras ----
    _fig_pca_by(scores, evr, cond, "condition", COND_COLOR,
                FIG / "20_fingerprint_pca_by_condition.png",
                "PCA de huellas · color = condición", reg, COMPLEXES, idx_of, reg_set)
    _fig_pca_by(scores, evr, rclass, "regulator_class", CLASS_COLOR,
                FIG / "21_fingerprint_pca_by_regulator_class.png",
                "PCA de huellas · color = clase", reg, COMPLEXES, idx_of, reg_set)
    _fig_heatmap(S, reg, clusters, source, FIG / "22_fingerprint_similarity_heatmap.png")
    _fig_network(scores, edges_df, reg, source, idx_of,
                 FIG / "23_fingerprint_neighbor_network.png")

    summary = {
        "matrix": args.matrix, "n_regulators": int(len(reg)),
        "n_genes": int(X.shape[1]), "top_genes": args.top_genes,
        "panel_sources": panel["source"].value_counts().to_dict(),
        "evr_pc1_5": [round(float(x), 4) for x in evr[:5]],
        "pc1_abs_vs_ndownstream_spearman": round(float(pc1_ndown), 3),
        "n_clusters": int(len(set(clusters))),
        "complex_validation": comp_recs,
    }
    (TAB / "fingerprint_summary.json").write_text(json.dumps(summary, indent=2))
    print("OK · tablas y figuras escritas")


# ------------------------------------------------------------------- figuras
def _annotate_complexes(ax, scores, reg, complexes, idx_of, reg_set):
    ring = {"SAGA": "#d62728", "Mediator": "#1f77b4", "TCR": "#111"}
    for name, members in complexes.items():
        pts = [idx_of[g] for g in members if g in reg_set]
        if not pts:
            continue
        ax.scatter(scores[pts, 0], scores[pts, 1], s=90, facecolor="none",
                   edgecolor=ring.get(name, "k"), linewidth=1.6, zorder=5)
        for gi in pts:
            ax.annotate(reg[gi], (scores[gi, 0], scores[gi, 1]), fontsize=6, zorder=6)


def _fig_pca_by(scores, evr, labels, kind, cmap, path, title, reg, complexes, idx_of, reg_set):
    fig, ax = plt.subplots(figsize=(11, 9))
    for lab in pd.unique(labels):
        m = labels == lab
        ax.scatter(scores[m, 0], scores[m, 1], s=34, alpha=0.7, edgecolor="none",
                   color=cmap.get(lab, "#bbbbbb"), label=str(lab))
    _annotate_complexes(ax, scores, reg, complexes, idx_of, reg_set)
    ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}%)")
    ax.set_title(title + "  (anillos = complejo conocido)")
    ax.legend(title=kind, loc="best")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def _fig_heatmap(S, reg, clusters, source, path):
    order = np.argsort(clusters)
    So = S[np.ix_(order, order)]
    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.imshow(So, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_title("Similitud coseno entre huellas (ordenado por cluster)")
    ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.046, label="cosine similarity")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def _fig_network(scores, edges_df, reg, source, idx_of, path):
    # nodos en sus coords PCA; aristas = nearest neighbors (top similitudes)
    fig, ax = plt.subplots(figsize=(11, 9))
    top = edges_df.head(400)
    for r in top.itertuples():
        a, b = idx_of.get(r.source_gene), idx_of.get(r.target_gene)
        if a is None or b is None:
            continue
        ax.plot([scores[a, 0], scores[b, 0]], [scores[a, 1], scores[b, 1]],
                color="#cccccc", lw=0.4, alpha=0.6, zorder=1)
    for s in pd.unique(source):
        m = source == s
        ax.scatter(scores[m, 0], scores[m, 1], s=36, color=SOURCE_COLOR.get(s, "#999"),
                   label=str(s), zorder=3, edgecolor="white", linewidth=0.3)
    ax.set_title("Red de vecinos transcriptómicos (layout PCA · color = fuente del panel)")
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(loc="best")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


if __name__ == "__main__":
    main()
