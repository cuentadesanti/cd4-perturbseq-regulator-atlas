#!/usr/bin/env python3
"""Side analysis — transcriptional fingerprints / perturbation program similarity.

Moves from "which regulators are strong" to "what transcriptomic programs each
perturbation induces and which regulators resemble each other". Does NOT replace the
core ranking, adds NO heavy models, does NOT download the full h5ad.

A perturbation's fingerprint = its vector of effects over the ~10k measured genes
(layers/zscore or layers/log_fc). Regulators with similar fingerprints act on the
same program. Over a balanced panel we run PCA, similarity, and clustering.

BALANCED panel (not pure top-EB, which would give pure chromatin/SAGA/Mediator):
    top global + top context-specific + promoted/demoted by the reproducibility audit.

    python scripts/analyze_fingerprints.py --n 200 --matrix zscore --top-genes 2000

Data source:
    --matrix log_fc  → local cache data/cache/log_fc.f32.npy (instant)
    --matrix zscore  → remote slice of layers/zscore for the panel rows (~few min),
                       cached in data/cache/panel_<matrix>_<hash>.npy

Outputs: docs/tables/fingerprint_*.csv|json · docs/figures/2{0..3}_fingerprint_*.png ·
         docs/FINGERPRINT_ANALYSIS.md (written separately).
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

# Known complexes — validation that the space recovers real structure.
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
    obs["row"] = np.arange(len(obs))  # position == row in cache/h5ad
    var = pd.read_csv(CACHE / "fingerprint_var.csv")
    return obs, var


def cosine_sim(M):
    n = np.linalg.norm(M, axis=1, keepdims=True)
    n[n == 0] = 1.0
    U = M / n
    return U @ U.T


# --------------------------------------------------------- panel selection
def build_panel(obs, n_cap):
    """Balanced panel: global + context-specific + promoted + demoted."""
    classes = _read("regulator_classes.csv")
    hub = _read("hub_ranking_bayes.csv")
    repro = _read("hub_ranking_bayes_reproducibility_aware.csv")
    audit = _read("reproducibility_audit.csv")
    if classes is None or hub is None:
        raise SystemExit("Missing regulator_classes.csv / hub_ranking_bayes.csv (run `make audit`).")

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

    # promoted / demoted: first the explicit audit labels,
    # then the largest rank_shift from the reproducibility-aware ranking
    if audit is not None:
        for g in audit.loc[audit.status == "promoted", "gene"]:
            picks.setdefault(g, "promoted")
        for g in audit.loc[audit.status == "demoted", "gene"]:
            picks.setdefault(g, "demoted")
    if repro is not None and "rank_shift" in repro.columns:
        rs = repro.set_index("gene")["rank_shift"]
        for g in rs.sort_values(ascending=False).head(25).index:
            picks.setdefault(g, "promoted")
        for g in rs.sort_values().head(25).index:
            picks.setdefault(g, "demoted")

    # (gene, condition) — peak condition = max n_downstream among significant KD in obs
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
    # cap: prioritize balanced coverage by source
    if len(panel) > n_cap:
        panel = (panel.sort_values(["source", "regpower_eb_mean"], ascending=[True, False])
                      .groupby("source", group_keys=False)
                      .head(int(np.ceil(n_cap / panel["source"].nunique())))
                      .head(n_cap))
    panel = panel.reset_index(drop=True)
    print(f"panel: {len(panel)} regulators · sources {panel['source'].value_counts().to_dict()}")
    return panel


# ------------------------------------------------------------ data reading
def read_matrix(panel, matrix):
    """Returns (M float32 [n_panel × n_genes]) for the requested metric."""
    rows = panel["row"].values
    if matrix == "log_fc" and (CACHE / "log_fc.f32.npy").exists():
        full = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
        return np.asarray(full[rows, :], dtype=np.float32)

    # remote slice of layers/<matrix>, cached by hash of the row set
    key = hashlib.md5((matrix + ",".join(map(str, sorted(rows)))).encode()).hexdigest()[:10]
    cpath = CACHE / f"panel_{matrix}_{key}.npy"
    if cpath.exists():
        print(f"  panel {matrix} from cache {cpath.name}")
        return np.load(cpath)
    import h5py, fsspec  # noqa
    print(f"  reading {len(rows)} rows from layers/{matrix} (remote slice)…")
    f = fsspec.open(S3_URL, anon=True, default_cache_type="readahead").open()
    h5 = h5py.File(f, "r")
    lay = h5["layers"][matrix]
    order = np.argsort(rows)
    M = np.empty((len(rows), lay.shape[1]), dtype=np.float32)
    M[order] = lay[np.sort(rows), :].astype(np.float32)
    np.save(cpath, M)
    return M


# --------------------------------------------------------------------- analysis
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200, help="cap on regulators in the panel")
    ap.add_argument("--matrix", choices=["zscore", "log_fc"], default="zscore")
    ap.add_argument("--top-genes", type=int, default=2000, help="genes by variance in the panel")
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
    print(f"raw matrix: {M.shape} ({args.matrix})")

    # genes (columns): top by variance in the panel
    v = M.var(axis=0)
    keep = np.argsort(v)[::-1][:args.top_genes]
    keep = keep[v[keep] > 0]
    Mg = M[:, keep]
    genes_v = genes_all[keep]
    # standardize columns (comparable across genes) before PCA/similarity
    mu, sd = Mg.mean(0), Mg.std(0)
    sd[sd == 0] = 1.0
    X = (Mg - mu) / sd
    print(f"standardized matrix: {X.shape} (top-{args.top_genes} genes by variance)")

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
    print(f"explained variance PC1..5: {np.round(evr[:5], 3)}")
    print(f"|PC1| vs n_downstream (spearman): {pc1_ndown:.3f}  (low = PC1 is not magnitude)")

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

    # ---- complex validation (permutation) ----
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

    # ---- program assignment: nearest known-complex centroid per regulator ----
    # Programs are ANCHORED to the permutation-validated complexes (SAGA/Mediator/TCR): each
    # regulator is assigned to the complex whose leave-one-out fingerprint centroid it most
    # resembles, IF cosine >= PROGRAM_COS_MIN and it beats the 2nd-best complex by >= MARGIN_MIN;
    # otherwise "mixed". This uses the same standardized space X as the validated cosine
    # similarity, and does NOT rely on the agnostic flat clustering (which merges the activation
    # programs into one blob). It is a nearest-prototype classifier against curated complexes —
    # transparent, conservative, and auditable via program_label_evidence.csv. "mixed" is expected
    # for the majority (we only have prototypes for 3 complexes).
    PROGRAM_OF = {"SAGA": "SAGA/chromatin", "Mediator": "Mediator/transcription", "TCR": "TCR signaling"}
    PROGRAM_COS_MIN, MARGIN_MIN, MARKER_CONSIST_MIN, K_MARK = 0.45, 0.05, 0.70, 12

    def _cos(a, b):
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        return float(a @ b / (na * nb)) if na > 0 and nb > 0 else 0.0

    complex_members = {n: [g for g in mem if g in reg_set] for n, mem in COMPLEXES.items()}
    complex_members = {n: m for n, m in complex_members.items() if len(m) >= 2}
    complex_member_set = set().union(*complex_members.values()) if complex_members else set()

    prog_label, near_complex, near_cos, near_margin = [], [], [], []
    for i, g in enumerate(reg):
        cosd = {}
        for n, mem in complex_members.items():
            idxs = [idx_of[m] for m in mem if m != g]     # leave-one-out for members
            cosd[n] = _cos(X[i], X[idxs].mean(0)) if idxs else -1.0
        ranked = sorted(cosd.items(), key=lambda kv: -kv[1])
        best, bestc = ranked[0]
        margin = bestc - (ranked[1][1] if len(ranked) > 1 else -1.0)
        near_complex.append(best); near_cos.append(round(bestc, 3)); near_margin.append(round(margin, 3))
        prog_label.append(PROGRAM_OF[best] if (bestc >= PROGRAM_COS_MIN and margin >= MARGIN_MIN) else "mixed")
    prog_label = np.array(prog_label, dtype=object)

    # ---- program fingerprint markers: convergent downstream RESPONSE genes per program ----
    # SEMANTICS: genes whose perturbation-RESPONSE z-scores (raw zscore Mg, relative to the panel)
    # are consistently high/low across the program's regulators — NOT baseline cell-type expression.
    prog_markers = {}
    for p in PROGRAM_OF.values():
        mask = prog_label == p
        if mask.sum() < 2:
            prog_markers[p] = []
            continue
        mg = Mg[mask].mean(0)
        consist = (np.sign(Mg[mask]) == np.sign(mg)).mean(0)
        cand = [(genes_v[gi], float(mg[gi]), float(consist[gi])) for gi in range(len(genes_v))
                if consist[gi] >= MARKER_CONSIST_MIN and mg[gi] != 0]
        cand.sort(key=lambda t: -abs(t[1]))
        prog_markers[p] = [t for t in cand if t[1] > 0][:K_MARK] + [t for t in cand if t[1] < 0][:K_MARK]
    pd.DataFrame([{"program": p, "direction": "up" if z > 0 else "down", "gene": g,
                   "mean_z": round(z, 3), "consistency": round(cons, 3)}
                  for p, ms in prog_markers.items() for g, z, cons in ms]
                 ).to_csv(TAB / "fingerprint_program_markers.csv", index=False)

    # ---- auditable label evidence (per program) ----
    # `assigned_neighbors` = non-curated genes PLACED in a program by fingerprint similarity —
    # candidate program neighbors, NOT a claim of physical complex membership.
    ev_rows = []
    for cx, p in PROGRAM_OF.items():
        mask = prog_label == p
        members_here = [reg[i] for i in range(len(reg)) if mask[i]]
        known = [g for g in members_here if g in complex_members.get(cx, [])]
        assigned = [g for g in members_here if g not in complex_member_set]
        cos_here = [near_cos[i] for i in range(len(reg)) if mask[i]]
        ev_rows.append({
            "program_label": p, "anchor_complex": cx, "n_regulators": int(mask.sum()),
            "n_known_complex_members": len(known), "known_members": ";".join(known),
            "assigned_neighbors": ";".join(assigned[:15]),
            "mean_centroid_cosine": round(float(np.mean(cos_here)), 3) if cos_here else None,
            "top_marker_genes": ";".join(f"{g}{'+' if z > 0 else '-'}" for g, z, _ in prog_markers[p][:8])})
    nmix = int((prog_label == "mixed").sum())
    ev_rows.append({"program_label": "mixed", "anchor_complex": "-", "n_regulators": nmix,
                    "n_known_complex_members": int(sum(prog_label[i] == "mixed" and reg[i] in complex_member_set
                                                       for i in range(len(reg)))),
                    "known_members": "", "assigned_neighbors": "", "mean_centroid_cosine": None,
                    "top_marker_genes": ""})
    pd.DataFrame(ev_rows).to_csv(TAB / "program_label_evidence.csv", index=False)
    progs = sorted(set(p for p in prog_label if p != "mixed"))
    n_lab = int((prog_label != "mixed").sum())
    print("  programs: " + ", ".join(f"{p}({int((prog_label == p).sum())})" for p in progs)
          + f" · mixed={nmix}")

    # ---- promoted/demoted neighborhood coherence (reported honestly, vs global reference) ----
    coh_rows = []
    for i, g in enumerate(reg):
        order = [j for j in np.argsort(S[i])[::-1] if j != i][:args.k_neighbors]
        top = order[0]
        coh_rows.append({
            "gene": g, "source": source[i], "regulator_class": rclass[i], "program_label": prog_label[i],
            "mean_knn_sim": round(float(np.mean([S[i, j] for j in order])), 3),
            "top_neighbor": reg[top], "top_neighbor_source": source[top],
            "top_neighbor_is_known_complex": bool(reg[top] in complex_member_set),
            "in_program": bool(prog_label[i] != "mixed")})
    coh_df = pd.DataFrame(coh_rows)
    coh_df.to_csv(TAB / "fingerprint_audit_coherence.csv", index=False)

    def _agg(src):
        d = coh_df[coh_df.source == src]
        return {} if d.empty else {"n": int(len(d)),
                                    "mean_knn_sim": round(float(d.mean_knn_sim.mean()), 3),
                                    "pct_in_program": round(float(d.in_program.mean()) * 100)}
    coherence_summary = {s: _agg(s) for s in ("global", "context-specific", "promoted", "demoted")}
    print("  coherence (mean_knn_sim · %in-program): "
          + " · ".join(f"{s}={coherence_summary[s]}" for s in coherence_summary if coherence_summary[s]))

    # ---- consolidated findings (one row per panel regulator) ----
    nn_by_gene = {}
    for r in nn_recs:
        nn_by_gene.setdefault(r["gene"], []).append((r["neighbor"], r["similarity"]))
    find_rows = []
    for i, g in enumerate(reg):
        p = prog_label[i]
        marks = prog_markers.get(p, [])[:6] if p != "mixed" else []
        find_rows.append({
            "gene": g, "condition": cond[i], "regulator_class": rclass[i], "source": source[i],
            "program_label": p, "nearest_complex": near_complex[i], "nearest_complex_cosine": near_cos[i],
            "margin_over_next": near_margin[i],
            "nearest_neighbors": ";".join(f"{n}:{s:.2f}" for n, s in nn_by_gene.get(g, [])[:3]),
            "program_markers": ";".join(f"{gn}{'+' if z > 0 else '-'}" for gn, z, _ in marks)})
    pd.DataFrame(find_rows).to_csv(TAB / "fingerprint_findings.csv", index=False)

    # ---- tablas ----
    pca_df = panel.copy()
    pca_df["cluster"] = clusters
    pca_df["program_label"] = prog_label
    for k in range(min(6, scores.shape[1])):
        pca_df[f"PC{k+1}"] = np.round(scores[:, k], 4)
    pca_df.to_csv(TAB / "fingerprint_pca_scores.csv", index=False)
    panel.to_csv(TAB / "fingerprint_panel.csv", index=False)
    nn_df.to_csv(TAB / "fingerprint_neighbors.csv", index=False)
    edges_df.to_csv(TAB / "fingerprint_similarity_edges.csv", index=False)

    # clusters: agnostic structural summary (members, dominant panel source)
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

    # ---- figures ----
    _fig_pca_by(scores, evr, cond, "condition", COND_COLOR,
                FIG / "20_fingerprint_pca_by_condition.png",
                "Fingerprint PCA · color = condition", reg, COMPLEXES, idx_of, reg_set)
    _fig_pca_by(scores, evr, rclass, "regulator_class", CLASS_COLOR,
                FIG / "21_fingerprint_pca_by_regulator_class.png",
                "Fingerprint PCA · color = class", reg, COMPLEXES, idx_of, reg_set)
    _fig_heatmap(S, reg, clusters, source, FIG / "22_fingerprint_similarity_heatmap.png")
    _fig_network(scores, edges_df, reg, source, idx_of,
                 FIG / "23_fingerprint_neighbor_network.png")
    PROGRAM_COLOR = {"SAGA/chromatin": "#d62728", "Mediator/transcription": "#1f77b4",
                     "TCR signaling": "#2ca02c", "mixed": "#c7c7c7"}
    _fig_pca_by(scores, evr, prog_label, "program", PROGRAM_COLOR,
                FIG / "24_fingerprint_pca_by_program.png",
                "Fingerprint similarity organizes perturbations into programs",
                reg, COMPLEXES, idx_of, reg_set)

    summary = {
        "matrix": args.matrix, "n_regulators": int(len(reg)),
        "n_genes": int(X.shape[1]), "top_genes": args.top_genes,
        "panel_sources": panel["source"].value_counts().to_dict(),
        "evr_pc1_5": [round(float(x), 4) for x in evr[:5]],
        "pc1_abs_vs_ndownstream_spearman": round(float(pc1_ndown), 3),
        "n_clusters": int(len(set(clusters))),
        "n_regulators_in_programs": int(n_lab),
        "programs": progs,
        "complex_validation": comp_recs,
        "audit_coherence": coherence_summary,
    }
    (TAB / "fingerprint_summary.json").write_text(json.dumps(summary, indent=2))
    print("OK · tables and figures written")


# ------------------------------------------------------------------- figures
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
    ax.set_title(title + "  (rings = known complex)")
    ax.legend(title=kind, loc="best")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def _fig_heatmap(S, reg, clusters, source, path):
    order = np.argsort(clusters)
    So = S[np.ix_(order, order)]
    fig, ax = plt.subplots(figsize=(11, 10))
    im = ax.imshow(So, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_title("Cosine similarity between fingerprints (ordered by cluster)")
    ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=ax, fraction=0.046, label="cosine similarity")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def _fig_network(scores, edges_df, reg, source, idx_of, path):
    # nodes at their PCA coords; edges = nearest neighbors (top similarities)
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
    ax.set_title("Transcriptomic neighbor network (PCA layout · color = panel source)")
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(loc="best")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


if __name__ == "__main__":
    main()
