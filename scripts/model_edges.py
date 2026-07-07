#!/usr/bin/env python3
"""Model 1 — uncertainty-aware effect network (normal-normal EB). OPTIONAL / bonus.

Runs only if the spike (model_edges_spike.py) declared remote access VIABLE.
Reads by *slice* from S3 the rows of the candidate regulators (top of Model 2) —
it does NOT download the 17 GB h5ad — and runs EB over log_fc/lfcSE to produce edges
with effect probability.

Capped by default at N=8 regulators (~40s at 4.5s/row measured by the spike). Scale with --n.

    pip install h5py s3fs fsspec
    python scripts/model_edges.py --n 8

Output: docs/tables/robust_edges.csv
"""
import argparse
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "docs" / "tables"
S3_URL = "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"
LOG2_1P5 = np.log2(1.5)


def read_col(grp, name, h5py):
    """Decodes an AnnData obs/var column (categorical or string dataset)."""
    node = grp[name]
    if isinstance(node, h5py.Group):                       # categorical
        cats = [c.decode() if isinstance(c, bytes) else c for c in node["categories"][:]]
        codes = node["codes"][:]
        return np.array([cats[i] if i >= 0 else None for i in codes], dtype=object)
    arr = node[:]
    return np.array([x.decode() if isinstance(x, bytes) else x for x in arr], dtype=object)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=8, help="number of candidate regulators (rows to read)")
    args = ap.parse_args()

    hub = TAB / "hub_ranking_bayes.csv"
    if not hub.exists():
        sys.exit("Missing docs/tables/hub_ranking_bayes.csv — run first: python scripts/model_hubs.py")
    ranking = pd.read_csv(hub)

    try:
        import h5py, fsspec, s3fs  # noqa: F401
    except Exception as e:
        sys.exit(f"[optional] missing libraries ({e}). pip install h5py s3fs fsspec")
    import h5py, fsspec

    # candidates: peak condition per gene, significant KD only
    de = pd.read_csv(ROOT / "data" / "suppl_tables" / "DE_stats.suppl_table.csv")
    de = de[de["ontarget_significant"].astype(bool)]
    peak = (de.sort_values("n_downstream", ascending=False)
              .drop_duplicates("target_contrast_gene_name")
              .set_index("target_contrast_gene_name"))
    cand_genes = [g for g in ranking["target_contrast_gene_name"] if g in peak.index][:args.n]

    print(f"== Model 1 · EB over edges · {len(cand_genes)} regulators ==")
    t0 = time.time()
    f = fsspec.open(S3_URL, anon=True).open()
    h5 = h5py.File(f, "r")
    obs, var = h5["obs"], h5["var"]
    gene_obs = read_col(obs, "target_contrast_gene_name", h5py)
    cond_obs = read_col(obs, "culture_condition", h5py)
    var_names = read_col(var, "gene_name" if "gene_name" in var else var.attrs.get("_index", "_index"), h5py)
    row_of = {(g, c): i for i, (g, c) in enumerate(zip(gene_obs, cond_obs))}
    log_fc, lfcSE = h5["layers"]["log_fc"], h5["layers"]["lfcSE"]
    print(f"  h5ad opened + indexes in {time.time()-t0:.1f}s (streaming)")

    recs = []
    for g in cand_genes:
        c = peak.loc[g, "culture_condition"]
        i = row_of.get((g, c))
        if i is None:
            continue
        t = time.time()
        y = log_fc[i, :].astype(float)
        se = lfcSE[i, :].astype(float)
        print(f"  {g:10s} ({c:8s}) fila {i:>6} en {time.time()-t:4.1f}s")
        recs.append((g, c, y, se))
    h5.close()

    # tau2 by method of moments over the candidate edges (documented: sample biased toward high signal)
    allY = np.concatenate([r[2] for r in recs])
    allSE = np.concatenate([r[3] for r in recs])
    ok = np.isfinite(allY) & np.isfinite(allSE) & (allSE > 0)
    tau2 = max(float(np.var(allY[ok]) - np.median(allSE[ok] ** 2)), 1e-6)
    print(f"  tau2 (prior, from candidate edges) = {tau2:.4f}")

    rows = []
    for g, c, y, se in recs:
        v = 1.0 / (1.0 / tau2 + 1.0 / se ** 2)
        m = v * (y / se ** 2)
        sd = np.sqrt(v)
        p_pos = norm.sf(0, loc=m, scale=sd)
        p_big = norm.sf(LOG2_1P5, loc=m, scale=sd) + norm.cdf(-LOG2_1P5, loc=m, scale=sd)
        keep = p_big > 0.8
        for j in np.where(keep)[0]:
            rows.append({
                "perturbed_gene": g, "condition": c, "measured_gene": var_names[j],
                "log_fc": round(float(y[j]), 3), "lfcSE": round(float(se[j]), 3),
                "theta_post_mean": round(float(m[j]), 3), "theta_post_sd": round(float(sd[j]), 3),
                "p_effect_positive": round(float(p_pos[j]), 3),
                "p_abs_effect_gt_1p5x": round(float(p_big[j]), 3),
            })
    edges = pd.DataFrame(rows).sort_values("p_abs_effect_gt_1p5x", ascending=False)
    out = TAB / "robust_edges.csv"
    edges.to_csv(out, index=False)
    print(f"\n  robust edges (P(|effect|>1.5x)>0.8): {len(edges):,}")
    print(f"  table → {out.relative_to(ROOT)}")
    print(edges.head(10).to_string(index=False))
    print("\n✓ Model 1 (bonus) complete.")


if __name__ == "__main__":
    main()
