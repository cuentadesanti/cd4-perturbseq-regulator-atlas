#!/usr/bin/env python3
"""Effect-size ranking (critique #1): rank perturbations on a continuous, power-robust
effect magnitude instead of n_downstream (a thresholded DE-call count).

Why this exists
---------------
n_downstream conflates biological magnitude with statistical power. But the naive fix
(L2 / sum|log2FC| over ALL genes) is WORSE: raw per-gene log-FC noise scales inversely
with cell number, so the all-gene norm correlates ~-0.7 with n_cells_target -- it measures
imprecision, not biology. The honest decomposition is breadth x intensity:

    breadth   = n_sig            (# genes with adj_p < FDR)  -- FDR-controlled, mild power dep.
    intensity = mean|log2FC| over the FDR-significant genes  -- power-decoupled
    magnitude = sum|log2FC| over the FDR-significant genes   -- total real displacement
    l2_z      = ||zscore||_2                                 -- precision-weighted (test-stat like)

We compute all of these by streaming layers/{zscore,adj_p_value} from the remote h5ad in
row-blocks (a few minutes; only small per-contrast summaries are kept, not the 2.8 GB layers),
using the locally cached log_fc.f32.npy for the log-FC values. We then report, honestly, how
each metric correlates with n_cells_target (the power-decoupling audit) and with the existing
n_downstream ranking (does the swap matter?).

    python scripts/rank_effect_size.py [--fdr 0.1] [--block 2000]

Outputs:
    docs/tables/effect_size_metric_audit.csv   per-metric rho vs n_downstream and n_cells
    docs/tables/effect_size_ranking.csv        per-gene ranking on the chosen magnitude
"""
import argparse
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
CACHE = Path("/Users/cuentadesanti/code/hackaton/data/cache")
TAB = ROOT / "docs" / "tables"
S3_URL = "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"


def stream_summaries(fdr, block):
    """Stream zscore + adj_p_value in row-blocks; return per-contrast magnitude summaries.
    Uses cached log_fc for the log-FC values (rows aligned 1:1 with the h5ad)."""
    import h5py, fsspec
    lfc = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
    n, g = lfc.shape
    cols = {k: np.empty(n, np.float64) for k in
            ["n_sig", "mag_sig", "mean_sig", "l2_sig", "l2_z", "sabs_z", "l2_all"]}
    t0 = time.time()
    f = fsspec.open(S3_URL, anon=True, default_cache_type="readahead").open()
    h5 = h5py.File(f, "r")
    Z, P = h5["layers"]["zscore"], h5["layers"]["adj_p_value"]
    for s in range(0, n, block):
        e = min(s + block, n)
        z = Z[s:e]                       # (b, g) float64
        p = P[s:e]
        L = np.abs(lfc[s:e].astype(np.float64))
        sig = p < fdr                    # FDR-significant genes (the real DE set)
        nsig = sig.sum(1)
        magsig = (L * sig).sum(1)
        cols["n_sig"][s:e] = nsig
        cols["mag_sig"][s:e] = magsig
        cols["mean_sig"][s:e] = np.where(nsig > 0, magsig / np.maximum(nsig, 1), np.nan)
        cols["l2_sig"][s:e] = np.sqrt(((L * sig) ** 2).sum(1))
        cols["l2_z"][s:e] = np.sqrt(np.nansum(z ** 2, axis=1))
        cols["sabs_z"][s:e] = np.nansum(np.abs(z), axis=1)
        cols["l2_all"][s:e] = np.sqrt((L ** 2).sum(1))
        print(f"  [{time.time()-t0:5.1f}s] rows {e}/{n}", end="\r", flush=True)
    print()
    return pd.DataFrame(cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fdr", type=float, default=0.1)
    ap.add_argument("--block", type=int, default=2000)
    ap.add_argument("--rank-on", default="mag_sig",
                    help="magnitude column to rank the per-gene table on")
    args = ap.parse_args()
    TAB.mkdir(parents=True, exist_ok=True)

    obs = pd.read_csv(CACHE / "fingerprint_obs.csv")
    summ = stream_summaries(args.fdr, args.block)
    df = pd.concat([obs.reset_index(drop=True), summ], axis=1)

    # sanity: our streamed DE set must equal theirs. n_sig(adj_p<fdr) matches n_total_de_genes
    # EXACTLY; n_downstream = n_total_de_genes - 1 (drops the on-target gene), so it is rank-identical
    # to n_sig (spearman 1.0) but offset by one per row. We verify both explicitly.
    if "n_total_de_genes" in df:
        ok = df.n_total_de_genes.notna() & df.n_sig.notna()
        frac_exact = float((df.n_sig[ok] == df.n_total_de_genes[ok]).mean())
        print(f"sanity  n_sig == n_total_de_genes: frac_equal={frac_exact:.4f} "
              f"(expect 1.0000 -> our FDR DE set is exactly theirs)")
    m = df.n_downstream.notna() & df.n_sig.notna() & (df.ontarget_significant == True)
    rho_nsig = spearmanr(df.n_sig[m], df.n_downstream[m]).correlation
    print(f"sanity  spearman(n_sig, n_downstream | sig-KD) = {rho_nsig:.6f} "
          f"(expect 1.000000 -> rank-identical over the ranking population)")

    # ---- power-decoupling audit: contrast-level, within significant-KD contrasts ----
    sig = df[df.ontarget_significant == True].copy()
    metrics = ["n_downstream", "n_sig", "mag_sig", "mean_sig", "l2_sig", "l2_z", "sabs_z", "l2_all"]
    audit = []
    for k in metrics:
        x = sig[k].to_numpy(float); ok = np.isfinite(x)
        audit.append({
            "metric": k,
            "rho_vs_n_downstream": round(spearmanr(x[ok], sig.n_downstream.to_numpy(float)[ok]).correlation, 3),
            "rho_vs_n_cells": round(spearmanr(x[ok], sig.n_cells_target.to_numpy(float)[ok]).correlation, 3),
        })
    audit = pd.DataFrame(audit)
    audit.to_csv(TAB / "effect_size_metric_audit.csv", index=False)
    print("\n=== metric audit (want |rho_vs_n_cells| small = power-decoupled) ===")
    print(audit.to_string(index=False))

    # ---- per-gene ranking on the chosen magnitude (KD-gated, peak across conditions) ----
    g = sig.groupby("target_contrast_gene_name")
    per = pd.DataFrame({
        "n_downstream_peak": g.n_downstream.max(),
        "mag_sig_peak": g.mag_sig.max(),
        "mean_sig_peak": g.mean_sig.max(),
        "l2_z_peak": g.l2_z.max(),
    }).dropna(subset=[args.rank_on + "_peak"])
    per["rank_magnitude"] = per[args.rank_on + "_peak"].rank(ascending=False).astype(int)
    per["rank_n_downstream"] = per["n_downstream_peak"].rank(ascending=False).astype(int)
    per = per.sort_values("rank_magnitude")
    per.to_csv(TAB / "effect_size_ranking.csv")

    rho = spearmanr(per["mag_sig_peak"], per["n_downstream_peak"]).correlation
    t_new = set(per.sort_values(args.rank_on + "_peak", ascending=False).head(30).index)
    t_old = set(per.sort_values("n_downstream_peak", ascending=False).head(30).index)
    print(f"\n=== ranking swap (per gene, n={len(per)}) ===")
    print(f"spearman(mag_sig, n_downstream) = {rho:.3f}")
    print(f"top-30 overlap = {len(t_new & t_old)}/30")
    print("dropped by magnitude:", sorted(t_old - t_new))
    print("promoted by magnitude:", sorted(t_new - t_old))
    print(f"\nwrote {TAB/'effect_size_metric_audit.csv'} and {TAB/'effect_size_ranking.csv'}")


if __name__ == "__main__":
    main()
