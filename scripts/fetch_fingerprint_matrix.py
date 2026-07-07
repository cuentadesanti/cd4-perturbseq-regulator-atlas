#!/usr/bin/env python3
"""Downloads and caches the transcriptional fingerprint matrix (layers/log_fc).

The remote h5ad (17 GB) stores `layers/log_fc` as a DENSE, CONTIGUOUS, UNCOMPRESSED
matrix of shape (33983 perturbation×condition) x (10282 measured genes). Each row
is the *fingerprint* of a perturbation: the log fold-change of each downstream measured gene.

Because it is contiguous, a sequential read yields ~10 MB/s -> the full matrix
(~2.8 GB in float64) downloads in ~5 min. We cache it in float32 (1.4 GB) so the
downstream analysis (similarity / PCA / clustering) runs instantly and offline.

    pip install h5py s3fs fsspec
    python scripts/fetch_fingerprint_matrix.py

Outputs (in the main checkout; data/ is gitignored):
    data/cache/log_fc.f32.npy        (33983, 10282) float32
    data/cache/fingerprint_obs.csv   per-row metadata (perturbation, condition, QC)
    data/cache/fingerprint_var.csv   measured gene names (columns)
"""
import time
from pathlib import Path
import numpy as np
import pandas as pd

S3_URL = "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"
# data/ is gitignored and lives in the main checkout, not the worktree
CACHE = Path("/Users/cuentadesanti/code/hackaton/data/cache")
BLOCK = 2000  # rows per sequential read

# obs columns useful for the analysis (QC, reproducibility, counts)
OBS_COLS = [
    "index", "target_contrast_gene_name", "culture_condition", "target_contrast",
    "n_cells_target", "n_up_genes", "n_down_genes", "n_total_de_genes", "n_downstream",
    "ontarget_effect_size", "ontarget_significant", "offtarget_flag", "distal_offtarget_flag",
    "low_target_gex", "neighboring_gene_KD", "single_guide_estimate", "n_guides",
    "guide_correlation_all", "guide_correlation_signif",
    "donor_correlation_all_mean", "donor_correlation_hits_mean",
]


def read_col(grp, name, h5py):
    node = grp[name]
    if isinstance(node, h5py.Group):  # AnnData categorical
        cats = [c.decode() if isinstance(c, bytes) else c for c in node["categories"][:]]
        codes = node["codes"][:]
        return np.array([cats[i] if i >= 0 else None for i in codes], dtype=object)
    arr = node[:]
    if arr.dtype.kind in ("S", "O"):
        return np.array([x.decode() if isinstance(x, bytes) else x for x in arr], dtype=object)
    return arr


def main():
    import h5py, fsspec
    CACHE.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    f = fsspec.open(S3_URL, anon=True, default_cache_type="readahead").open()
    h5 = h5py.File(f, "r")
    print(f"[{time.time()-t0:5.1f}s] h5ad opened (streaming)")

    # --- obs (per-row metadata) ---
    obs = h5["obs"]
    avail = [c for c in OBS_COLS if c in obs]
    obs_df = pd.DataFrame({c: read_col(obs, c, h5py) for c in avail})
    obs_df.to_csv(CACHE / "fingerprint_obs.csv", index=False)
    print(f"[{time.time()-t0:5.1f}s] obs -> {obs_df.shape} saved")

    # --- var (measured genes = columns) ---
    var = h5["var"]
    vkey = "gene_name" if "gene_name" in var else "_index"
    gene_names = read_col(var, vkey, h5py)
    gene_ids = read_col(var, "gene_ids", h5py) if "gene_ids" in var else gene_names
    pd.DataFrame({"gene_name": gene_names, "gene_id": gene_ids}).to_csv(
        CACHE / "fingerprint_var.csv", index=False)
    print(f"[{time.time()-t0:5.1f}s] var -> {len(gene_names)} genes saved")

    # --- log_fc (dense matrix, sequential block read) ---
    lf = h5["layers"]["log_fc"]
    n, m = lf.shape
    out = np.empty((n, m), dtype=np.float32)
    for i in range(0, n, BLOCK):
        j = min(i + BLOCK, n)
        out[i:j] = lf[i:j, :].astype(np.float32)
        mb = out[:j].nbytes / 1e6
        print(f"[{time.time()-t0:5.1f}s]   rows {j:5d}/{n}  ({mb:6.0f} MB float32)")
    np.save(CACHE / "log_fc.f32.npy", out)
    print(f"[{time.time()-t0:5.1f}s] log_fc -> {out.shape} float32 saved "
          f"({out.nbytes/1e9:.2f} GB)")
    print(f"OK · total {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
