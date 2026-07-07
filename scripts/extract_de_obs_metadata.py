#!/usr/bin/env python3
"""Extracts ONLY the .obs of GWCD4i.DE_stats.h5ad (reproducibility metadata).

Closes the ranking caveat: replaces the cross-condition proxy with REAL cross-guide/
cross-donor reproducibility. The .obs is small (~34k metadata rows); the 17 GB .layers
are NOT read. Remote via S3 if not local; idempotent.

    pip install h5py s3fs fsspec
    python scripts/extract_de_obs_metadata.py [--force]

Output: docs/tables/de_obs_reproducibility_metadata.csv
"""
import argparse
import sys
import time
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "docs" / "tables"
TAB.mkdir(parents=True, exist_ok=True)
OUT = TAB / "de_obs_reproducibility_metadata.csv"
S3_URL = "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"

COLUMNS = [
    "target_contrast", "target_contrast_gene_name", "culture_condition",
    "single_guide_estimate", "n_guides",
    "guide_correlation_all", "guide_correlation_signif", "guide_n_signif_ontarget",
    "donor_correlation_all_mean", "donor_correlation_hits_mean", "donor_correlation_hits_min",
]


def read_col(grp, name, h5py):
    node = grp[name]
    if isinstance(node, h5py.Group):                       # AnnData categorical
        cats = [c.decode() if isinstance(c, bytes) else c for c in node["categories"][:]]
        codes = node["codes"][:]
        return np.array([cats[i] if i >= 0 else None for i in codes], dtype=object)
    arr = node[:]
    if arr.dtype.kind == "S" or (arr.dtype == object):
        return np.array([x.decode() if isinstance(x, bytes) else x for x in arr], dtype=object)
    return arr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="re-download even if the CSV exists")
    args = ap.parse_args()

    if OUT.exists() and not args.force:
        print(f"{OUT.relative_to(ROOT)} already exists — use --force to re-download. (idempotent)")
        sys.exit(0)

    try:
        import h5py, fsspec, s3fs  # noqa: F401
    except Exception as e:
        print(f"NOT VIABLE: missing libraries ({e}). pip install h5py s3fs fsspec")
        sys.exit(0)
    import h5py, fsspec

    print(f"== Extracting .obs (reproducibility metadata) ==\n  source: {S3_URL}")
    t0 = time.time()
    try:
        f = fsspec.open(S3_URL, anon=True).open()
        h5 = h5py.File(f, "r")
        obs = h5["obs"]
    except Exception as e:
        print(f"NOT VIABLE: could not open the remote h5ad ({type(e).__name__}: {e})")
        sys.exit(0)

    present = [c for c in COLUMNS if c in obs]
    missing = [c for c in COLUMNS if c not in obs]
    if missing:
        print(f"  [warning] columns missing from .obs: {missing}")
    data = {}
    for c in present:
        t = time.time()
        data[c] = read_col(obs, c, h5py)
        print(f"    {c:32s} read in {time.time()-t:4.1f}s")
    h5.close()

    df = pd.DataFrame(data)
    df.to_csv(OUT, index=False)
    dt = time.time() - t0
    size_kb = OUT.stat().st_size / 1024
    print(f"\n  rows={len(df):,} · cols={len(present)} · {size_kb:.0f} KB · {dt:.1f}s (.obs only, no .layers)")
    print(f"  VIABLE — table → {OUT.relative_to(ROOT)}")
    # coverage summary of the key metrics
    for c in ["single_guide_estimate", "guide_correlation_all", "donor_correlation_hits_mean"]:
        if c in df:
            nn = df[c].notna().sum()
            print(f"    {c:32s} non-null: {nn:,}/{len(df):,} ({nn/len(df)*100:.0f}%)")


if __name__ == "__main__":
    main()
