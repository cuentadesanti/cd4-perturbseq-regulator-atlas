#!/usr/bin/env python3
"""Inspects a .h5ad / .h5mu without loading the whole matrix into memory.

Usage:
    python scripts/inspect.py data/GWCD4i.DE_stats.h5ad
"""
import sys
import h5py


def summarize(path: str) -> None:
    print(f"== {path} ==")
    with h5py.File(path, "r") as f:
        def show(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(f"  {name:60s} {str(obj.shape):>18s} {obj.dtype}")
        # Only the first level so as not to flood the output
        for key in f.keys():
            print(f"[{key}]")
            f[key].visititems(show) if isinstance(f[key], h5py.Group) else None


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    summarize(sys.argv[1])
