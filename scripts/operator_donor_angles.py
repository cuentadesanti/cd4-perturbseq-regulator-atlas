#!/usr/bin/env python3
"""Step 4 — are the gene programs donor-reproducible AS SUBSPACES?

Principal angles between top-k gene-program subspaces of DISJOINT donor pairs only
(the 6 pairwise modalities share donors, so overlapping pairs inflate overlap).
Donors {1,2,3,4}: (1,2)vs(3,4), (1,3)vs(2,4), (1,4)vs(2,3). Per-donor matrices are
NOT in the local cache; absent => NEEDS-DATA + exit 0.

    python scripts/operator_donor_angles.py --k 5
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _opkernels as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"; FIG = ROOT / "docs" / "figures"
DISJOINT = [((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3))]
COLS = ["pair_a", "pair_b", "k", "mean_cos2", "null_mean_cos2", "null_p95", "above_null"]


def load_pairs():
    if not (CACHE / "by_donors_index.csv").exists():
        return None
    out = {}
    for a in range(1, 5):
        for b in range(a + 1, 5):
            p = CACHE / f"by_donors_pair_{a}{b}.npy"
            if p.exists():
                out[(a, b)] = np.load(p)
    return out or None


def top_subspace(M, k):
    Mc = M - M.mean(axis=0, keepdims=True)
    return np.linalg.svd(Mc, full_matrices=False)[2][:k].T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--null-n", type=int, default=500)
    args = ap.parse_args()
    TAB.mkdir(exist_ok=True, parents=True)
    mats = load_pairs()
    if mats is None:
        print("[NEEDS-DATA] per-donor-pair matrices absent. Expected:\n"
              "  data/cache/by_donors_index.csv (regulator,gene) AND\n"
              "  data/cache/by_donors_pair_12.npy ... by_donors_pair_34.npy (R x G, z-score)\n"
              "Produce via a fetch analogous to scripts/fetch_fingerprint_matrix.py over the\n"
              "per-donor layers. Writing empty table, exit 0.")
        pd.DataFrame(columns=COLS).to_csv(TAB / "operator_donor_angles.csv", index=False)
        return
    G = next(iter(mats.values())).shape[1]
    nm, np95 = op.random_subspace_null(G, args.k, n=args.null_n)
    rows = []
    for a, b in DISJOINT:
        if a in mats and b in mats:
            c2 = op.principal_angles(top_subspace(mats[a], args.k), top_subspace(mats[b], args.k))
            rows.append(dict(pair_a=f"{a[0]}{a[1]}", pair_b=f"{b[0]}{b[1]}", k=args.k,
                             mean_cos2=float(c2.mean()), null_mean_cos2=nm, null_p95=np95,
                             above_null=bool(c2.mean() > np95)))
    df = pd.DataFrame(rows, columns=COLS)
    df.to_csv(TAB / "operator_donor_angles.csv", index=False)
    FIG.mkdir(exist_ok=True, parents=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(df["pair_a"] + " vs " + df["pair_b"], df["mean_cos2"])
    ax.axhline(np95, ls="--", c="crimson", label="random-subspace 95th pct")
    ax.set_ylabel(f"mean cos²θ (top-{args.k})"); ax.legend()
    ax.set_title("Program-subspace overlap, disjoint donor pairs")
    fig.tight_layout(); fig.savefig(FIG / "36_operator_donor_angles.png", dpi=150)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
