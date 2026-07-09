#!/usr/bin/env python3
"""Step 5 (stretch) — square-block deconvolution + asymmetric subsumption.

HYPOTHESIS-GENERATING ONLY. A ≈ I - L^{-1} is valid only on the square block where
perturbed regulators are also readout genes (verified 91 here), is a linear approx
to a nonlinear system, is ill-conditioned (every edge carries the block condition
number), and in z-score space the steady-state reading is a further stretch.

    python scripts/operator_deconvolution.py --n-robust 50 --ridge 1e-2
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def deconvolve(L, ridge):
    n = L.shape[0]
    return np.eye(n) - np.linalg.inv(L + ridge * np.eye(n)), float(np.linalg.cond(L + ridge * np.eye(n)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-robust", type=int, default=50)
    ap.add_argument("--ridge", type=float, default=1e-2)
    ap.add_argument("--sig-thr", type=float, default=2.0)   # |z| threshold on the panel
    args = ap.parse_args()
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    tensor = d["tensor"]; regs, genes = d["regulators"].astype(str), d["genes"].astype(str)
    # Square block = regulators that are ALSO readout genes. Prefer the donor-robust
    # shortlist (deduplicated — top_robust_regulators.csv has repeated rows per
    # condition/contrast, which would otherwise seed duplicate self-pairs), then top up
    # from the panel's regulator-genes to reach n_robust so the block is a meaningful size.
    rp = TAB / "top_robust_regulators.csv"
    rset, gset = set(regs), set(genes)
    shortlist = (list(dict.fromkeys(pd.read_csv(rp).iloc[:, 0].astype(str).tolist()))
                 if rp.exists() else [])
    block = [g for g in shortlist if g in rset and g in gset]
    for g in regs:                       # top up from panel regulator-genes (regs are unique)
        if len(block) >= args.n_robust:
            break
        if g in gset and g not in block:
            block.append(g)
    block = block[: args.n_robust]
    if len(block) < 10:
        print(f"[deconv] square block too small ({len(block)}); not supported by data.")
    ri = [list(regs).index(g) for g in block]; gi = [list(genes).index(g) for g in block]
    print(f"[deconv] square block size = {len(block)}")
    edges = []
    for c, cond in enumerate(COND_ORDER):
        L = np.nan_to_num(tensor[np.ix_(ri, gi, [c])][:, :, 0], nan=0.0)
        A, cn = deconvolve(L, args.ridge)
        for i, s in enumerate(block):
            for j, t in enumerate(block):
                if i != j and A[i, j] != 0:
                    edges.append(dict(source=s, target=t, weight=float(A[i, j]),
                                      block_condition=cond, cond_number=cn))
    pd.DataFrame(edges).to_csv(TAB / "operator_deconvolution_edges.csv", index=False)
    sub = []
    for c in range(3):
        sets = {g: set(np.where(np.abs(np.nan_to_num(tensor[i, :, c])) > args.sig_thr)[0])
                for i, g in enumerate(regs)}
        for i in block:
            for j in block:
                if i != j and len(sets[j]):
                    sub.append(dict(source=i, target=j, condition=COND_ORDER[c],
                                    subsumption_frac=len(sets[i] & sets[j]) / len(sets[j])))
    pd.DataFrame(sub).to_csv(TAB / "operator_subsumption.csv", index=False)
    print(f"[deconv] {len(edges)} edges; condition numbers per condition reported.")


if __name__ == "__main__":
    main()
