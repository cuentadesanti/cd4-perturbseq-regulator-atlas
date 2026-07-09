#!/usr/bin/env python3
"""Step 2 — CP of the z-score operator tensor: regulator (x) gene (x) condition.

Gating ('TCR gated, chromatin constitutive') is the SHAPE of c_k, but only after
the RMS condition-scale control, and only when a bootstrap CI on c_k's range
excludes 'flat'. Precision is handled by the z-score representation (Step 0), so
there is NO fake weight tensor here.

    python scripts/decompose_operator_cp.py --rank auto --scale-control rms
    python scripts/decompose_operator_cp.py --rank 4 --scale-control none   # counterfactual

Outputs: docs/tables/operator_cp_{factors,stability,cosine,enrichment}.csv,
         docs/figures/33_*, docs/figures/34_*
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import _opkernels as op
from decompose_operator_svd import load_genesets

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"; FIG = ROOT / "docs" / "figures"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]
FLAT = 0.15   # normalized-range threshold for "flat"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank", default="auto")
    ap.add_argument("--scale-control", choices=["rms", "none"], default="rms")
    ap.add_argument("--max-rank", type=int, default=8)
    ap.add_argument("--stab-threshold", type=float, default=0.7)
    ap.add_argument("--stab-subsample", type=int, default=400)
    ap.add_argument("--boot-n", type=int, default=100)
    ap.add_argument("--sensitivity", action="store_true")
    args = ap.parse_args()

    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    tensor, mask = d["tensor"].astype(np.float64), d["mask"]
    genes, regs = d["genes"].astype(str), d["regulators"].astype(str)
    n_cells = d["n_cells"].astype(np.float64)
    print(f"[confound] standing meter = {op.spearman_power(tensor, mask, n_cells):.3f} (expect ~0)")

    tensor_fit = op.rms_normalize_conditions(tensor, mask)[0] if args.scale_control == "rms" else tensor
    if args.scale_control == "rms":
        print("[scale-control] RMS-per-observed-entry normalization applied")

    stab_rows = []
    for r in range(2, args.max_rank + 1):
        s = op.split_half_stability(tensor_fit, mask, r, n_splits=6, random_state=0,
                                    subsample=args.stab_subsample)
        stab_rows.append(dict(rank=r, mean_matched_cosine=s))
        print(f"[stability] rank={r} cos={s:.3f}")
    stab = pd.DataFrame(stab_rows)
    if args.rank == "auto":
        ok = stab[stab["mean_matched_cosine"] > args.stab_threshold]
        rank = int(ok["rank"].max()) if len(ok) else int(stab.loc[stab["mean_matched_cosine"].idxmax(), "rank"])
        print(f"[rank] auto={rank}" + ("" if len(ok) else " (argmax; none cleared threshold)"))
    else:
        rank = int(args.rank)

    lam, factors = op.fix_cp_gauge(*op.cp_fit_masked(tensor_fit, mask, rank, n_init=10, random_state=0))
    A, B, C = factors
    degen = op.cp_degeneracy(factors)
    cosmat = op.gene_mode_cosine(factors)
    reg_ncell = np.nanmean(n_cells, axis=1)

    # bootstrap CI on the condition profiles
    _, boot = op.bootstrap_cp_conditions(tensor_fit, mask, rank, n_boot=args.boot_n,
                                          random_state=0, subsample=args.stab_subsample)
    rng_k = boot.max(axis=1) - boot.min(axis=1)                 # (n_boot, rank)
    scale_k = np.abs(boot).max(axis=1) + 1e-9
    norm_range = rng_k / scale_k
    range_lo95 = np.nanpercentile(norm_range, 2.5, axis=0)      # (rank,)

    genesets, bg = load_genesets(), list(genes)
    frows, erows = [], []
    off = cosmat.copy(); np.fill_diagonal(off, 0.0)
    for k in range(rank):
        rho = float(spearmanr(A[:, k], reg_ncell).statistic)
        gsel = np.argsort(-np.abs(B[:, k]))[:50]
        enr = op.hypergeometric_enrichment(list(genes[gsel]), genesets, bg)
        enr.insert(0, "factor", k + 1); erows.append(enr)
        label = enr.iloc[0]["set_name"] if (len(enr) and enr.iloc[0]["fdr"] < 0.05) else "unlabeled"
        peak = COND_ORDER[int(np.argmax(C[:, k]))]
        gated_ci = bool(range_lo95[k] > FLAT)
        frows.append(dict(
            factor=k + 1, lambda_=float(lam[k]),
            cond_Rest=float(C[0, k]), cond_Stim8hr=float(C[1, k]), cond_Stim48hr=float(C[2, k]),
            gating_shape=(f"gated(peak={peak})" if gated_ci else "constitutive(flat)"),
            gated_ci=gated_ci, range_lo95=float(range_lo95[k]),
            power_rho=rho, power_confounded=bool(abs(rho) > 0.3),
            degeneracy=float(degen[k]), max_cofactor_cosine=float(off[k].max()),
            top_regulators=";".join(map(str, regs[np.argsort(-np.abs(A[:, k]))[:8]])),
            top_genes=";".join(map(str, genes[gsel[:10]])), program_label=label))

    TAB.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(frows).to_csv(TAB / "operator_cp_factors.csv", index=False)
    stab.to_csv(TAB / "operator_cp_stability.csv", index=False)
    pd.DataFrame(cosmat, columns=[f"f{k+1}" for k in range(rank)]).to_csv(
        TAB / "operator_cp_cosine.csv", index=False)
    pd.concat(erows, ignore_index=True).to_csv(TAB / "operator_cp_enrichment.csv", index=False)

    FIG.mkdir(exist_ok=True, parents=True)
    f1, a1 = plt.subplots(figsize=(6, 4))
    a1.plot(stab["rank"], stab["mean_matched_cosine"], "o-"); a1.axhline(args.stab_threshold, ls="--", c="grey")
    a1.set_xlabel("CP rank"); a1.set_ylabel("split-half matched cosine")
    a1.set_title("CP rank by split-half stability"); f1.tight_layout()
    f1.savefig(FIG / "33_operator_cp_stability.png", dpi=150)

    f2, a2 = plt.subplots(figsize=(7, 4))
    for k in range(rank):
        lo = np.nanpercentile(boot[:, :, k], 2.5, axis=0); hi = np.nanpercentile(boot[:, :, k], 97.5, axis=0)
        a2.plot(COND_ORDER, C[:, k], "o-", label=f"f{k+1}: {frows[k]['program_label']} "
                f"({'gated' if frows[k]['gated_ci'] else 'flat'})")
        a2.fill_between(COND_ORDER, lo, hi, alpha=0.15)
    a2.set_ylabel("condition modulation c_k (RMS-controlled)"); a2.legend(fontsize=7)
    a2.set_title("CP condition factors with bootstrap 95% CI"); f2.tight_layout()
    f2.savefig(FIG / "34_operator_cp_condition_factors.png", dpi=150)

    clean = [r for r in frows if not r["power_confounded"] and r["degeneracy"] < 0.9
             and r["max_cofactor_cosine"] < 0.7]
    gated = [r["factor"] for r in clean if r["gated_ci"]]
    flat = [r["factor"] for r in clean if not r["gated_ci"]]
    print(f"[result] rank={rank}; clean-gated={gated}; clean-constitutive={flat}")
    if args.sensitivity:
        print("[sensitivity] rebuild Step 0 with an alt panel/gene rule, refit, and compare "
              "gating_shape per matched factor; flips => selection artifact.")


if __name__ == "__main__":
    main()
