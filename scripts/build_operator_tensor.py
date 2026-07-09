#!/usr/bin/env python3
"""Step 0 — build the regulator x gene x condition Z-SCORE tensor (expanded panel).

Panel = original 200-regulator fingerprint panel UNION top ~600 NEW regulators
(KD-significant in all 3 conditions) by ontarget_effect_size. Representation is
the pooled remote layers/zscore (precision-decoupled), fetched once and cached.
Three fail-closed assertions guarantee the representation is pooled z-score before
anything is cached.

    python scripts/build_operator_tensor.py --n-total 800 --top-genes 2000

Outputs: data/cache/operator_tensor.npz, docs/tables/operator_tensor_summary.json
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import _opkernels as op
from analyze_fingerprints import read_matrix, build_panel

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]
TCR_ANCHORS = ["ZAP70", "LCK", "LAT", "CD3D", "CD3E", "CD3G", "CD247", "FYN",
               "ITK", "PLCG1", "LCP2", "VAV1", "PRKCQ", "CARD11", "BCL10", "MALT1"]


def load_obs():
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv").reset_index(drop=True)
    obs["row"] = np.arange(len(obs))
    obs["ontarget_significant_bool"] = (
        obs["ontarget_significant"].astype(str).str.strip().eq("True"))
    return obs


def select_expanded_panel(obs, n_total):
    """original-200 fingerprint panel UNION top-new by effect, all sig in 3 conds."""
    sig = obs[obs["ontarget_significant_bool"]]
    conds = sig.groupby("target_contrast_gene_name")["culture_condition"].apply(
        lambda s: set(s) & set(COND_ORDER))
    all3 = set(conds[conds.apply(lambda s: set(COND_ORDER).issubset(s))].index)
    original = set(build_panel(obs, 200)["gene"]) & all3
    eff = (sig[sig["target_contrast_gene_name"].isin(all3)]
           .groupby("target_contrast_gene_name")["ontarget_effect_size"].max()
           .sort_values(ascending=False))
    keep = list(original)
    for g in eff.index:
        if len(keep) >= n_total:
            break
        if g not in original:
            keep.append(g)
    return set(keep), original


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-total", type=int, default=800)
    ap.add_argument("--top-genes", type=int, default=2000)
    args = ap.parse_args()

    obs = load_obs()
    keep, original = select_expanded_panel(obs, args.n_total)

    # rows for the panel (all sig rows of kept regulators in the 3 conditions)
    panel = obs[obs["ontarget_significant_bool"]
                & obs["target_contrast_gene_name"].isin(keep)
                & obs["culture_condition"].isin(COND_ORDER)].copy()
    panel = panel.reset_index(drop=True)
    Z = read_matrix(panel[["row"]], "zscore")        # (n_panel, 10282), row-aligned
    Z = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)
    rownorm_cv = float(np.linalg.norm(Z, axis=1).std() / (np.linalg.norm(Z, axis=1).mean() + 1e-12))

    # gene axis: top variance in Z-SCORE space (power-decoupled; note this in writeup)
    v = Z.var(axis=0)
    gene_idx = np.argsort(v)[::-1][:args.top_genes]
    var = pd.read_csv(CACHE / "fingerprint_var.csv")
    genes = var.loc[gene_idx, "gene_name"].to_numpy()

    # assemble with LOCAL row index (Z is row-aligned to `panel`)
    panel["row"] = np.arange(len(panel))
    tensor, mask, regulators, n_cells = op.assemble_tensor(Z, panel, gene_idx, COND_ORDER)
    in_orig = np.array([g in original for g in regulators], dtype=bool)

    # ---- fail-closed representation assertions ----
    rho = op.spearman_power(tensor, mask, n_cells)
    anchor_cvs = []
    for i, g in enumerate(regulators):
        if g in TCR_ANCHORS:
            norms = [np.linalg.norm(tensor[i, :, c][mask[i, :, c]])
                     for c in range(3) if mask[i, :, c].any()]
            if len(norms) == 3 and np.mean(norms) > 0:
                anchor_cvs.append(np.std(norms) / np.mean(norms))
    anchor_cv = float(np.median(anchor_cvs)) if anchor_cvs else float("nan")

    errs = []
    if not (rownorm_cv > 0.30):
        errs.append(f"row-norm CV {rownorm_cv:.3f} <= 0.30 (representation may be within-condition z)")
    if not (anchor_cv > 0.10):
        errs.append(f"TCR-anchor cross-condition CV {anchor_cv:.3f} <= 0.10 (gating would be dead)")
    if not (abs(rho) < 0.15):
        errs.append(f"confound meter |rho| {abs(rho):.3f} >= 0.15 (power re-entered)")
    if errs:
        raise SystemExit("[REFUSE-TO-CACHE] representation assertions failed:\n  - "
                         + "\n  - ".join(errs))

    CACHE.mkdir(exist_ok=True, parents=True)
    np.savez(CACHE / "operator_tensor.npz",
             tensor=tensor, mask=mask, regulators=regulators.astype(str),
             genes=genes.astype(str), conditions=np.array(COND_ORDER),
             n_cells=n_cells, in_original_panel=in_orig)
    summary = dict(n_regulators=int(len(regulators)),
                   n_new_regulators=int((~in_orig).sum()),
                   n_genes=int(len(genes)), conditions=COND_ORDER,
                   n_cells_confound_rho=rho, rownorm_cv=rownorm_cv,
                   anchor_cross_cond_cv=anchor_cv,
                   observed_cells=int(mask.any(axis=1).sum()),
                   representation="pooled_zscore", passed_assertions=True)
    TAB.mkdir(exist_ok=True, parents=True)
    (TAB / "operator_tensor_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
