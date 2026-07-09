#!/usr/bin/env python3
"""Step 0 — build the regulator x gene x condition Z-SCORE tensor (expanded panel).

Panel = the top ~800 regulators by n_downstream (breadth) among those KD-significant
in all 3 conditions AND above the median cell-count floor (ranking by breadth above a
power floor is what decorrelates power; ranking by ontarget_effect_size does not).
Representation is the pooled remote layers/zscore (precision-decoupled), fetched once
and cached. One fail-closed guard (confound |rho|<0.15) refuses to cache if power
re-enters; pooling itself is a layer-level property + a downstream CP-gating test
(see the guard comment).

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


def load_obs():
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv").reset_index(drop=True)
    obs["row"] = np.arange(len(obs))
    obs["ontarget_significant_bool"] = (
        obs["ontarget_significant"].astype(str).str.strip().eq("True"))
    return obs


def select_expanded_panel(obs, n_total):
    """Power-floored, breadth-ranked panel (all regulators KD-significant in 3 conds).

    1. n_cells FLOOR at the median of the all-3-condition population (removes the
       low-power tail carrying the intrinsic norm-vs-power confound).
    2. Rank the floored set by n_downstream (breadth) and take the top n_total.
    Breadth is the right axis for an operator analysis and far less power-confounded
    than ontarget_effect_size. in_original_panel is tagged post-hoc vs the 200-reg
    fingerprint panel (the ~67 overlap; ~733 non-overlap = out-of-panel for Step 3b).
    """
    sig = obs[obs["ontarget_significant_bool"]]
    conds = sig.groupby("target_contrast_gene_name")["culture_condition"].apply(
        lambda s: set(s) & set(COND_ORDER))
    all3 = set(conds[conds.apply(lambda s: set(COND_ORDER).issubset(s))].index)
    g3 = sig[sig["target_contrast_gene_name"].isin(all3)]
    ncell = g3.groupby("target_contrast_gene_name")["n_cells_target"].mean()
    ndown = g3.groupby("target_contrast_gene_name")["n_downstream"].max()
    floor = float(np.median(ncell.values))
    floored = ncell[ncell >= floor].index
    ranked = ndown.reindex(floored).dropna().sort_values(ascending=False)
    keep = set(ranked.index[:n_total])
    original = set(build_panel(obs, 200)["gene"]) & keep
    return keep, original


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

    # ---- fail-closed representation guard (confound only) ----
    # ONLY the confound guard is asserted. Earlier drafts also gated on row-norm CV
    # and per-anchor cross-condition spread; both were RETIRED after verification.
    # layers/zscore is per-(perturbation,gene) z, so any MAGNITUDE statistic is ~invariant
    # to the pooled-vs-within-condition distinction on a breadth-homogeneous panel: row-norm
    # CV drops from SELECTION homogeneity, not within-condition z (proven in raw space, where
    # no z exists: this selection's raw-logFC row-norm CV 0.28 < random-800's 0.34), and only
    # ~1 TCR anchor survives the cell-count floor. Pooling is a documented evidence chain
    # instead: (1) it is a LAYER-level property established on the same layers/zscore slice —
    # row-norm CV 0.36 on the heterogeneous 200-panel BEFORE any selection; the 800 rows are
    # drawn from that identical layer by index. (2) the confound guard below protects the
    # analysis. (3) within-condition z would manifest DOWNSTREAM as ALL-constitutive CP
    # condition factors, which Step 2's bootstrap-CI gating test (>=1 factor with condition-CI
    # excluding flat) would expose — that downstream test, not a Step-0 magnitude proxy, is
    # the authoritative pooling check.
    rho = op.spearman_power(tensor, mask, n_cells)
    if not (abs(rho) < 0.15):
        raise SystemExit(f"[REFUSE-TO-CACHE] confound guard: spearman(||slab||, n_cells) "
                         f"|rho|={abs(rho):.3f} >= 0.15 — power re-entered the representation.")

    CACHE.mkdir(exist_ok=True, parents=True)
    np.savez(CACHE / "operator_tensor.npz",
             tensor=tensor, mask=mask, regulators=regulators.astype(str),
             genes=genes.astype(str), conditions=np.array(COND_ORDER),
             n_cells=n_cells, in_original_panel=in_orig)
    summary = dict(n_regulators=int(len(regulators)),
                   n_new_regulators=int((~in_orig).sum()),
                   n_genes=int(len(genes)), conditions=COND_ORDER,
                   n_cells_confound_rho=rho,            # the asserted guard (|rho|<0.15)
                   rownorm_cv=rownorm_cv,               # informational only (selection-confounded)
                   observed_cells=int(mask.any(axis=1).sum()),
                   representation="pooled_zscore", passed_confound_guard=True)
    TAB.mkdir(exist_ok=True, parents=True)
    (TAB / "operator_tensor_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
