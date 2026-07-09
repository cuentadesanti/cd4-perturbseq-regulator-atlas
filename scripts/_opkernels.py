#!/usr/bin/env python3
"""Pure numerical kernels for the empirical regulatory operator analysis.

No I/O, no argparse, no plotting — unit-tested in tests/test_opkernels.py.
Mirrors the shared-helper convention of scripts/_figstyle.py.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import spearmanr


def assemble_tensor(matrix, obs, gene_idx, cond_order):
    """Stack per-(regulator, condition) rows of `matrix` into a dense-with-mask tensor.

    obs columns: target_contrast_gene_name, culture_condition, n_cells_target,
    ontarget_significant_bool, row. `matrix[row, gene_idx]` supplies each cell;
    unobserved cells are NaN in `tensor`, False in `mask`.
    """
    cond_pos = {c: k for k, c in enumerate(cond_order)}
    sig = obs[obs["ontarget_significant_bool"]]
    regulators = sorted(sig["target_contrast_gene_name"].unique())
    reg_pos = {g: i for i, g in enumerate(regulators)}
    R, G, C = len(regulators), len(gene_idx), len(cond_order)
    tensor = np.full((R, G, C), np.nan, dtype=np.float32)
    mask = np.zeros((R, G, C), dtype=bool)
    n_cells = np.full((R, C), np.nan, dtype=np.float32)
    gene_idx = np.asarray(gene_idx)
    for _, rec in sig.iterrows():
        c = rec["culture_condition"]
        if c not in cond_pos:
            continue
        i, k = reg_pos[rec["target_contrast_gene_name"]], cond_pos[c]
        tensor[i, :, k] = matrix[int(rec["row"]), gene_idx]
        mask[i, :, k] = True
        n_cells[i, k] = rec["n_cells_target"]
    return tensor, mask, np.array(regulators, dtype=object), n_cells


def rms_normalize_conditions(tensor, mask):
    """Scale each condition slab to unit RMS over OBSERVED entries.

    RMS (not raw Frobenius) so the control is a pure per-entry scale control and
    not also a sqrt(coverage) re-weight when conditions have different #observed.
    """
    out = tensor.copy()
    C = tensor.shape[2]
    scales = np.zeros(C, dtype=np.float64)
    for c in range(C):
        vals = tensor[:, :, c][mask[:, :, c]].astype(np.float64)
        rms = float(np.sqrt(np.mean(vals ** 2))) if vals.size else 0.0
        scales[c] = rms
        if rms > 0:
            out[:, :, c] = out[:, :, c] / rms
    return out, scales


def spearman_power(tensor, mask, n_cells):
    norms, ncell = [], []
    R, _, C = tensor.shape
    for i in range(R):
        for c in range(C):
            if mask[i, :, c].any():
                norms.append(float(np.linalg.norm(tensor[i, :, c][mask[i, :, c]])))
                ncell.append(float(n_cells[i, c]))
    if len(norms) < 3:
        return float("nan")
    return float(spearmanr(norms, ncell).statistic)
