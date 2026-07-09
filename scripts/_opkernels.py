#!/usr/bin/env python3
"""Pure numerical kernels for the empirical regulatory operator analysis.

No I/O, no argparse, no plotting — unit-tested in tests/test_opkernels.py.
Mirrors the shared-helper convention of scripts/_figstyle.py.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom
from statsmodels.stats.multitest import multipletests


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


def varimax(loadings, gamma=1.0, max_iter=100, tol=1e-6):
    """Kaiser-normalized VARIMAX rotation of a (G, k) loading matrix.

    gamma=1.0 is varimax (maximizes the per-column variance of squared loadings ->
    simple structure ACROSS factors). Do NOT set gamma=0 (that is quartimax, a
    different criterion that tends to a general factor). Convergence uses the
    canonical Kaiser ratio criterion on the singular-value sum (the objective is
    monotone non-decreasing under this algorithm); a rotation-change tolerance can
    fire on the first step and exit at a non-simple-structure solution.
    """
    L = np.asarray(loadings, dtype=np.float64)
    G, k = L.shape
    if k < 2:
        return L.copy(), np.eye(k)
    h = np.sqrt((L ** 2).sum(axis=1, keepdims=True)); h[h == 0] = 1.0
    Ln = L / h
    Rm = np.eye(k)
    d_old = 0.0
    for _ in range(max_iter):
        Lr = Ln @ Rm
        u, s, vt = np.linalg.svd(
            Ln.T @ (Lr ** 3 - (gamma / G) * Lr @ np.diag((Lr ** 2).sum(axis=0))))
        Rm = u @ vt
        d = s.sum()
        if d_old != 0 and d / d_old < 1 + tol:
            break
        d_old = d
    return (Ln @ Rm) * h, Rm


def hypergeometric_enrichment(gene_list, gene_sets, background):
    bg = set(background); drawn = set(gene_list) & bg
    M, n_drawn = len(bg), len(drawn)
    rows = []
    for name, members in gene_sets.items():
        setg = set(members) & bg; K = len(setg)
        overlap = drawn & setg; x = len(overlap)
        p = hypergeom.sf(x - 1, M, K, n_drawn) if K > 0 and n_drawn > 0 else 1.0
        rows.append(dict(set_name=name, n_overlap=x, set_size=K, n_drawn=n_drawn,
                         background_size=M, pvalue=float(p),
                         overlap_genes=",".join(sorted(overlap))))
    df = pd.DataFrame(rows)
    if len(df):
        df["fdr"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    return df.sort_values("pvalue").reset_index(drop=True)
