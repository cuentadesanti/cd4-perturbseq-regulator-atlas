import numpy as np
import pandas as pd
import pytest
import _opkernels as op


def _toy_obs():
    rows, r = [], 0
    plan = {"A": ["Rest", "Stim8hr", "Stim48hr"],
            "B": ["Rest", "Stim8hr", "Stim48hr"],
            "C": ["Rest", "Stim8hr"]}   # C missing Stim48hr
    for g, conds in plan.items():
        for c in conds:
            rows.append(dict(target_contrast_gene_name=g, culture_condition=c,
                             n_cells_target=100.0 + r, ontarget_significant_bool=True,
                             row=r)); r += 1
    return pd.DataFrame(rows)


def test_assemble_tensor_shapes_and_mask():
    obs = _toy_obs()
    G_all = 5
    matrix = np.arange(len(obs) * G_all, dtype=np.float32).reshape(len(obs), G_all)
    tensor, mask, regs, n_cells = op.assemble_tensor(
        matrix, obs, np.array([0, 2, 4]), ["Rest", "Stim8hr", "Stim48hr"])
    assert tensor.shape == (3, 3, 3) and mask.shape == (3, 3, 3)
    assert list(regs) == ["A", "B", "C"]
    ci = list(regs).index("C")
    assert mask[ci, :, 2].sum() == 0 and np.isnan(tensor[ci, :, 2]).all()
    assert mask[list(regs).index("A"), :, 2].sum() == 3


def test_rms_normalize_conditions_unit_rms_over_observed():
    rng = np.random.default_rng(0)
    tensor = rng.normal(size=(4, 6, 3)).astype(np.float32)
    mask = np.ones_like(tensor, dtype=bool)
    mask[0, :, 2] = False                     # ragged: fewer observed in slab 2
    tn, scales = op.rms_normalize_conditions(tensor, mask)
    for c in range(3):
        vals = tn[:, :, c][mask[:, :, c]]
        assert np.isclose(np.sqrt(np.mean(vals ** 2)), 1.0, atol=1e-5)
    assert scales.shape == (3,)


def test_spearman_power_detects_planted_confound():
    rng = np.random.default_rng(1)
    R, G = 20, 8
    n_cells = rng.uniform(50, 5000, size=(R, 3))
    tensor = np.zeros((R, G, 3), dtype=np.float32)
    for c in range(3):
        for i in range(R):
            tensor[i, :, c] = n_cells[i, c] / 1000.0
    mask = np.ones((R, G, 3), dtype=bool)
    assert op.spearman_power(tensor, mask, n_cells) > 0.9


def _varimax_crit(X):
    # varimax objective: sum over factors of the variance of squared loadings
    return float(((X ** 2).var(axis=0)).sum())


def test_varimax_recovers_structure_and_is_varimax():
    # asymmetric simple structure (a balanced/symmetric toy is a saddle for varimax),
    # then mixed by a random orthogonal rotation for varimax to undo.
    rng = np.random.default_rng(0)
    L = np.zeros((9, 3))
    L[0:4, 0] = [0.9, 0.8, 0.7, 0.6]
    L[4:7, 1] = [0.85, 0.75, 0.65]
    L[7:9, 2] = [0.9, 0.5]
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    M = L @ Q
    rot, R = op.varimax(M)
    # (1) recovers simple structure: nearly every variable dominant on one factor
    dom = np.abs(rot).max(axis=1) / (np.abs(rot).sum(axis=1) + 1e-9)
    assert (dom > 0.85).mean() >= 8 / 9
    # (2) it is VARIMAX (not quartimax) and did NOT early-stop: the varimax criterion
    #     strictly increases vs the mixed input
    assert _varimax_crit(rot) > _varimax_crit(M) + 1e-6
    # (3) the objective is monotone non-decreasing across iterations (probe via max_iter)
    crits = [_varimax_crit(op.varimax(M, max_iter=t)[0]) for t in range(1, 8)]
    assert all(crits[i + 1] >= crits[i] - 1e-9 for i in range(len(crits) - 1))
    # (4) orthonormal rotation
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-6)


def test_hypergeometric_enrichment_flags_planted_set():
    gene_list = ["ISG15", "MX1", "OAS1", "IFIT1", "STAT1"]
    gene_sets = {"IFN": ["ISG15", "MX1", "OAS1", "IFIT1", "STAT1", "IRF7", "IFI6"],
                 "RANDOM": ["AAA", "BBB", "CCC", "DDD"]}
    background = gene_list + [f"BG{i}" for i in range(500)] + gene_sets["IFN"]
    res = op.hypergeometric_enrichment(gene_list, gene_sets, background)
    top = res.sort_values("pvalue").iloc[0]
    assert top["set_name"] == "IFN" and top["fdr"] < 0.05
