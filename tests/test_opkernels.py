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
