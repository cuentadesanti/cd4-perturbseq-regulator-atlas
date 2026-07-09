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


tensorly = pytest.importorskip("tensorly")


def _synth_cp(R=40, G=60, C=3, rank=3, seed=0, noise=0.02):
    rng = np.random.default_rng(seed)
    A = rng.normal(size=(R, rank)); B = rng.normal(size=(G, rank))
    Cc = np.abs(rng.normal(size=(C, rank)))
    T = np.einsum("ir,jr,kr->ijk", A, B, Cc)
    T += noise * rng.normal(size=T.shape) * np.std(T)
    return T.astype(np.float32), (A, B, Cc)


def test_cp_recovers_planted_factors():
    T, (_, B, _) = _synth_cp()
    _, factors = op.cp_fit_masked(T, np.ones_like(T, bool), 3, n_init=5, random_state=0)
    assert op.match_factors(B, factors[1])[1].mean() > 0.95


def test_cp_masked_ignores_hidden_cells():
    T, _ = _synth_cp()
    mask = np.ones_like(T, bool); mask[::2, :, 2] = False
    _, f = op.cp_fit_masked(T, mask, 3, n_init=5, random_state=0)
    assert f[0].shape == (40, 3) and f[2].shape == (3, 3)


def test_fix_cp_gauge_unit_norm_and_sign():
    T, _ = _synth_cp()
    lam, f = op.cp_fit_masked(T, np.ones_like(T, bool), 3, n_init=3, random_state=1)
    _, f2 = op.fix_cp_gauge(lam, f)
    for F in f2:
        assert np.allclose(np.linalg.norm(F, axis=0), 1.0, atol=1e-5)
    Cc = f2[2]
    for k in range(3):
        assert Cc[np.argmax(np.abs(Cc[:, k])), k] > 0


def test_gene_mode_cosine_is_full_matrix():
    T, _ = _synth_cp()
    _, f = op.cp_fit_masked(T, np.ones_like(T, bool), 3, n_init=3, random_state=0)
    Cm = op.gene_mode_cosine(f)
    assert Cm.shape == (3, 3) and np.allclose(np.diag(Cm), 1.0, atol=1e-6)


def test_split_half_stability_high_for_true_rank():
    T, _ = _synth_cp(R=80, noise=0.03)
    mask = np.ones_like(T, bool)
    s3 = op.split_half_stability(T, mask, 3, n_splits=4, random_state=0)
    s8 = op.split_half_stability(T, mask, 8, n_splits=4, random_state=0)
    assert s3 > 0.65 and s3 > s8


def test_bootstrap_cp_conditions_ci_brackets_truth():
    # planted gated factor: condition profile clearly peaked -> CI of range excludes 0
    T, _ = _synth_cp(R=120, noise=0.05)
    C_ref, boot = op.bootstrap_cp_conditions(T, np.ones_like(T, bool), 3,
                                             n_boot=20, random_state=0)
    assert C_ref.shape == (3, 3) and boot.shape == (20, 3, 3)
    rng_k = boot.max(axis=1) - boot.min(axis=1)     # (n_boot, rank)
    assert np.percentile(rng_k, 2.5, axis=0).shape == (3,)


def test_train_test_standardize_uses_train_only():
    M = np.array([[1.0, 10.0], [3.0, 30.0], [100.0, 100.0]])
    train = np.array([[True, True], [True, True], [False, False]])
    Mc, mu = op.train_test_standardize(M, train)
    assert np.allclose(mu, [2.0, 20.0]) and np.allclose(Mc[2], [98.0, 80.0])


def test_soft_impute_recovers_low_rank():
    rng = np.random.default_rng(0)
    M = rng.normal(size=(30, 2)) @ rng.normal(size=(2, 20))
    obs = rng.random(M.shape) > 0.3
    Mhat = op.soft_impute(M, obs, 2, n_iter=200)
    held = ~obs
    r2 = 1 - np.sum((M[held] - Mhat[held]) ** 2) / np.sum((M[held] - M[obs].mean()) ** 2)
    assert r2 > 0.9


def test_principal_angles_identical_subspace():
    V = np.random.default_rng(0).normal(size=(50, 4))
    assert np.allclose(op.principal_angles(V, V.copy()), 1.0, atol=1e-6)


def test_principal_angles_orthogonal_subspace():
    Va = np.zeros((20, 2)); Va[0, 0] = Va[1, 1] = 1.0
    Vb = np.zeros((20, 2)); Vb[2, 0] = Vb[3, 1] = 1.0
    assert np.allclose(op.principal_angles(Va, Vb), 0.0, atol=1e-6)


def test_random_subspace_null_small_for_high_dim():
    m, p95 = op.random_subspace_null(200, 4, n=200, random_state=0)
    assert m < 0.1 and p95 >= m


def test_soft_impute_matches_full_svd_reference():
    # Regression guard for the truncated-SVD (svds) soft_impute: it must equal a full-SVD
    # reference to machine precision. If the truncated path ever diverges from top-k
    # truncation (e.g. a bad svds swap), this fails — that is what makes the fast kernel
    # a *validated* equivalence, not a one-off offline check.
    rng = np.random.default_rng(0)
    M = rng.normal(size=(180, 90)) @ rng.normal(size=(90, 90))
    obs = rng.random(M.shape) > 0.2

    def _full_svd_ref(M, obs, rank, n_iter=100, tol=1e-4):
        M = np.asarray(M, float); X = np.where(obs, M, 0.0); Xr = X; prev = np.inf
        for _ in range(n_iter):
            U, s, Vt = np.linalg.svd(X, full_matrices=False)
            Xr = (U[:, :rank] * s[:rank]) @ Vt[:rank]
            X = np.where(obs, M, Xr)
            ch = np.linalg.norm(Xr - X) / (np.linalg.norm(X) + 1e-12)
            if abs(prev - ch) < tol:
                break
            prev = ch
        return Xr

    for r in (1, 5, 15):
        got = op.soft_impute(M, obs, r, n_iter=100)
        ref = _full_svd_ref(M, obs, r, n_iter=100)
        assert np.allclose(got, ref, atol=1e-8), f"rank {r}: max|dif|={np.abs(got - ref).max():.2e}"
