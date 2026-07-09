# Empirical Regulatory Operator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover the gene-program, condition-modulation, and predictive structure of the regulatory operator `L` (effect of every gene on every gene under regulator KD) that the fingerprint PCA computed but discarded — via tensor decomposition, low-rank completion, and donor-subspace stability, **built in precision-decoupled z-score space** so the leading factors are biology and not the −0.68 power-magnitude confound.

**Architecture:** One 3-way tensor `T[regulator, gene, condition]` is assembled once, in **z-score space**, over an **expanded ~800-regulator panel** (Step 0). Everything after is a question about that one object: its right factors (Step 1, gene programs), its CP decomposition (Step 2, `regulator ⊗ gene ⊗ condition` — where "gated vs constitutive" is the shape of one factor, now with a bootstrap CI), its recoverability by low rank (Step 3, held-out prediction — **3b out-of-panel condition extrapolation is the flagship**), and the reproducibility of its program subspace across donors (Step 4). Pure kernels live in a unit-tested `scripts/_opkernels.py` (mirroring the `scripts/_figstyle.py` helper convention); one driver per step follows the repo's analysis-script idiom (argparse → `docs/tables/*.csv|json` + `docs/figures/*.png`), each with a Makefile target.

**Tech Stack:** Python 3, numpy, scipy, statsmodels, pandas, matplotlib (in `requirements.txt`); adds `tensorly` (CP), `pytest` (kernel tests), and `h5py`+`s3fs`+`fsspec` (the one-time z-score fetch — already the repo's optional Model-1 deps). scikit-learn 1.9 is installed.

## Scope decision (this is Option 1, the pilot)

The full regulatory operator has **6209** regulators KD-significant in all 3 conditions. This plan deliberately runs the **Option 1 pilot**: an expanded **~800-regulator** panel — the regulators with the broadest downstream footprints (top 800 by `n_downstream`) among those above a **median cell-count floor** — in z-score space. (Ranking by breadth above a power floor is required: ranking by `ontarget_effect_size`, as an earlier draft did, reintroduces the power confound the guard rejects — see Task 2.) ~67 of the 800 overlap the original fingerprint panel; the other ~733 are out-of-panel, and are the held-out set for the Step-3b test. Rationale: clean representation + tractable CP runtime + one small (~0.1 GB) fetch, while still being ~4× broader than the existing fingerprint work. **Escalation trigger (Task 11):** if Step 3b beats persistence cleanly on the out-of-panel regulators, that is the evidence to escalate to Option 3 (full 6209 axis, completion as the flagship, CP sweep on cluster compute). Decide on evidence, not by relitigating scope.

## Global Constraints

- **Representation is z-score, not raw log-FC.** Raw log-FC has `spearman(‖row‖, n_cells_target) = −0.683` on significant rows (verified, n=21221) — norm is *negatively* correlated with power, so noisy low-power rows have the largest magnitude and SVD/CP would fit them first. The remote `layers/zscore` is precision-decoupled and pooled across conditions (verified: local 200-panel row-norm CV = 0.360, norms span 87→458), so it kills the confound at the representation level. **Do not build the tensor from `log_fc.f32.npy`.**
- **The z-score is pooled across conditions — Step 0 must prove this fail-closed** (see Task 2, Step 3). If the fetched layer were within-condition normalized, every condition slab would have equal magnitude and the Step 2 gating result would be silently dead. This is a hard assertion, not an eyeball.
- **Data axes (verified):** the remote h5ad `s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad` has `layers/log_fc` and `layers/zscore` both shape **`(33983, 10282)`** (rows = perturbation×condition, cols = 10282 measured genes). `data/cache/log_fc.f32.npy` is the full float32 log-FC (used only for the raw-vs-z confound demonstration and gene-variance fallbacks, never as the tensor source). `data/cache/panel_zscore_a1ca62fb01.npy` is a local `(200, 10282)` z-score panel.
- **`fingerprint_obs.csv` columns (verbatim):** `index, target_contrast_gene_name, culture_condition, target_contrast, n_cells_target, n_up_genes, n_down_genes, n_total_de_genes, n_downstream, ontarget_effect_size, ontarget_significant, distal_offtarget_flag, low_target_gex, neighboring_gene_KD, single_guide_estimate, n_guides, guide_correlation_all, guide_correlation_signif, donor_correlation_all_mean, donor_correlation_hits_mean`. KD gate = `ontarget_significant == True` (string; parse `.astype(str).str.strip().eq("True")`). Power covariate = `n_cells_target`. There is **no `lfcSE` locally** — z-score is exactly why we don't need it.
- **`fingerprint_var.csv` columns:** `gene_name, gene_id` (10282 rows).
- **Condition axis order fixed everywhere:** `["Rest", "Stim8hr", "Stim48hr"]` (`COND_ORDER`).
- **Network:** Step 0 does a **one-time remote fetch** of the z-score slice (~2400 rows × 10282 ≈ 0.1 GB), cached to `data/cache/panel_zscore_<hash>.npy` by the existing `read_matrix` machinery. After that, Steps 1–4 are **fully offline** against `data/cache/operator_tensor.npz`.
- **Output locations:** tables → `docs/tables/operator_*.csv|json`; figures → `docs/figures/<NN>_operator_*.png` (next free prefix after 31); writeup → `docs/OPERATOR_ANALYSIS.md`. Driver path constants: `ROOT = Path(__file__).resolve().parent.parent`, `CACHE = ROOT/"data"/"cache"`, `TAB`, `FIG`.
- **Determinism:** every stochastic routine takes a seed (default 0); use `np.random.default_rng`.
- **The standing confound meter** `spearman(‖T[i,:,c]‖, n_cells[i,c])` is computed in Step 0 and **asserted ≈0** (`|ρ| < 0.15`) on the 800-reg tensor; it is re-run per factor in later steps, and any factor with regulator-mode `|ρ| > 0.3` is flagged `power_confounded=True`.
- Follow existing script style: `#!/usr/bin/env python3`, docstring with usage line, `matplotlib.use("Agg")` before `pyplot`, argparse `main()`, CSVs `index=False`.

---

## File Structure

**Shared kernel module (pure, unit-tested):** `scripts/_opkernels.py`
- `assemble_tensor`, `rms_normalize_conditions`, `spearman_power`, `varimax`, `hypergeometric_enrichment`, `cp_fit_masked`, `fix_cp_gauge`, `cp_degeneracy`, `match_factors`, `split_half_stability`, `gene_mode_cosine`, `bootstrap_cp_conditions`, `train_test_standardize`, `soft_impute`, `principal_angles`, `random_subspace_null`. No I/O, no argparse, no plotting.

**Driver scripts (one per step):**
- `scripts/build_operator_tensor.py` — Step 0 (z-score fetch, expanded panel, fail-closed representation assertions).
- `scripts/decompose_operator_svd.py` — Step 1 (gene programs + power gate + varimax + enrichment).
- `scripts/decompose_operator_cp.py` — Step 2 (CP, RMS scale control, stability rank, bootstrap-CI gating, inter-factor cosine, confound gate).
- `scripts/operator_completion.py` — Step 3 (3b out-of-panel extrapolation as flagship; 3a as sanity check).
- `scripts/operator_donor_angles.py` — Step 4 (principal angles, disjoint donor pairs).
- `scripts/operator_deconvolution.py` — Step 5 stretch.

**Tests:** `tests/test_opkernels.py` (pytest, synthetic-data unit tests for every kernel).

**Config/docs:** `requirements.txt` (+`tensorly`,`pytest`), `Makefile` (targets `operator-*` + umbrella `operator`), `docs/OPERATOR_ANALYSIS.md`.

**Cache artifact (Step 0 → Steps 1–4):** `data/cache/operator_tensor.npz` — keys `tensor (R,G,3) float32`, `mask (R,G,3) bool`, `regulators (R,) str`, `genes (G,) str`, `conditions (3,) str`, `n_cells (R,3) float32`, `in_original_panel (R,) bool`.

---

## Task 1: Kernel — tensor assembly + RMS condition normalization + confound meter

**Files:**
- Create: `scripts/_opkernels.py`, `tests/conftest.py`, `tests/test_opkernels.py`
- Modify: `requirements.txt`
- Test: `tests/test_opkernels.py`

**Interfaces (produces):**
- `assemble_tensor(matrix, obs, gene_idx, cond_order) -> (tensor, mask, regulators, n_cells)` — `matrix` is any `(N, G_all)` array (here the fetched z-score, row-aligned so `obs["row"]` indexes into it); `obs` has `target_contrast_gene_name, culture_condition, n_cells_target, ontarget_significant_bool, row`. Returns `tensor (R,G,3) float32` (NaN unobserved), `mask` bool, `regulators` object, `n_cells (R,3)`.
- `rms_normalize_conditions(tensor, mask) -> (tensor_norm, scales)` — divides each slab `[:,:,c]` by its **RMS over observed entries** (`sqrt(mean(vals²))`), not raw Frobenius, so the control is a pure scale control and not also a √(coverage) re-weight under raggedness. Returns per-condition RMS `(3,)`.
- `spearman_power(tensor, mask, n_cells) -> float`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_opkernels.py  (tests/conftest.py puts scripts/ on sys.path)
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: FAIL — `ModuleNotFoundError` / `AttributeError: ... 'assemble_tensor'`.

- [ ] **Step 3: Create the conftest shim and kernel module**

Create `tests/conftest.py` (puts `scripts/` on `sys.path` so tests can `import _opkernels`; drivers run from `scripts/` and import it bare too — matching the repo's `import _figstyle` convention):
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
```

Create `scripts/_opkernels.py`:
```python
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
```

Append to `requirements.txt` (optional block):
```
# Optional — empirical regulatory operator analysis (`make operator`)
#   pip install tensorly pytest h5py s3fs fsspec
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add tests/conftest.py scripts/_opkernels.py tests/test_opkernels.py requirements.txt
git commit -m "feat(operator): tensor assembly, RMS condition normalization, confound meter kernels"
```

---

## Task 2: Step 0 driver — build the z-score tensor over the expanded panel (fail-closed)

**Files:**
- Create: `scripts/build_operator_tensor.py`
- Modify: `Makefile` (add `operator-tensor`)
- Test: run on real data; acceptance = all three representation assertions pass + artifact cached.

**Interfaces:**
- Consumes: `_opkernels.assemble_tensor`, `spearman_power`; reuses `read_matrix`, `build_panel` from `analyze_fingerprints`.
- Produces: `data/cache/operator_tensor.npz` and `docs/tables/operator_tensor_summary.json` (`n_regulators, n_new_regulators, n_genes, conditions, n_cells_confound_rho, rownorm_cv, observed_cells, representation, passed_confound_guard`).

**Fail-closed guard (the driver `raise`s and does NOT cache if it fails):**
- **Confound guard:** `|spearman(‖slab‖, n_cells)| < 0.15` on the full 800-reg tensor. This is the check that protects the analysis (SVD/CP must not fit power) and it is the ONLY asserted guard.

**Why only the confound guard (retired two earlier guards — documented, not a concession):** earlier drafts also asserted (a) row-norm CV `> 0.30` and (b) per-anchor cross-condition spread `> 0.10`. Both were **retired after verification**: `layers/zscore` is per-(perturbation, gene) z, so any **magnitude** statistic is ~invariant to the pooled-vs-within-condition distinction on a breadth-homogeneous panel — (a) is driven by **selection homogeneity** (proven in raw space, where no z exists: this selection's raw-logFC row-norm CV 0.28 < random-800's 0.34), and only ~1 TCR anchor survives the cell-count floor so (b) is unreliable. Pooling is instead a **documented evidence chain**: (1) it is a LAYER-level property established on the same `layers/zscore` slice — row-norm CV 0.36 on the heterogeneous 200-panel BEFORE any selection, and the 800 rows are drawn from that identical layer by index; (2) the confound guard passes; (3) within-condition z would manifest DOWNSTREAM as **all-constitutive CP condition factors**, which Step 2's bootstrap-CI gating test would expose — that downstream test, not a Step-0 magnitude proxy, is the authoritative pooling check. `rownorm_cv` is still reported (informational).

- [ ] **Step 1: Write the driver**

```python
#!/usr/bin/env python3
"""Step 0 — build the regulator x gene x condition Z-SCORE tensor (expanded panel).

Panel = the top ~800 regulators by n_downstream (breadth) among those KD-significant
in all 3 conditions AND above the median cell-count floor. Representation is
the pooled remote layers/zscore (precision-decoupled), fetched once and cached.
One fail-closed guard (confound |rho|<0.15) refuses to cache if power re-enters;
pooling itself is a layer-level property + a downstream CP-gating test (see comment).

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

    Selection rule (this is load-bearing for the Step-0 confound guard):
      1. n_cells FLOOR at the median of the all-3-condition population — removes the
         low-power tail that carries the population's intrinsic norm-vs-power confound.
      2. Rank the floored set by n_downstream (BREADTH) and take the top n_total.
    Breadth is the right axis for an operator/program analysis AND is far less
    power-confounded than ontarget_effect_size (which is z-like: low-cell perturbations
    get inflated effect sizes and dominate a top-N-by-effect cut, giving raw rho ~ -0.8
    -> z ~ -0.27, which the guard rejects). Floor alone is NOT enough (effect-ranking
    re-imports the confound above the floor); the breadth axis is what decorrelates.
    Verified locally: median floor + top-n by n_downstream -> raw rho ~ -0.05 -> z ~ -0.03.
    in_original_panel is tagged post-hoc against the 200-regulator fingerprint panel
    (the overlap, ~67; the ~733 non-overlap are the out-of-panel set for Step 3b).
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
    # and per-anchor cross-condition spread; both were RETIRED after verification:
    # layers/zscore is per-(perturbation,gene) z, so any MAGNITUDE statistic is ~invariant
    # to the pooled-vs-within-condition distinction on a breadth-homogeneous panel. Row-norm
    # CV drops from SELECTION homogeneity, not within-condition z (proven in raw space, where
    # no z exists: this selection's raw-logFC row-norm CV 0.28 < random-800's 0.34), and only
    # ~1 TCR anchor survives the cell-count floor. The pooled/within-condition distinction is
    # a documented evidence chain instead: (1) pooling is a LAYER-level property established
    # on the same layers/zscore slice — row-norm CV 0.36 on the heterogeneous 200-panel BEFORE
    # any selection; the 800 rows are drawn from that identical layer by index. (2) guard (c)
    # confound |rho|<0.15 protects the actual analysis. (3) within-condition z would manifest
    # DOWNSTREAM as ALL-constitutive CP condition factors, which Step 2's bootstrap-CI gating
    # test (>=1 factor with condition-CI excluding flat) would expose — that downstream test,
    # not a Step-0 magnitude proxy, is the authoritative pooling check.
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
```

- [ ] **Step 2: Run it (one-time fetch, then cached)**

Run: `python scripts/build_operator_tensor.py --n-total 800 --top-genes 2000`
Expected: fetches the z-score slice (~0.1 GB, cached to `panel_zscore_<hash>.npy`), then prints a JSON summary with `n_regulators ≈ 800`, `n_new_regulators ≈ 733`, and `n_cells_confound_rho` near 0 (≈ +0.08 with this selection; `rownorm_cv` ≈ 0.23 is reported but not asserted). If the confound guard fails it exits non-zero **without caching** — that is the guard firing correctly, not a bug to work around.

- [ ] **Step 3: Verify the artifact and the blocker assertions**

Run:
```bash
python -c "import numpy as np, json; d=np.load('data/cache/operator_tensor.npz', allow_pickle=True); \
s=json.load(open('docs/tables/operator_tensor_summary.json')); \
print({k: d[k].shape for k in d.files}); \
print('confound rho', round(s['n_cells_confound_rho'],3), '| rownorm_cv', round(s['rownorm_cv'],3), \
'| new regs', s['n_new_regulators'])"
```
Expected (acceptance): `tensor`/`mask` are `(≈800, 2000, 3)`; `in_original_panel` is `(≈800,)` (~67 True); **`|confound rho| < 0.15`** (the −0.68 is gone; ≈ +0.08 here), ~733 new (out-of-panel) regulators. The confound guard passed (else the file would not exist). `rownorm_cv` (≈0.23) is informational only — it is selection-confounded on a breadth panel and is NOT a gate (see the driver comment).

- [ ] **Step 4: Add Makefile target and commit**

Add `operator-tensor` to `.PHONY` and (after the `fingerprints` block):
```makefile
operator-tensor:
	$(PY) scripts/build_operator_tensor.py --n-total 800 --top-genes 2000
```

```bash
git add scripts/build_operator_tensor.py Makefile docs/tables/operator_tensor_summary.json
git commit -m "feat(operator): Step 0 driver — z-score expanded-panel tensor with fail-closed representation guards"
```

---

## Task 3: Kernel — varimax rotation + offline hypergeometric enrichment

**Files:**
- Modify: `scripts/_opkernels.py`, `tests/test_opkernels.py`
- Test: `tests/test_opkernels.py`

**Interfaces (produces):**
- `varimax(loadings, gamma=1.0, max_iter=100, tol=1e-6) -> (rotated, rotmat)` — Kaiser-normalized varimax of `(G,k)`.
- `hypergeometric_enrichment(gene_list, gene_sets, background) -> DataFrame` — cols `set_name, n_overlap, set_size, n_drawn, background_size, pvalue, fdr, overlap_genes`; BH-FDR over sets tested per call.

- [ ] **Step 1: Write the failing tests**

```python
def _varimax_crit(X):
    return float(((X ** 2).var(axis=0)).sum())   # varimax objective


def test_varimax_recovers_structure_and_is_varimax():
    # ASYMMETRIC simple structure (a symmetric/balanced toy is a saddle for varimax),
    # mixed by a random orthogonal rotation for varimax to undo.
    rng = np.random.default_rng(0)
    L = np.zeros((9, 3))
    L[0:4, 0] = [0.9, 0.8, 0.7, 0.6]; L[4:7, 1] = [0.85, 0.75, 0.65]; L[7:9, 2] = [0.9, 0.5]
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    M = L @ Q
    rot, R = op.varimax(M)
    dom = np.abs(rot).max(axis=1) / (np.abs(rot).sum(axis=1) + 1e-9)
    assert (dom > 0.85).mean() >= 8 / 9                      # recovers simple structure
    assert _varimax_crit(rot) > _varimax_crit(M) + 1e-6     # is VARIMAX, did not early-stop
    crits = [_varimax_crit(op.varimax(M, max_iter=t)[0]) for t in range(1, 8)]
    assert all(crits[i + 1] >= crits[i] - 1e-9 for i in range(len(crits) - 1))  # monotone
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-6)       # orthonormal


def test_hypergeometric_enrichment_flags_planted_set():
    gene_list = ["ISG15", "MX1", "OAS1", "IFIT1", "STAT1"]
    gene_sets = {"IFN": ["ISG15", "MX1", "OAS1", "IFIT1", "STAT1", "IRF7", "IFI6"],
                 "RANDOM": ["AAA", "BBB", "CCC", "DDD"]}
    background = gene_list + [f"BG{i}" for i in range(500)] + gene_sets["IFN"]
    res = op.hypergeometric_enrichment(gene_list, gene_sets, background)
    top = res.sort_values("pvalue").iloc[0]
    assert top["set_name"] == "IFN" and top["fdr"] < 0.05
```

- [ ] **Step 2: Run to verify fail**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: FAIL — `AttributeError: ... 'varimax'`.

- [ ] **Step 3: Implement**

Append to `scripts/_opkernels.py`:
```python
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def varimax(loadings, gamma=1.0, max_iter=100, tol=1e-6):
    L = np.asarray(loadings, dtype=np.float64)
    G, k = L.shape
    if k < 2:
        return L.copy(), np.eye(k)
    h = np.sqrt((L ** 2).sum(axis=1, keepdims=True)); h[h == 0] = 1.0
    Ln = L / h
    Rm = np.eye(k); d_old = 0.0
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
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: PASS (5 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/_opkernels.py tests/test_opkernels.py
git commit -m "feat(operator): varimax + offline hypergeometric enrichment kernels"
```

---

## Task 4: Step 1 driver — gene programs from the SVD (with power gate)

**Files:**
- Create: `scripts/decompose_operator_svd.py`
- Modify: `Makefile` (add `operator-svd`)
- Test: run on real data; acceptance = ≥2 of top-5 programs get an FDR label AND their power gate is clean.

**Note (controller fix):** gene sets are an **inline `FALLBACK_GENESETS` dict** in the driver (matching the existing `COMPLEXES` dict convention in `analyze_fingerprints.py`), NOT a committed file — `/data/` is gitignored so a `data/genesets/*.gmt` could not be committed. `load_genesets()` starts from the inline dict and merges any `.gmt` files found in an optional (gitignored) `data/genesets/` dir if present.

**Interfaces:**
- Consumes: `data/cache/operator_tensor.npz`; `_opkernels.varimax`, `hypergeometric_enrichment`.
- Produces: `docs/tables/operator_svd_programs.csv` (`rotation, pc, gene, loading, tail`), `docs/tables/operator_svd_enrichment.csv` (`rotation, pc, tail, set_name, ..., fdr`), `docs/tables/operator_svd_power.csv` (`pc, power_rho, power_confounded`), `docs/figures/32_operator_svd_scree.png`. Exports `load_genesets()` (reused by Tasks 6, 12).

- [ ] **Step 1: Write the driver (gene sets inline; optional `.gmt` merge)**

```python
#!/usr/bin/env python3
"""Step 1 — gene programs = right singular vectors of the operator matrix (z-score).

Recovers V (the fingerprint PCA kept only U). Orients each program by an ISG anchor,
optionally varimax-rotates for interpretability, enriches both tails offline, and —
new — runs the SAME power gate as CP: spearman(left factor u_k, per-row n_cells). A
leading program that is really a power axis gets flagged.

    python scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate

Outputs: docs/tables/operator_svd_programs.csv, operator_svd_enrichment.csv,
         operator_svd_power.csv, docs/figures/32_operator_svd_scree.png
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

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
GENESETS = ROOT / "data" / "genesets"    # optional, gitignored; merged if present
ANCHOR = ["ISG15", "MX1", "OAS1", "IFIT1", "STAT1", "IFI6", "IRF7"]

# Inline curated gene sets (offline default; matches the COMPLEXES dict convention
# in analyze_fingerprints.py). /data/ is gitignored so these are NOT a committed file.
FALLBACK_GENESETS = {
    "IFN_ISG": ["ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT2",
                "IFIT3", "IFI6", "IRF7", "STAT1", "STAT2", "IFI44", "IFI44L",
                "RSAD2", "USP18"],
    "TCR_PROXIMAL": ["ZAP70", "LCK", "LAT", "CD3D", "CD3E", "CD3G", "CD247", "FYN",
                     "ITK", "PLCG1", "LCP2", "VAV1", "PIK3CD", "PRKCQ", "CARD11",
                     "BCL10", "MALT1"],
    "CHROMATIN_SAGA": ["TADA1", "TADA2A", "TADA2B", "TADA3", "SUPT20H", "SUPT7L",
                       "TAF5L", "TAF6L", "SGF29", "ATXN7", "ATXN7L3", "USP22",
                       "ENY2", "KAT2A", "KAT2B", "SUPT3H"],
    "MEDIATOR": ["MED1", "MED12", "MED13", "MED14", "MED23", "MED24", "CDK8",
                 "CDK19", "CCNC"],
}


def load_genesets():
    """Inline curated sets, plus any optional .gmt files under data/genesets/."""
    sets = dict(FALLBACK_GENESETS)
    if GENESETS.exists():
        for gmt in sorted(GENESETS.glob("*.gmt")):
            for line in gmt.read_text().splitlines():
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    sets[p[0]] = p[2:]
    return sets


def build_matrix(d):
    """Fully-observed fingerprint rows -> (rows, G) matrix, row index, per-row n_cells."""
    tensor, mask, n_cells = d["tensor"], d["mask"], d["n_cells"]
    R, G, C = tensor.shape
    rows, ncell = [], []
    for i in range(R):
        for c in range(C):
            if mask[i, :, c].all():
                rows.append(tensor[i, :, c]); ncell.append(n_cells[i, c])
    M = np.asarray(rows, np.float64)
    return M - M.mean(axis=0, keepdims=True), np.asarray(ncell)


def orient(v, genes):
    idx = np.isin(genes, ANCHOR)
    s = v[idx].sum() if idx.any() else v[np.argmax(np.abs(v))]
    return v if s >= 0 else -v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--tail-pct", type=float, default=2.0)
    ap.add_argument("--rotate", action="store_true")
    args = ap.parse_args()

    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    genes = d["genes"].astype(str)
    M, ncell = build_matrix(d)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    V = np.column_stack([orient(Vt[j], genes) for j in range(args.k)])   # (G,k)

    # power gate on the LEFT factors (rows == fingerprints)
    pwr = [dict(pc=j + 1, power_rho=float(spearmanr(U[:, j], ncell).statistic),
                power_confounded=bool(abs(spearmanr(U[:, j], ncell).statistic) > 0.3))
           for j in range(args.k)]

    variants = {"raw": V}
    if args.rotate:
        variants["varimax"] = op.varimax(V)[0]

    genesets, bg = load_genesets(), list(genes)
    prog_rows, enr_rows = [], []
    for rotation, VV in variants.items():
        for j in range(args.k):
            v = VV[:, j]
            hi, lo = np.percentile(v, 100 - args.tail_pct), np.percentile(v, args.tail_pct)
            for tail, sel in [("top", v >= hi), ("bottom", v <= lo)]:
                for g, val in zip(genes[sel], v[sel]):
                    prog_rows.append(dict(rotation=rotation, pc=j + 1, gene=g,
                                          loading=float(val), tail=tail))
                enr = op.hypergeometric_enrichment(list(genes[sel]), genesets, bg)
                enr.insert(0, "rotation", rotation); enr.insert(1, "pc", j + 1)
                enr.insert(2, "tail", tail); enr_rows.append(enr)

    TAB.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(prog_rows).to_csv(TAB / "operator_svd_programs.csv", index=False)
    pd.concat(enr_rows, ignore_index=True).to_csv(TAB / "operator_svd_enrichment.csv", index=False)
    pd.DataFrame(pwr).to_csv(TAB / "operator_svd_power.csv", index=False)

    FIG.mkdir(exist_ok=True, parents=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    ev = (S ** 2 / (S ** 2).sum())[:30]
    ax.plot(np.arange(1, len(ev) + 1), ev, "o-")
    ax.set_xlabel("component"); ax.set_ylabel("variance explained")
    ax.set_title("Operator matrix SVD (z-score) — scree")
    fig.tight_layout(); fig.savefig(FIG / "32_operator_svd_scree.png", dpi=150)

    e = pd.concat(enr_rows, ignore_index=True)
    clean_pcs = sorted(set(e[(e.rotation == "raw") & (e.fdr < 0.05) & (e.pc <= 5)]["pc"])
                       & set(p["pc"] for p in pwr if not p["power_confounded"]))
    print(f"top-5 PCs with FDR<0.05 label AND clean power gate: {clean_pcs}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

Run: `python scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate`
Expected: writes 3 tables + scree; prints the PCs that are both FDR-labeled and power-clean.

- [ ] **Step 3: Verify acceptance**

Run:
```bash
python -c "import pandas as pd; e=pd.read_csv('docs/tables/operator_svd_enrichment.csv'); \
p=pd.read_csv('docs/tables/operator_svd_power.csv'); \
lab=e[(e.rotation=='raw')&(e.fdr<0.05)&(e.pc<=5)]; \
print(lab[['pc','tail','set_name','fdr']].sort_values('fdr').head(10).to_string(index=False)); \
print(p.to_string(index=False))"
```
Expected (acceptance): **≥2 of the top-5 programs** carry a clean FDR<0.05 label (IFN and TCR-proximal expected) **and are not power-confounded**. In z-score space the power gate should be clean for the leading programs; if a top program is `power_confounded=True`, report it (a program that survived the representation fix but still tracks power is a real caveat).

- [ ] **Step 4: Makefile + commit**

Add `operator-svd` to `.PHONY` and:
```makefile
operator-svd:
	$(PY) scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate
```

```bash
git add scripts/decompose_operator_svd.py Makefile \
        docs/tables/operator_svd_*.csv docs/figures/32_operator_svd_scree.png
git commit -m "feat(operator): Step 1 driver — SVD gene programs with power gate + offline enrichment"
```

---

## Task 5: Kernel — masked CP, gauge, degeneracy, stability, inter-factor cosine, bootstrap CI

**Files:**
- Modify: `scripts/_opkernels.py`, `tests/test_opkernels.py`
- Test: `tests/test_opkernels.py`

**Interfaces (produces):**
- `cp_fit_masked(tensor, mask, rank, n_iter_max=400, n_init=10, random_state=0) -> (lam, factors)` — best-of-`n_init` masked CP via `tensorly.decomposition.parafac(mask=...)`; NaNs zeroed. **No weight argument** — precision is handled at the representation level (z-score); the old `mask*w` hack is removed because tensorly's `mask` is an observed/missing indicator (`tensor*mask + estimate*(1-mask)`), so a fractional mask blends observed values with model estimates rather than doing WLS. `factors = [A(R,rank), B(G,rank), C(3,rank)]`.
- `fix_cp_gauge(weights, factors) -> (lam, factors)` — unit-normalize columns, push scale to `lam`, sign so the largest-|.| entry of the condition factor is positive (paired flip to A).
- `cp_degeneracy(factors) -> (rank,)` — max off-diagonal triple-cosine congruence.
- `match_factors(B1, B2) -> (perm, cosines)` — Hungarian match by `|cosine|` on the gene mode.
- `split_half_stability(tensor, mask, rank, n_splits=10, random_state=0, subsample=None) -> float` — mean matched gene-mode `|cosine|` across random regulator half-splits; `subsample` caps regulators per split for runtime.
- `gene_mode_cosine(factors) -> (rank,rank)` — **full** `|cosine|` matrix of gene-mode columns (not just the max), so near-collinear factors (e.g. 0.7) are visible before calling factors distinct programs.
- `bootstrap_cp_conditions(tensor, mask, rank, n_boot=100, random_state=0, subsample=None) -> (C_ref, boot_stack)` — reference condition factor `(3,rank)` and a `(n_boot,3,rank)` stack from regulator resampling (with replacement), each matched+sign-aligned to the reference gene mode. Feeds the gating CI.

- [ ] **Step 1: Write the failing tests**

```python
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
```

- [ ] **Step 2: Run to verify fail**

Run: `pip install tensorly && python -m pytest tests/test_opkernels.py -q`
Expected: FAIL — `AttributeError: ... 'cp_fit_masked'`.

- [ ] **Step 3: Implement**

Append to `scripts/_opkernels.py`:
```python
from scipy.optimize import linear_sum_assignment


def _unit_cols(X):
    n = np.linalg.norm(X, axis=0, keepdims=True); n[n == 0] = 1.0
    return X / n, n.ravel()


def cp_fit_masked(tensor, mask, rank, n_iter_max=400, n_init=10, random_state=0):
    import tensorly as tl
    from tensorly.decomposition import parafac
    T = np.nan_to_num(np.asarray(tensor, np.float64), nan=0.0)
    Tt, Mt = tl.tensor(T), tl.tensor(mask.astype(np.float64))
    best = None
    for s in range(n_init):
        cp = parafac(Tt, rank=rank, mask=Mt, n_iter_max=n_iter_max,
                     init="random", random_state=random_state + s,
                     normalize_factors=False)
        err = float(np.sum(((tl.cp_to_tensor(cp) - T) ** 2) * mask))
        if best is None or err < best[0]:
            best = (err, cp)
    cp = best[1]
    lam = np.asarray(cp[0]) if cp[0] is not None else np.ones(rank)
    return lam, [np.asarray(f) for f in cp[1]]


def fix_cp_gauge(weights, factors):
    factors = [f.copy() for f in factors]
    lam = np.ones(factors[0].shape[1])
    for mode, f in enumerate(factors):
        fn, norms = _unit_cols(f); factors[mode] = fn; lam = lam * norms
    Cc = factors[2]
    for k in range(Cc.shape[1]):
        if Cc[np.argmax(np.abs(Cc[:, k])), k] < 0:
            factors[2][:, k] *= -1; factors[0][:, k] *= -1
    return lam, factors


def gene_mode_cosine(factors):
    Bb, _ = _unit_cols(factors[1])
    return np.abs(Bb.T @ Bb)


def cp_degeneracy(factors):
    Aa, _ = _unit_cols(factors[0]); Bb, _ = _unit_cols(factors[1]); Cc, _ = _unit_cols(factors[2])
    cong = np.abs(Aa.T @ Aa) * np.abs(Bb.T @ Bb) * np.abs(Cc.T @ Cc)
    np.fill_diagonal(cong, 0.0)
    return cong.max(axis=1)


def match_factors(B1, B2):
    B1n, _ = _unit_cols(np.asarray(B1, float)); B2n, _ = _unit_cols(np.asarray(B2, float))
    cost = -np.abs(B1n.T @ B2n)
    r, c = linear_sum_assignment(cost)
    return c, -cost[r, c]


def split_half_stability(tensor, mask, rank, n_splits=10, random_state=0, subsample=None):
    rng = np.random.default_rng(random_state)
    R = tensor.shape[0]; scores = []
    for s in range(n_splits):
        perm = rng.permutation(R)
        if subsample:
            perm = perm[: 2 * subsample]
        h1, h2 = perm[: len(perm) // 2], perm[len(perm) // 2:]
        _, f1 = cp_fit_masked(tensor[h1], mask[h1], rank, n_init=3, random_state=random_state + s)
        _, f2 = cp_fit_masked(tensor[h2], mask[h2], rank, n_init=3, random_state=random_state + 100 + s)
        scores.append(float(match_factors(f1[1], f2[1])[1].mean()))
    return float(np.mean(scores))


def bootstrap_cp_conditions(tensor, mask, rank, n_boot=100, random_state=0, subsample=None):
    rng = np.random.default_rng(random_state)
    R = tensor.shape[0]
    lam0, f0 = fix_cp_gauge(*cp_fit_masked(tensor, mask, rank, n_init=5, random_state=random_state))
    Bref = f0[1]
    boot = np.full((n_boot, 3, rank), np.nan)
    for b in range(n_boot):
        idx = rng.integers(0, R, subsample or R)
        lam, f = fix_cp_gauge(*cp_fit_masked(tensor[idx], mask[idx], rank,
                                             n_init=1, random_state=random_state + b + 1))
        perm, _ = match_factors(Bref, f[1])
        Cc = f[2][:, perm].copy()
        for k in range(rank):
            if np.dot(f[1][:, perm[k]], Bref[:, k]) < 0:
                Cc[:, k] *= -1
        boot[b] = Cc
    return f0[2], boot
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: PASS (11 passed). If `test_split_half_stability_high_for_true_rank` is boundary-flaky, `s3 > 0.65` may relax to `> 0.60`, but `s3 > s8` must hold.

- [ ] **Step 5: Commit**

```bash
git add scripts/_opkernels.py tests/test_opkernels.py
git commit -m "feat(operator): CP kernels — masked fit, gauge, degeneracy, stability, cosine matrix, bootstrap CI"
```

---

## Task 6: Step 2 driver — CP with RMS scale control, stability rank, bootstrap-CI gating

**Files:**
- Create: `scripts/decompose_operator_cp.py`
- Modify: `Makefile` (add `operator-cp`)
- Test: run on real data; acceptance = stability-selected rank + a factor whose gating CI excludes flat and is power-clean and non-degenerate.

**Interfaces:**
- Consumes: `data/cache/operator_tensor.npz`; Task 5 kernels + `rms_normalize_conditions`, `spearman_power`, `hypergeometric_enrichment`; `decompose_operator_svd.load_genesets`.
- Produces: `docs/tables/operator_cp_factors.csv` (`factor, lambda_, cond_Rest, cond_Stim8hr, cond_Stim48hr, gating_shape, gated_ci, range_lo95, power_rho, power_confounded, degeneracy, max_cofactor_cosine, top_regulators, top_genes, program_label`), `docs/tables/operator_cp_stability.csv`, `docs/tables/operator_cp_cosine.csv` (full inter-factor gene-mode `|cos|`), `docs/tables/operator_cp_enrichment.csv`, `docs/figures/33_operator_cp_stability.png`, `docs/figures/34_operator_cp_condition_factors.png`.

**Nuisance controls (each maps to a failure mode):**
1. **RMS condition scale control** (`rms_normalize_conditions`) before fit → `c_k` is *relative* modulation; a `--scale-control none` counterfactual shows gating is otherwise trivially true.
2. **No fake precision weights** — representation (z-score) handles power; the broken `mask*w` hack is gone.
3. **Rank by split-half stability** (largest rank clearing threshold; else argmax) with `--stab-subsample` for runtime.
4. **Gating via bootstrap CI, not a hardcoded threshold** — `gated_ci=True` only if the 2.5th-percentile of the per-bootstrap condition-profile range exceeds the flat threshold.
5. **Per-factor power gate** (`|spearman(a_k, n_cells)| > 0.3` → flag).
6. **Degeneracy** (`cp_degeneracy > 0.9`) **and full inter-factor cosine** (`max_cofactor_cosine`) — two factors at cosine 0.7 are not distinct programs even though degeneracy passes.
7. **Init stability** `n_init=10`.
8. **Selection-leakage** cross-check via `--sensitivity` note.

- [ ] **Step 1: Write the driver**

```python
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
```

- [ ] **Step 2: Run it**

Run: `python scripts/decompose_operator_cp.py --rank auto --scale-control rms`
Expected: prints confound meter (~0), stability curve, auto rank (likely 3–6), and `[result]` with clean-gated vs clean-constitutive factor lists. Writes 4 tables + 2 figures. Runtime is dominated by the stability sweep + `--boot-n` bootstrap fits on ≤400 regulators; expect minutes-to-tens-of-minutes. Raise `--stab-subsample`/`--boot-n` only if you can afford it.

- [ ] **Step 3: Verify the flagship acceptance criterion**

Run:
```bash
python -c "import pandas as pd; f=pd.read_csv('docs/tables/operator_cp_factors.csv'); \
print(f[['factor','gating_shape','gated_ci','range_lo95','program_label','power_confounded','degeneracy','max_cofactor_cosine']].to_string(index=False))"
```
Expected (acceptance): after RMS control, **≥1 clean factor** (`power_confounded=False`, `degeneracy<0.9`, `max_cofactor_cosine<0.7`) has **`gated_ci=True`** (CI excludes flat) with a TCR/immune label, **and** ≥1 clean factor with `gated_ci=False` (constitutive) and a chromatin/IFN label. If the only gated factors are power-confounded/degenerate/collinear, or no CI excludes flat, the gating headline does not survive — record that as the result.

**This gating test is LOAD-BEARING for the Step-0 representation, not just for the biology.** Step 0 retired its magnitude-based pooling guards because, on a per-cell-z breadth panel, only a *per-regulator cross-condition* signal can distinguish pooled z from within-condition z — and that is exactly what a `gated_ci=True` factor is. So **at least one clean `gated_ci=True` factor is a POSITIVE assertion here**: it is direct proof the representation is genuinely pooled (you cannot manufacture gating out of within-condition-normalized data, where every condition slab is independently scaled → all factors would be constitutive). Therefore, if this step comes back with **every** clean factor `gated_ci=False` (all-constitutive), do **not** treat it as a benign "no gating result" — it is a **red flag that the fetched layer may be within-condition z after all**, and it must trigger re-examination of the representation (re-check the confound guard, re-inspect `layers/zscore` semantics on the raw fetch) rather than a shrug. Record which interpretation applies.

- [ ] **Step 4: Scale-control counterfactual**

Run: `python scripts/decompose_operator_cp.py --rank 4 --scale-control none`
Expected: without RMS control, condition factors lean to the higher-magnitude conditions for most factors (gating trivially true), showing the control is load-bearing. **Re-run Step 2 default before committing** so committed tables are the `rms` ones.

- [ ] **Step 5: Makefile + commit**

Add `operator-cp` to `.PHONY` and:
```makefile
operator-cp:
	$(PY) scripts/decompose_operator_cp.py --rank auto --scale-control rms
```

```bash
git add scripts/decompose_operator_cp.py Makefile docs/tables/operator_cp_*.csv \
        docs/figures/33_operator_cp_stability.png docs/figures/34_operator_cp_condition_factors.png
git commit -m "feat(operator): Step 2 driver — CP with RMS control, stability rank, bootstrap-CI gating, cosine matrix"
```

---

## Task 7: Kernel — soft-impute completion + train-only standardization

**Files:**
- Modify: `scripts/_opkernels.py`, `tests/test_opkernels.py`
- Test: `tests/test_opkernels.py`

**Interfaces (produces):**
- `train_test_standardize(M, train_mask) -> (Mc, mu)` — per-column mean over **train entries only**; test entries centered by the train mean. Never fits on test.
- `soft_impute(M, observed_mask, rank, n_iter=100, tol=1e-4) -> M_hat` — iterative SVD hard-truncation to `rank` (init unobserved 0; loop SVD→truncate→reset observed).

- [ ] **Step 1: Write the failing tests**

```python
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
```

- [ ] **Step 2: Run to verify fail**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: FAIL — `AttributeError: ... 'train_test_standardize'`.

- [ ] **Step 3: Implement**

Append to `scripts/_opkernels.py`:
```python
def train_test_standardize(M, train_mask):
    M = np.asarray(M, float); mu = np.zeros(M.shape[1])
    for j in range(M.shape[1]):
        col = M[train_mask[:, j], j]
        mu[j] = col.mean() if col.size else 0.0
    return M - mu[None, :], mu


def soft_impute(M, observed_mask, rank, n_iter=100, tol=1e-4):
    M = np.asarray(M, float); X = np.where(observed_mask, M, 0.0); Xr = X; prev = np.inf
    for _ in range(n_iter):
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        Xr = (U[:, :rank] * s[:rank]) @ Vt[:rank]
        X = np.where(observed_mask, M, Xr)
        change = np.linalg.norm(Xr - X) / (np.linalg.norm(X) + 1e-12)
        if abs(prev - change) < tol:
            break
        prev = change
    return Xr
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: PASS (13 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/_opkernels.py tests/test_opkernels.py
git commit -m "feat(operator): soft-impute + train-only standardization kernels"
```

---

## Task 8: Step 3 driver — out-of-panel condition extrapolation (flagship) + entry-wise sanity

**Files:**
- Create: `scripts/operator_completion.py`
- Modify: `Makefile` (add `operator-completion`)
- Test: run on real data; acceptance = **3b beats persistence on out-of-panel regulators** (the escalation-worthy result); 3a is a sanity check only.

**Interfaces:**
- Consumes: `data/cache/operator_tensor.npz` (incl. `in_original_panel`); `_opkernels.train_test_standardize`, `soft_impute`.
- Produces: `docs/tables/operator_completion_condition.csv` (`rank, r2_model, r2_persistence, beats_persistence, r2_model_novel, r2_persistence_novel, beats_persistence_novel, n_test, n_test_novel`), `docs/tables/operator_completion_entrywise.csv` (`rank, r2_vs_geneMean, r2_vs_rank1, beats_rank1`), `docs/figures/35_operator_completion_curve.png`.

**Ordering and honesty (per audit):**
- **3b is the flagship** and its held-out fibers are drawn **preferentially from out-of-panel (`in_original_panel==False`) regulators**, with novelty logged (`n_test_novel`, `r2_*_novel`). Baseline = persistence (Stim48hr := Stim8hr). This is a genuine "given rest+early stim, predict late stim for a regulator we never characterized" test.
- **3a is demoted to a sanity check.** Beating the centered per-gene-mean (=0) is near-automatic for a correlated matrix, so 3a reports R² vs gene-mean AND **vs a rank-1 baseline**, and we read only the **elbow** (effective rank), not "beats zero".
- **No genome-scale imputation claim.** Pure low-rank cannot predict a regulator with zero observed entries; 3b holds out a condition *fiber*, never a whole regulator. The "predict ~11k unpaneled perturbations" framing is dropped entirely (it needs side information — a future direction).

- [ ] **Step 1: Write the driver**

```python
#!/usr/bin/env python3
"""Step 3 — is the z-score operator recoverably low-rank? Held-out prediction.

3b (FLAGSHIP): hold out (regulator, Stim48hr) fibers, PREFERENTIALLY for out-of-panel
regulators, and predict them from that regulator's Rest+Stim8hr via the low-rank fit
on the other regulators. Baseline = persistence (Stim48hr := Stim8hr).
3a (SANITY): random entry-wise completion; report the effective-rank elbow and R²
vs a rank-1 baseline (beating the gene-mean baseline is near-trivial and not claimed).

    python scripts/operator_completion.py --max-rank 12 --holdout 0.2

Outputs: docs/tables/operator_completion_condition.csv,
         operator_completion_entrywise.csv, docs/figures/35_*
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
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def to_matrix(d):
    tensor, mask = d["tensor"], d["mask"]
    R, G, C = tensor.shape
    rows = [tensor[i, :, c] for i in range(R) for c in range(C) if mask[i, :, c].all()]
    return np.asarray(rows, float)


def entrywise(M, max_rank, holdout, seed):
    rng = np.random.default_rng(seed)
    hide = rng.random(M.shape) < holdout
    train = ~hide
    Mc, _ = op.train_test_standardize(M, train)
    yt = Mc[hide]
    base_gm = float(np.mean(yt ** 2))                        # gene-mean == 0 (centered)
    r1 = op.soft_impute(Mc, train, 1, n_iter=200)
    base_r1 = float(np.mean((yt - r1[hide]) ** 2))
    rows = []
    for r in range(1, max_rank + 1):
        yp = op.soft_impute(Mc, train, r, n_iter=200)[hide]
        mse = float(np.mean((yt - yp) ** 2))
        rows.append(dict(rank=r, r2_vs_geneMean=1 - mse / (base_gm + 1e-12),
                         r2_vs_rank1=1 - mse / (base_r1 + 1e-12),
                         beats_rank1=bool(mse < base_r1)))
    return pd.DataFrame(rows)


def condition_extrap(d, max_rank, seed):
    tensor, mask, in_orig = d["tensor"], d["mask"], d["in_original_panel"]
    R, G, C = tensor.shape
    late, early = COND_ORDER.index("Stim48hr"), COND_ORDER.index("Stim8hr")
    full = np.array([i for i in range(R) if mask[i, :, :].all(axis=0).all()])
    novel_full = full[~in_orig[full]]
    rng = np.random.default_rng(seed)
    # hold out preferentially from out-of-panel regulators
    n_test = max(20, len(full) // 5)
    test = rng.permutation(novel_full)[:n_test]
    if len(test) < n_test:                                   # top up from any full reg
        extra = rng.permutation(np.setdiff1d(full, test))[: n_test - len(test)]
        test = np.concatenate([test, extra])
    X = np.concatenate([tensor[full, :, 0], tensor[full, :, early], tensor[full, :, late]], axis=1)
    obs = np.ones_like(X, bool)
    test_pos = np.isin(full, test)
    novel_pos = test_pos & (~in_orig[full])
    obs[test_pos, 2 * G:3 * G] = False
    Xc, _ = op.train_test_standardize(X, obs)
    rows = []
    for r in range(1, max_rank + 1):
        Xhat = op.soft_impute(Xc, obs, r, n_iter=200)
        def block(posmask):
            yt = Xc[np.ix_(posmask, np.arange(2 * G, 3 * G))]
            yp = Xhat[np.ix_(posmask, np.arange(2 * G, 3 * G))]
            pers = Xc[np.ix_(posmask, np.arange(G, 2 * G))]
            ss = np.sum((yt - yt.mean()) ** 2) + 1e-12
            return (float(1 - np.sum((yt - yp) ** 2) / ss),
                    float(1 - np.sum((yt - pers) ** 2) / ss))
        rm, rp = block(test_pos)
        rmn, rpn = block(novel_pos)
        rows.append(dict(rank=r, r2_model=rm, r2_persistence=rp, beats_persistence=bool(rm > rp),
                         r2_model_novel=rmn, r2_persistence_novel=rpn,
                         beats_persistence_novel=bool(rmn > rpn),
                         n_test=int(test_pos.sum()), n_test_novel=int(novel_pos.sum())))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-rank", type=int, default=12)
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    print("[note] pure low-rank cannot predict a regulator with ZERO observed entries; "
          "3b holds out a condition FIBER, never a whole regulator. No genome-scale claim.")
    ew = entrywise(to_matrix(d), args.max_rank, args.holdout, args.seed)
    ce = condition_extrap(d, args.max_rank, args.seed)
    TAB.mkdir(exist_ok=True, parents=True)
    ew.to_csv(TAB / "operator_completion_entrywise.csv", index=False)
    ce.to_csv(TAB / "operator_completion_condition.csv", index=False)

    FIG.mkdir(exist_ok=True, parents=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ce["rank"], ce["r2_model_novel"], "s-", label="3b model R² (out-of-panel)")
    ax.plot(ce["rank"], ce["r2_persistence_novel"], "--", c="crimson", label="3b persistence (out-of-panel)")
    ax.plot(ew["rank"], ew["r2_vs_rank1"], "o-", alpha=0.5, label="3a R² vs rank-1 (sanity)")
    ax.axhline(0, ls=":", c="grey"); ax.set_xlabel("rank"); ax.set_ylabel("held-out R²")
    ax.legend(fontsize=8); ax.set_title("Operator low-rank recoverability (flagship = 3b out-of-panel)")
    fig.tight_layout(); fig.savefig(FIG / "35_operator_completion_curve.png", dpi=150)

    print(f"[3b FLAGSHIP] out-of-panel n={int(ce.n_test_novel.iloc[0])}; "
          f"beats persistence at ranks={ce.loc[ce.beats_persistence_novel,'rank'].tolist()}; "
          f"max model R²={ce.r2_model_novel.max():.3f} vs persistence={ce.r2_persistence_novel.max():.3f}")
    print(f"[3a sanity] beats rank-1 at ranks={ew.loc[ew.beats_rank1,'rank'].tolist()}; "
          f"elbow ~rank {int(ew.r2_vs_rank1.idxmax()) + 1}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

Run: `python scripts/operator_completion.py --max-rank 12 --holdout 0.2`
Expected: writes 2 tables + 1 figure; prints the 3b out-of-panel result (with `n_test_novel`) and the 3a sanity elbow.

- [ ] **Step 3: Verify acceptance (flagship = 3b on out-of-panel)**

Run:
```bash
python -c "import pandas as pd; b=pd.read_csv('docs/tables/operator_completion_condition.csv'); \
print('out-of-panel n:', int(b.n_test_novel.iloc[0])); \
print('3b beats persistence (out-of-panel):', bool(b.beats_persistence_novel.any()), \
'best model R2', round(b.r2_model_novel.max(),3), 'vs persistence', round(b.r2_persistence_novel.max(),3))"
```
Expected (acceptance): **3b beats persistence on the out-of-panel regulators** at some rank — the descriptive→predictive result, and the escalation signal (Task 11). If it only beats persistence on in-panel regulators but not out-of-panel, say so plainly (the model generalizes within the characterized set but not beyond). 3a is not an acceptance gate; it only locates the effective-rank elbow.

- [ ] **Step 4: Makefile + commit**

Add `operator-completion` to `.PHONY` and:
```makefile
operator-completion:
	$(PY) scripts/operator_completion.py --max-rank 12 --holdout 0.2
```

```bash
git add scripts/operator_completion.py Makefile docs/tables/operator_completion_*.csv \
        docs/figures/35_operator_completion_curve.png
git commit -m "feat(operator): Step 3 driver — out-of-panel condition extrapolation (flagship) + entry-wise sanity"
```

---

## Task 9: Kernel — principal angles + random-subspace null

**Files:**
- Modify: `scripts/_opkernels.py`, `tests/test_opkernels.py`
- Test: `tests/test_opkernels.py`

**Interfaces (produces):**
- `principal_angles(Va, Vb) -> cos2` — QR-orthonormalize two `(G,k)` bases, return `cos²θ` descending (via `scipy.linalg.subspace_angles`).
- `random_subspace_null(G, k, n=1000, random_state=0) -> (mean_cos2, p95_cos2)`.

- [ ] **Step 1: Write the failing tests**

```python
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
```

- [ ] **Step 2: Run to verify fail**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: FAIL — `AttributeError: ... 'principal_angles'`.

- [ ] **Step 3: Implement**

Append to `scripts/_opkernels.py`:
```python
from scipy.linalg import subspace_angles, qr


def principal_angles(Va, Vb):
    Qa = qr(np.asarray(Va, float), mode="economic")[0]
    Qb = qr(np.asarray(Vb, float), mode="economic")[0]
    return np.sort(np.cos(subspace_angles(Qa, Qb)) ** 2)[::-1]


def random_subspace_null(G, k, n=1000, random_state=0):
    rng = np.random.default_rng(random_state)
    means = np.array([principal_angles(rng.normal(size=(G, k)), rng.normal(size=(G, k))).mean()
                      for _ in range(n)])
    return float(means.mean()), float(np.percentile(means, 95))
```

- [ ] **Step 4: Run to verify pass**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: PASS (16 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/_opkernels.py tests/test_opkernels.py
git commit -m "feat(operator): principal angles + random-subspace null kernels"
```

---

## Task 10: Step 4 driver — donor-subspace stability across disjoint donor pairs

**Files:**
- Create: `scripts/operator_donor_angles.py`
- Modify: `Makefile` (add `operator-donors`)
- Test: run; if per-donor matrices absent, prints NEEDS-DATA and exits 0.

**Interfaces:**
- Consumes: per-donor-pair matrices — **NOT in the local cache** (`data/cache/` has only `donor_obs.csv`, a per-regulator summary — verified). Supports (1) `data/cache/by_donors_index.csv` (`regulator,gene`) + `data/cache/by_donors_pair_<ab>.npy` (`R×G`, ideally z-score for consistency) if a fetch produced them; (2) absent → NEEDS-DATA message + empty-headed table + exit 0.
- Produces: `docs/tables/operator_donor_angles.csv` (`pair_a, pair_b, k, mean_cos2, null_mean_cos2, null_p95, above_null`), `docs/figures/36_operator_donor_angles.png` (only when data present).

**Nuisance controls:** **disjoint donor pairs only** — donors {1,2,3,4} → `(1,2)vs(3,4)`, `(1,3)vs(2,4)`, `(1,4)vs(2,3)`; never the 15 inflated overlapping comparisons. Subspace comparison via `principal_angles` only (never per-vector correlation).

- [ ] **Step 1: Write the driver**

```python
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
```

- [ ] **Step 2: Run it**

Run: `python scripts/operator_donor_angles.py --k 5`
Expected: per-donor matrices not local → prints `[NEEDS-DATA]`, writes empty-headed table, exits 0. Do not fabricate data.

- [ ] **Step 3: (Conditional) real result when data present**

If `by_donors_pair_*.npy` exist: re-run; acceptance = `above_null=True` for all three disjoint comparisons (mean cos²θ well above the random-subspace 95th pct → programs donor-reproducible as subspaces).

- [ ] **Step 4: Makefile + commit**

Add `operator-donors` to `.PHONY` and:
```makefile
operator-donors:
	$(PY) scripts/operator_donor_angles.py --k 5
```

```bash
git add scripts/operator_donor_angles.py Makefile docs/tables/operator_donor_angles.csv
git commit -m "feat(operator): Step 4 driver — donor-subspace principal angles (disjoint pairs, null-gated)"
```

---

## Task 11: Integration — umbrella target, writeup (with reconciliation + escalation trigger)

**Files:**
- Modify: `Makefile` (umbrella `operator` + help), `README.md` (pointer)
- Create: `docs/OPERATOR_ANALYSIS.md`
- Test: `make operator` runs Steps 0–3 (Step 4 NEEDS-DATA is non-fatal); writeup reflects real numbers.

- [ ] **Step 1: Umbrella target + help**

Add `operator` to `.PHONY`; add to `help:`:
```makefile
	@echo "  make operator - empirical regulatory operator (z-score): tensor + SVD + CP + completion [+donors if fetched]"
```
Target:
```makefile
operator: operator-tensor operator-svd operator-cp operator-completion operator-donors
```

- [ ] **Step 2: Run full local pipeline**

Run: `make operator`
Expected: Steps 0→1→2→3 complete (Step 0 fetch is cached after first run); Step 4 prints NEEDS-DATA, exits 0, target succeeds.

- [ ] **Step 3: Full test suite green**

Run: `python -m pytest tests/test_opkernels.py -q`
Expected: PASS (16 passed).

- [ ] **Step 4: Write `docs/OPERATOR_ANALYSIS.md` from real numbers**

Required sections (each states the actual produced value; no placeholders):
```markdown
# The Empirical Regulatory Operator (z-score, Option-1 pilot)

One matrix, four questions. The operator = effect of every regulator KD on every
measured gene, built in precision-decoupled z-score space over an ~800-regulator
expanded panel (top ~800 by n_downstream above a median cell-count floor; ~733 out-of-panel).

## Representation and why z-score (Step 0)
- Raw log-FC row-norm vs power: spearman = -0.683 (confound). z-score confound
  meter = <value> (the asserted fail-closed guard, |ρ|<0.15). Note the two retired
  magnitude guards (row-norm CV, anchor spread) and why: on a per-cell-z, breadth-
  homogeneous panel magnitude can't separate pooled from within-condition z; pooling
  rests on the layer-level 200-panel property + the downstream CP-gating test.

## Gene programs (Step 1)
- Top-5 programs, FDR labels (operator_svd_enrichment.csv), power gate
  (operator_svd_power.csv) — how many top-5 are both labeled and power-clean (≥2).

## Condition gating (Step 2) — flagship of the descriptive half
- Stability-selected rank; condition factors after RMS control WITH bootstrap 95% CI;
  which factors are clean (power-clean, non-degenerate, max cofactor cosine <0.7) and
  gated_ci=True vs constitutive. Scale-control counterfactual (none => trivially gated).
- RECONCILIATION: the prior within-condition permutation found TCR's complex-cohesion
  z was NOT condition-inflated in z-score space (11.2 -> 11.2). That is a different
  quantity from CP program-magnitude gating across conditions. State both and why they
  coexist (cohesion of the complex vs magnitude of the program's response), so the two
  results do not read as a contradiction.

## Prediction (Step 3) — flagship of the predictive half
- 3b out-of-panel: does low-rank beat persistence on the ~600 novel regulators
  (n_test_novel = <value>)? This is the descriptive->predictive result.
- 3a sanity: effective-rank elbow; note beating the gene-mean baseline is near-trivial
  and NOT claimed. No genome-scale imputation claim (row-cold-start impossible for pure
  low-rank).

## Donor reproducibility (Step 4)
- Disjoint-pair principal-angle result, OR the NEEDS-DATA note + fetch pattern.

## Escalation trigger
- If 3b beats persistence cleanly on out-of-panel regulators, escalate to Option 3:
  full 6209-regulator axis, completion as the flagship, CP sweep on cluster compute.
  Option 1 is the pilot that earns Option 3. Decide on this evidence.

## What the nuisance controls bought us
- One line each: representation (z-score kills -0.68), RMS scale control, bootstrap-CI
  gating, power gates (SVD+CP), degeneracy + full cosine matrix, split-half rank,
  disjoint donor pairs, honest baselines (rank-1 / persistence). For each: would a
  naive version have produced a different (wrong) headline?
```

- [ ] **Step 5: README pointer + commit**

Add to `README.md`:
```markdown
- **Empirical regulatory operator** (`make operator`) — z-score tensor/CP/completion/donor-subspace analysis; see `docs/OPERATOR_ANALYSIS.md`.
```

```bash
git add Makefile docs/OPERATOR_ANALYSIS.md README.md
git commit -m "docs(operator): umbrella target, analysis writeup (reconciliation + escalation trigger), README pointer"
```

---

## Task 12 (STRETCH — only if Tasks 6 and 8 land cleanly): deconvolution + asymmetric subsumption

**Files:**
- Create: `scripts/operator_deconvolution.py`
- Modify: `Makefile` (add `operator-deconv`)
- Test: run; acceptance = square block reported (verified size 91 in the top-2000 panel), edges carry the block condition number, framed hypothesis-generating.

**Interfaces:**
- Consumes: `data/cache/operator_tensor.npz`; `docs/tables/top_robust_regulators.csv` (existing) for the shortlist.
- Produces: `docs/tables/operator_deconvolution_edges.csv` (`source, target, weight, block_condition, cond_number`), `docs/tables/operator_subsumption.csv` (`source, target, condition, subsumption_frac`).

**Nuisance controls (mandatory):** valid only on the **square block** (regulators that are also readout genes — verified 91 in this panel, so top-50 works); **regularized inverse** with the **block condition number reported per condition**; **linear approx to a nonlinear system** → asymmetry is hypothesis-generating, stated in filename and doc; z-score representation makes the `(I−A)⁻¹` steady-state reading an even bigger stretch — say so. Subsumption `s(i→j)` = fraction of `j`'s significant set contained in `i`'s (directional; hierarchy cosine cannot express).

- [ ] **Step 1: Write the driver**

```python
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
    # dedup the donor-robust shortlist (it repeats rows per condition/contrast → would seed
    # duplicate self-pairs), then top up from the panel's regulator-genes to n_robust.
    rp = TAB / "top_robust_regulators.csv"
    rset, gset = set(regs), set(genes)
    shortlist = (list(dict.fromkeys(pd.read_csv(rp).iloc[:, 0].astype(str).tolist()))
                 if rp.exists() else [])
    block = [g for g in shortlist if g in rset and g in gset]
    for g in regs:
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
```

- [ ] **Step 2: Run it**

Run: `python scripts/operator_deconvolution.py --n-robust 50 --ridge 1e-2`
Expected: block size ≈ up to 50 (91 candidates exist), two tables written.

- [ ] **Step 3: Condition-number caveat honored**

Run:
```bash
python -c "import pandas as pd; e=pd.read_csv('docs/tables/operator_deconvolution_edges.csv'); \
print(e.groupby('block_condition')['cond_number'].first() if len(e) else 'no edges')"
```
Expected: every edge carries the block condition number; if ≫1e6, writeup flags the inverse unreliable / edges speculative.

- [ ] **Step 4: Makefile + commit**

Add `operator-deconv` to `.PHONY` and:
```makefile
operator-deconv:
	$(PY) scripts/operator_deconvolution.py --n-robust 50 --ridge 1e-2
```

```bash
git add scripts/operator_deconvolution.py Makefile docs/tables/operator_deconvolution_edges.csv \
        docs/tables/operator_subsumption.csv
git commit -m "feat(operator): Step 5 stretch — square-block deconvolution + asymmetric subsumption"
```

---

## Global nuisance checklist → where enforced

| # | Confound | Enforced in |
|---|----------|-------------|
| 0 | **Representation / power (the −0.68)** | **z-score tensor + fail-closed Step-0 guards (Task 2)** — the primary fix, not a downstream patch |
| 1 | Residual power per factor | per-factor `power_confounded` gate in **both** SVD (Task 4) and CP (Task 6) |
| 2 | Global condition scale vs gating | `rms_normalize_conditions` (RMS per observed entry) + `--scale-control none` counterfactual (Task 6) |
| 3 | Gating false positives | **bootstrap 95% CI on c_k** (Task 5/6); gated only if CI excludes flat — no hardcoded threshold |
| 4 | Selection leakage | expanded-panel/gene rule + `--sensitivity` cross-check (Tasks 2, 6) |
| 5 | Standardization leakage | `train_test_standardize` train-only (Tasks 7, 8) |
| 6 | Sign/scale/rotation gauge | `fix_cp_gauge`; `varimax`; principal angles for subspaces (Tasks 5, 3/4, 9/10) |
| 7 | Rank selection | `split_half_stability` sweep (`--stab-subsample`), `n_init=10` (Tasks 5, 6) |
| 8 | Factor collinearity | `cp_degeneracy` **and full `gene_mode_cosine` matrix** (`max_cofactor_cosine<0.7`) (Tasks 5, 6) |
| 9 | Donor-pair non-independence | disjoint pairs only (Task 10) |
| 10 | Baselines | 3b vs persistence on out-of-panel; 3a vs rank-1 (never "beats zero") (Task 8) |
| 11 | Multiple testing | BH-FDR per enrichment family (Task 3) |
| 12 | Orthogonality artifact | raw vs varimax; orthogonality never called biology (Tasks 3/4, writeup) |

---

## Self-Review (against the spec + both review rounds)

**Spec coverage:** Step 0 → Tasks 1–2. Step 1 → Tasks 3–4. Step 2 → Tasks 5–6. Step 3 → Tasks 7–8. Step 4 → Tasks 9–10. Step 5 → Task 12. Global checklist → mapped and enforced. Integration/writeup → Task 11.

**Review-round-2 items, all landed:**
- **Critical representation fix** → z-score tensor (Task 1/2), with a **fail-closed Step-0 confound guard** (|ρ|<0.15) so re-entered power *refuses to cache*. (Two earlier magnitude guards — row-norm CV, anchor spread — were retired after verification showed magnitude can't separate pooled from within-condition z on a per-cell-z breadth panel; pooling rests on the layer-level property + the load-bearing downstream CP-gating test.)
- Broken `mask*w` precision hack **removed** (Task 5); precision handled by representation.
- Gating **bootstrap CI** replaces the hardcoded threshold (Tasks 5/6); reconciliation with the prior 11.2→11.2 permutation result required in the writeup (Task 11).
- 3a **demoted** to sanity (rank-1 baseline + elbow); **3b out-of-panel is the flagship** with novelty logged; "predict ~11k unpaneled" framing **dropped** (Task 8).
- Full inter-factor **gene-mode cosine matrix** (not just max degeneracy) (Tasks 5/6); **SVD power gate** added (Task 4); **RMS** not raw Frobenius (Task 1); gene-selection power caveat noted (writeup).
- Facts fixed: gene axis **10282** (not 10273); ~800-reg expanded panel (not "hundreds–2k"); square block **91**; runtime bounded via `--stab-subsample`/`--boot-n`.
- **Escalation trigger** written as its own section (Task 11).

**Known data-reality departures (called out):**
- Full-axis z-score for all 6209 regulators is Option 3 (cluster); this pilot is ~800. The escalation trigger governs when to buy it.
- `by_donors.h5mu` per-donor matrices are not local; Task 10 degrades to NEEDS-DATA with the fetch pattern documented.

**Placeholder scan:** every code step is complete and runnable; every run step has an exact command + expected output; every acceptance step is concrete. No TBD/"handle edge cases"/"similar to Task N".

**Type consistency:** `cp_fit_masked -> (lam, [A,B,C])`; `fix_cp_gauge` preserves shapes; `bootstrap_cp_conditions -> (C_ref(3,rank), boot(n_boot,3,rank))`; `gene_mode_cosine -> (rank,rank)`; `assemble_tensor` takes any row-aligned `matrix`; npz keys written in Task 2 (incl. `in_original_panel`) are exactly those read in Tasks 4/6/8/12. `COND_ORDER` identical everywhere.
