# Empirical Regulatory Operator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover the gene-program, condition-modulation, and predictive structure of the regulatory operator matrix `L` (log-FC of every gene under every regulator KD) that the existing fingerprint PCA computed but discarded — via tensor decomposition, low-rank completion, and donor-subspace stability, with every result gated on the confounds that decide whether it is biology or an artifact.

**Architecture:** One 3-way tensor `T[regulator, gene, condition]` is assembled once from the local log-FC cache (Step 0). Everything after is a question about that one object: its right factors (Step 1, gene programs), its CP decomposition (Step 2, `regulator ⊗ gene ⊗ condition` — where "gated vs constitutive" becomes the shape of one factor), its recoverability by low rank (Step 3, held-out prediction), and the reproducibility of its program subspace across donors (Step 4). Pure numerical kernels live in a shared, unit-tested module `scripts/_operator.py` (mirroring the existing `scripts/_figstyle.py` helper convention); one driver script per step follows the repo's analysis-script idiom (argparse → `docs/tables/*.csv|json` + `docs/figures/*.png`), with a Makefile target each.

**Tech Stack:** Python 3, numpy, scipy, statsmodels, pandas, matplotlib (all already in `requirements.txt`); adds `tensorly` (CP/masked-CP) and `pytest` (kernel unit tests). scikit-learn 1.9 is already installed. Offline-first: no network for Steps 0–3; Step 4 needs per-donor matrices that are **not** in the local cache (see Global Constraints).

## Global Constraints

- **Data is fully local for Steps 0–3.** `data/cache/log_fc.f32.npy` is the entire contrast×gene log-FC matrix, shape `(33983, 10273)` float32, C-order. Row *i* corresponds to row *i* of `data/cache/fingerprint_obs.csv`; column *j* to row *j* of `data/cache/fingerprint_var.csv`. Do not stream or download anything for Steps 0–3.
- **`fingerprint_obs.csv` columns (verbatim):** `index, target_contrast_gene_name, culture_condition, target_contrast, n_cells_target, n_up_genes, n_down_genes, n_total_de_genes, n_downstream, ontarget_effect_size, ontarget_significant, distal_offtarget_flag, low_target_gex, neighboring_gene_KD, single_guide_estimate, n_guides, guide_correlation_all, guide_correlation_signif, donor_correlation_all_mean, donor_correlation_hits_mean`. The KD gate is `ontarget_significant == True`. The per-cell power covariate is `n_cells_target`. There is **no `lfcSE` column locally** — precision weighting (Step 2) falls back to `n_cells_target` (see that task).
- **`fingerprint_var.csv` columns:** `gene_name, gene_id`.
- **Condition axis order is fixed everywhere:** `["Rest", "Stim8hr", "Stim48hr"]` (call it `COND_ORDER`; matches `scripts/analyze_fingerprints.py`).
- **`ontarget_significant` is a string in the CSV** ("True"/"False"); parse with `.astype(str).str.strip().eq("True")`, never rely on truthiness of the raw column.
- **Output locations (repo convention):** tables → `docs/tables/operator_*.csv|json`; figures → `docs/figures/<NN>_operator_*.png` (next free two-digit prefix after the existing 31); writeup → `docs/OPERATOR_ANALYSIS.md`. Path constants in every driver: `ROOT = Path(__file__).resolve().parent.parent`, `CACHE = ROOT/"data"/"cache"`, `TAB = ROOT/"docs"/"tables"`, `FIG = ROOT/"docs"/"figures"`.
- **Determinism:** every stochastic routine takes an explicit `random_state` / seed; default `0`. Use `np.random.default_rng(seed)`.
- **The standing confound meter** — `spearman(||T[i,:,c]||, n_cells[i,c])` — is computed in Step 0 and re-run per extracted factor in every later step. Any factor whose regulator-mode loading correlates with `n_cells` at `|rho| > 0.3` is flagged `power_confounded=True` in its output row. This is not optional decoration; it is the accept/reject gate.
- **No network in Steps 0–3.** Step 4 and Step 5's stretch tasks are the only ones that may touch the network, and each states it explicitly and degrades gracefully when the data is absent.
- Follow the existing script style: `#!/usr/bin/env python3`, module docstring with a usage line, `matplotlib.use("Agg")` before `pyplot`, argparse `main()`, write CSVs with `index=False`.

---

## File Structure

**New shared kernel module (pure functions, unit-tested):**
- `scripts/_operator.py` — all pure numerical kernels: `assemble_tensor`, `frobenius_normalize_conditions`, `spearman_power`, `cp_fit_masked`, `fix_cp_gauge`, `cp_degeneracy`, `split_half_stability`, `match_factors`, `varimax`, `soft_impute`, `train_test_standardize`, `principal_angles`, `random_subspace_null`, `hypergeometric_enrichment`. No I/O, no argparse, no plotting.

**New driver scripts (one per step, repo analysis-script idiom):**
- `scripts/build_operator_tensor.py` — Step 0. Builds and caches the tensor + mask + `n_cells` + confound meter.
- `scripts/decompose_operator_svd.py` — Step 1. Gene programs from the matrix SVD (+ varimax + offline enrichment).
- `scripts/decompose_operator_cp.py` — Step 2. Masked/weighted CP, global-scale control, stability rank selection, per-factor confound test, gating readout.
- `scripts/operator_completion.py` — Step 3. Entry-wise completion (3a) and condition-mode extrapolation (3b) with honest baselines.
- `scripts/operator_donor_angles.py` — Step 4. Principal angles across **disjoint** donor pairs vs random-subspace null.
- `scripts/operator_deconvolution.py` — Step 5 (stretch, optional). Square-block network deconvolution + asymmetric subsumption.

**Tests (new; pytest introduced for the numerical kernels):**
- `tests/test_operator.py` — synthetic-data unit tests for every kernel in `scripts/_operator.py`.

**Config / docs:**
- `requirements.txt` — add `tensorly`, `pytest` (in the optional block).
- `Makefile` — add targets `operator-tensor`, `operator-svd`, `operator-cp`, `operator-completion`, `operator-donors`, and an umbrella `operator`.
- `docs/OPERATOR_ANALYSIS.md` — the writeup (final task).

**Cache artifact produced by Step 0 (consumed by Steps 1–4):**
- `data/cache/operator_tensor.npz` — keys: `tensor` `(R,G,3)` float32, `mask` `(R,G,3)` bool, `regulators` `(R,)` str, `genes` `(G,)` str, `conditions` `(3,)` str, `n_cells` `(R,3)` float32.

---

## Task 1: Kernel — tensor assembly + condition normalization + confound meter

**Files:**
- Create: `scripts/_operator.py`
- Create: `scripts/__init__.py` (empty package marker)
- Create: `tests/test_operator.py`
- Modify: `requirements.txt` (add `pytest`)
- Test: `tests/test_operator.py`

**Interfaces:**
- Produces:
  - `assemble_tensor(log_fc, obs, gene_idx, cond_order) -> (tensor, mask, regulators, n_cells)` where `log_fc` is the `(N, G_all)` memmap, `obs` a DataFrame with columns `target_contrast_gene_name, culture_condition, n_cells_target, ontarget_significant_bool, row`, `gene_idx` an int array of selected gene columns, `cond_order` a list of 3 condition names. Returns `tensor (R,G,3) float32` (NaN where unobserved), `mask (R,G,3) bool`, `regulators (R,) object`, `n_cells (R,3) float32`.
  - `frobenius_normalize_conditions(tensor, mask) -> (tensor_norm, scales)` — divides each condition slab `[:, :, c]` by its Frobenius norm over observed entries; returns per-condition scale `(3,)`.
  - `spearman_power(tensor, mask, n_cells) -> float` — Spearman rho between per-cell slab norm `||T[i,:,c]||` and `n_cells[i,c]` over observed `(i,c)` cells.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_operator.py
import numpy as np
import pandas as pd
import pytest
from scripts import _operator as op


def _toy_obs():
    # 3 regulators x 3 conditions; reg C missing Stim48hr.
    rows = []
    r = 0
    plan = {"A": ["Rest", "Stim8hr", "Stim48hr"],
            "B": ["Rest", "Stim8hr", "Stim48hr"],
            "C": ["Rest", "Stim8hr"]}
    for g, conds in plan.items():
        for c in conds:
            rows.append(dict(target_contrast_gene_name=g, culture_condition=c,
                             n_cells_target=100.0 + r, ontarget_significant_bool=True,
                             row=r))
            r += 1
    return pd.DataFrame(rows)


def test_assemble_tensor_shapes_and_mask():
    obs = _toy_obs()
    G_all = 5
    log_fc = np.arange(len(obs) * G_all, dtype=np.float32).reshape(len(obs), G_all)
    gene_idx = np.array([0, 2, 4])
    tensor, mask, regulators, n_cells = op.assemble_tensor(
        log_fc, obs, gene_idx, ["Rest", "Stim8hr", "Stim48hr"])
    assert tensor.shape == (3, 3, 3)      # 3 regs, 3 genes, 3 conds
    assert mask.shape == (3, 3, 3)
    assert list(regulators) == ["A", "B", "C"]
    ci = list(regulators).index("C")
    assert mask[ci, :, 2].sum() == 0                       # C has no Stim48hr
    assert mask[list(regulators).index("A"), :, 2].sum() == 3
    assert np.isnan(tensor[ci, :, 2]).all()


def test_frobenius_normalize_conditions_unit_norm():
    rng = np.random.default_rng(0)
    tensor = rng.normal(size=(4, 6, 3)).astype(np.float32)
    mask = np.ones_like(tensor, dtype=bool)
    tn, scales = op.frobenius_normalize_conditions(tensor, mask)
    for c in range(3):
        assert np.isclose(np.linalg.norm(tn[:, :, c]), 1.0, atol=1e-5)
    assert scales.shape == (3,)


def test_spearman_power_detects_planted_confound():
    rng = np.random.default_rng(1)
    R, G = 20, 8
    n_cells = rng.uniform(50, 5000, size=(R, 3))
    tensor = np.zeros((R, G, 3), dtype=np.float32)
    for c in range(3):
        for i in range(R):
            tensor[i, :, c] = n_cells[i, c] / 1000.0   # norm monotone in n_cells
    mask = np.ones((R, G, 3), dtype=bool)
    rho = op.spearman_power(tensor, mask, n_cells)
    assert rho > 0.9
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_operator.py -q`
Expected: FAIL — `ModuleNotFoundError` / `AttributeError: module 'scripts._operator' has no attribute 'assemble_tensor'`.

- [ ] **Step 3: Create the package marker and kernel module with these three functions**

Create empty `scripts/__init__.py` (so `from scripts import _operator` works):
```python
```

Create `scripts/_operator.py`:
```python
#!/usr/bin/env python3
"""Pure numerical kernels for the empirical regulatory operator analysis.

No I/O, no argparse, no plotting — everything here is unit-tested on synthetic
data in tests/test_operator.py. Drivers (build_operator_tensor.py etc.) import
these. Mirrors the shared-helper convention of scripts/_figstyle.py.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import spearmanr


def assemble_tensor(log_fc, obs, gene_idx, cond_order):
    """Stack per-(regulator, condition) log-FC rows into a dense-with-mask tensor.

    obs must have columns: target_contrast_gene_name, culture_condition,
    n_cells_target, ontarget_significant_bool, row. Only rows with
    ontarget_significant_bool == True fill cells; unobserved cells are NaN in
    `tensor` and False in `mask`.
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
        tensor[i, :, k] = log_fc[int(rec["row"]), gene_idx]
        mask[i, :, k] = True
        n_cells[i, k] = rec["n_cells_target"]
    return tensor, mask, np.array(regulators, dtype=object), n_cells


def frobenius_normalize_conditions(tensor, mask):
    """Scale each condition slab to unit Frobenius norm over observed entries.

    Removes the global 'Stim conditions simply have larger DE' size effect so a
    CP condition factor reads as *relative* modulation, not raw magnitude.
    """
    out = tensor.copy()
    C = tensor.shape[2]
    scales = np.zeros(C, dtype=np.float64)
    for c in range(C):
        vals = tensor[:, :, c][mask[:, :, c]]
        s = float(np.sqrt(np.sum(vals.astype(np.float64) ** 2)))
        scales[c] = s
        if s > 0:
            out[:, :, c] = out[:, :, c] / s
    return out, scales


def spearman_power(tensor, mask, n_cells):
    """Spearman rho between per-cell slab norm and n_cells over observed cells."""
    norms, ncell = [], []
    R, _, C = tensor.shape
    for i in range(R):
        for c in range(C):
            if mask[i, :, c].any():
                vals = tensor[i, :, c][mask[i, :, c]]
                norms.append(float(np.linalg.norm(vals)))
                ncell.append(float(n_cells[i, c]))
    if len(norms) < 3:
        return float("nan")
    return float(spearmanr(norms, ncell).statistic)
```

Append to `requirements.txt` under the optional section:
```
# Optional — empirical regulatory operator analysis (`make operator`)
# install only if you will run the operator tensor/CP/completion analysis:
#   pip install tensorly pytest
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_operator.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/__init__.py scripts/_operator.py tests/test_operator.py requirements.txt
git commit -m "feat(operator): tensor assembly, condition normalization, confound meter kernels"
```

---

## Task 2: Step 0 driver — build and cache the operator tensor

**Files:**
- Create: `scripts/build_operator_tensor.py`
- Modify: `Makefile` (add `operator-tensor` target)
- Test: run on real data; acceptance = printed confound meter + cached artifact.

**Interfaces:**
- Consumes: `_operator.assemble_tensor`, `frobenius_normalize_conditions`, `spearman_power`.
- Produces: `data/cache/operator_tensor.npz` (keys listed in File Structure) and `docs/tables/operator_tensor_summary.json` with fields `n_regulators, n_genes, conditions, n_cells_confound_rho, observed_cells, inclusion_rule, gene_selection`.

- [ ] **Step 1: Write the driver**

```python
#!/usr/bin/env python3
"""Step 0 — build the regulator x gene x condition log-FC tensor.

Inclusion rule (default): regulators KD-significant (ontarget_significant) in ALL
three conditions, so the tensor is dense along the condition mode. Gene axis: top
--top-genes by variance across the pooled selected slab. Alternative rules for the
selection-leakage sensitivity check are exposed as flags.

    python scripts/build_operator_tensor.py --top-genes 2000
    python scripts/build_operator_tensor.py --top-genes 2000 --inclusion any --gene-var-cond Rest

Outputs: data/cache/operator_tensor.npz, docs/tables/operator_tensor_summary.json
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scripts import _operator as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def load_obs():
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv").reset_index(drop=True)
    obs["row"] = np.arange(len(obs))
    obs["ontarget_significant_bool"] = (
        obs["ontarget_significant"].astype(str).str.strip().eq("True"))
    var = pd.read_csv(CACHE / "fingerprint_var.csv")
    return obs, var


def select_regulators(obs, inclusion):
    sig = obs[obs["ontarget_significant_bool"]]
    conds = (sig.groupby("target_contrast_gene_name")["culture_condition"]
             .apply(lambda s: set(s) & set(COND_ORDER)))
    if inclusion == "all":
        keep = conds[conds.apply(lambda s: set(COND_ORDER).issubset(s))].index
    else:  # "any" — significant in >=1 condition (ragged)
        keep = conds[conds.apply(lambda s: len(s) >= 1)].index
    return set(keep)


def select_genes(log_fc, obs, keep_regs, top_genes, gene_var_cond):
    m = obs["ontarget_significant_bool"] & obs["target_contrast_gene_name"].isin(keep_regs)
    if gene_var_cond != "pooled":
        m &= obs["culture_condition"].eq(gene_var_cond)
    rows = obs.loc[m, "row"].to_numpy()
    slab = log_fc[rows, :]
    var = slab.var(axis=0)
    return np.argsort(var)[::-1][:top_genes]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-genes", type=int, default=2000)
    ap.add_argument("--inclusion", choices=["all", "any"], default="all")
    ap.add_argument("--gene-var-cond", default="pooled",
                    choices=["pooled", "Rest", "Stim8hr", "Stim48hr"])
    args = ap.parse_args()

    obs, var = load_obs()
    log_fc = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
    keep_regs = select_regulators(obs, args.inclusion)
    gene_idx = select_genes(log_fc, obs, keep_regs, args.top_genes, args.gene_var_cond)
    genes = var.loc[gene_idx, "gene_name"].to_numpy()

    obs_sel = obs[obs["target_contrast_gene_name"].isin(keep_regs)].copy()
    tensor, mask, regulators, n_cells = op.assemble_tensor(
        log_fc, obs_sel, gene_idx, COND_ORDER)
    rho = op.spearman_power(tensor, mask, n_cells)

    CACHE.mkdir(exist_ok=True, parents=True)
    np.savez(CACHE / "operator_tensor.npz",
             tensor=tensor, mask=mask,
             regulators=regulators.astype(str), genes=genes.astype(str),
             conditions=np.array(COND_ORDER), n_cells=n_cells)
    summary = dict(n_regulators=int(len(regulators)), n_genes=int(len(genes)),
                   conditions=COND_ORDER, n_cells_confound_rho=rho,
                   observed_cells=int(mask.any(axis=1).sum()),
                   inclusion_rule=args.inclusion,
                   gene_selection=f"top{args.top_genes}var:{args.gene_var_cond}")
    TAB.mkdir(exist_ok=True, parents=True)
    (TAB / "operator_tensor_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python scripts/build_operator_tensor.py --top-genes 2000`
Expected: prints a JSON summary; `n_regulators` in the low hundreds to ~2k; writes `data/cache/operator_tensor.npz` and `docs/tables/operator_tensor_summary.json`. The spec's ρ≈−0.68 is the expectation for `n_cells_confound_rho`; record whatever it is.

- [ ] **Step 3: Verify the artifact loads and shapes agree**

Run:
```bash
python -c "import numpy as np; d=np.load('data/cache/operator_tensor.npz', allow_pickle=True); \
print({k: d[k].shape for k in d.files}); print('genes', d['genes'][:5]); \
print('conds', d['conditions'])"
```
Expected: `tensor` and `mask` are `(R, 2000, 3)`; `regulators (R,)`; `genes (2000,)`; `conditions ('Rest','Stim8hr','Stim48hr')`.

- [ ] **Step 4: Run the selection-leakage sensitivity build (used by Task 6's cross-check)**

Run: `python scripts/build_operator_tensor.py --top-genes 2000 --inclusion any --gene-var-cond Rest`
Expected: succeeds, overwrites the cache with the ragged/Rest-variance version. **Then re-run the Step 2 default build** to leave the `all`/`pooled` tensor as the committed default. (Task 6 re-runs both and compares factors; this step only proves both builds succeed.)

- [ ] **Step 5: Add the Makefile target and commit**

Add `operator-tensor` to the `.PHONY` line, and add this block after the `fingerprints` block:
```makefile
operator-tensor:
	$(PY) scripts/build_operator_tensor.py --top-genes 2000
```

```bash
git add scripts/build_operator_tensor.py Makefile docs/tables/operator_tensor_summary.json
git commit -m "feat(operator): Step 0 driver — build+cache regulator x gene x condition tensor"
```

---

## Task 3: Kernel — varimax rotation + offline hypergeometric enrichment

**Files:**
- Modify: `scripts/_operator.py` (add `varimax`, `hypergeometric_enrichment`)
- Modify: `tests/test_operator.py`
- Test: `tests/test_operator.py`

**Interfaces:**
- Produces:
  - `varimax(loadings, gamma=1.0, max_iter=100, tol=1e-6) -> (rotated, rotmat)` — Kaiser-normalized varimax on a `(G, k)` loading matrix; returns rotated loadings `(G, k)` and rotation `(k, k)`.
  - `hypergeometric_enrichment(gene_list, gene_sets, background) -> DataFrame` — columns `set_name, n_overlap, set_size, n_drawn, background_size, pvalue, fdr, overlap_genes`; BH-FDR (`statsmodels.stats.multitest.multipletests(method="fdr_bh")`) over all sets tested in the call.

- [ ] **Step 1: Write the failing tests**

```python
def test_varimax_recovers_axis_aligned_structure():
    L = np.zeros((6, 2))
    L[:3, 0] = 1.0
    L[3:, 1] = 1.0
    theta = np.pi / 6
    Rm = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    mixed = L @ Rm.T
    rot, _ = op.varimax(mixed)
    dominance = np.abs(rot).max(axis=1) / (np.abs(rot).sum(axis=1) + 1e-9)
    assert (dominance > 0.9).mean() >= 5 / 6


def test_hypergeometric_enrichment_flags_planted_set():
    gene_list = ["ISG15", "MX1", "OAS1", "IFIT1", "STAT1"]
    gene_sets = {"IFN": ["ISG15", "MX1", "OAS1", "IFIT1", "STAT1", "IRF7", "IFI6"],
                 "RANDOM": ["AAA", "BBB", "CCC", "DDD"]}
    background = gene_list + [f"BG{i}" for i in range(500)] + gene_sets["IFN"]
    res = op.hypergeometric_enrichment(gene_list, gene_sets, background)
    top = res.sort_values("pvalue").iloc[0]
    assert top["set_name"] == "IFN"
    assert top["fdr"] < 0.05
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_operator.py -q`
Expected: FAIL — `AttributeError: module 'scripts._operator' has no attribute 'varimax'`.

- [ ] **Step 3: Implement the two functions**

Append to `scripts/_operator.py`:
```python
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


def varimax(loadings, gamma=1.0, max_iter=100, tol=1e-6):
    """Kaiser-normalized varimax rotation of a (G, k) loading matrix."""
    L = np.asarray(loadings, dtype=np.float64)
    G, k = L.shape
    if k < 2:
        return L.copy(), np.eye(k)
    h = np.sqrt((L ** 2).sum(axis=1, keepdims=True))
    h[h == 0] = 1.0
    Ln = L / h                       # Kaiser normalization
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
    """Hypergeometric over-representation of gene_list in each gene_set, BH-FDR."""
    bg = set(background)
    drawn = set(gene_list) & bg
    M = len(bg)
    n_drawn = len(drawn)
    rows = []
    for name, members in gene_sets.items():
        setg = set(members) & bg
        K = len(setg)
        overlap = drawn & setg
        x = len(overlap)
        p = hypergeom.sf(x - 1, M, K, n_drawn) if K > 0 and n_drawn > 0 else 1.0
        rows.append(dict(set_name=name, n_overlap=x, set_size=K, n_drawn=n_drawn,
                         background_size=M, pvalue=float(p),
                         overlap_genes=",".join(sorted(overlap))))
    df = pd.DataFrame(rows)
    if len(df):
        df["fdr"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    return df.sort_values("pvalue").reset_index(drop=True)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_operator.py -q`
Expected: PASS (5 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/_operator.py tests/test_operator.py
git commit -m "feat(operator): varimax rotation + offline hypergeometric enrichment kernels"
```

---

## Task 4: Step 1 driver — gene programs from the matrix SVD (right factors)

**Files:**
- Create: `scripts/decompose_operator_svd.py`
- Create: `data/genesets/operator_genesets.gmt` (small offline curated fallback)
- Modify: `Makefile` (add `operator-svd`)
- Test: run on real data; acceptance = ≥2 of top-5 programs get an FDR-significant label.

**Interfaces:**
- Consumes: `data/cache/operator_tensor.npz`; `_operator.varimax`, `hypergeometric_enrichment`.
- Produces: `docs/tables/operator_svd_programs.csv` (`rotation, pc, gene, loading, tail`), `docs/tables/operator_svd_enrichment.csv` (`rotation, pc, tail, set_name, ...`), `docs/figures/32_operator_svd_scree.png`.
- Also produces the importable helper `load_genesets()` (reused by Tasks 6 and 12).

**Design note:** Build the matrix `M` from the tensor's fully-observed fingerprint rows flattened to `(rows, G)` (one row per observed regulator×condition slab), mean-centered per gene, then SVD. `V = Vt[:k].T` columns are the gene programs — exactly the right factor the fingerprint PCA discarded.

- [ ] **Step 1: Create the offline gene-set fallback file**

`data/genesets/operator_genesets.gmt` — tab-separated `set_name<TAB>description<TAB>gene1<TAB>gene2...`. Minimal curated sets so the step runs fully offline; extend by dropping any MSigDB `.gmt` into this directory.
```
IFN_ISG	interferon_stimulated	ISG15	MX1	MX2	OAS1	OAS2	OAS3	IFIT1	IFIT2	IFIT3	IFI6	IRF7	STAT1	STAT2	IFI44	IFI44L	RSAD2	USP18
TCR_PROXIMAL	tcr_signaling	ZAP70	LCK	LAT	CD3D	CD3E	CD3G	CD247	FYN	ITK	PLCG1	LCP2	VAV1	PIK3CD	PRKCQ	CARD11	BCL10	MALT1
CHROMATIN_SAGA	saga_complex	TADA1	TADA2A	TADA2B	TADA3	SUPT20H	SUPT7L	TAF5L	TAF6L	SGF29	ATXN7	ATXN7L3	USP22	ENY2	KAT2A	KAT2B	SUPT3H
MEDIATOR	mediator_complex	MED1	MED12	MED13	MED14	MED23	MED24	CDK8	CDK19	CCNC
```

- [ ] **Step 2: Write the driver**

```python
#!/usr/bin/env python3
"""Step 1 — gene programs = right singular vectors of the operator matrix.

The fingerprint PCA kept the LEFT factors (perturbation fingerprints) and
discarded V. This recovers V: rank genes by loading per program, orient by a
fixed anchor (ISGs load positive), optionally varimax-rotate the top-k for an
interpretable (non-orthogonal) view, and enrich both tails offline.

    python scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate

Outputs: docs/tables/operator_svd_programs.csv, operator_svd_enrichment.csv,
         docs/figures/32_operator_svd_scree.png
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts import _operator as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
GENESETS = ROOT / "data" / "genesets"
ANCHOR = ["ISG15", "MX1", "OAS1", "IFIT1", "STAT1", "IFI6", "IRF7"]  # orient +


def load_genesets():
    sets = {}
    for gmt in sorted(GENESETS.glob("*.gmt")):
        for line in gmt.read_text().splitlines():
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = parts[2:]
    return sets


def build_matrix(d):
    tensor, mask = d["tensor"], d["mask"]
    R, G, C = tensor.shape
    rows = []
    for i in range(R):
        for c in range(C):
            if mask[i, :, c].all():          # fully-observed fingerprint rows only
                rows.append(tensor[i, :, c])
    M = np.asarray(rows, dtype=np.float64)
    M = M - M.mean(axis=0, keepdims=True)    # per-gene centering
    return M


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
    M = build_matrix(d)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    V = Vt[:args.k].T                        # (G, k) gene programs
    V = np.column_stack([orient(V[:, j], genes) for j in range(args.k)])

    variants = {"raw": V}
    if args.rotate:
        variants["varimax"] = op.varimax(V)[0]

    genesets = load_genesets()
    background = list(genes)
    prog_rows, enr_rows = [], []
    for rotation, VV in variants.items():
        for j in range(args.k):
            v = VV[:, j]
            hi = np.percentile(v, 100 - args.tail_pct)
            lo = np.percentile(v, args.tail_pct)
            for tail, sel in [("top", v >= hi), ("bottom", v <= lo)]:
                tg = genes[sel]
                for g, val in zip(tg, v[sel]):
                    prog_rows.append(dict(rotation=rotation, pc=j + 1, gene=g,
                                          loading=float(val), tail=tail))
                enr = op.hypergeometric_enrichment(list(tg), genesets, background)
                enr.insert(0, "rotation", rotation)
                enr.insert(1, "pc", j + 1)
                enr.insert(2, "tail", tail)
                enr_rows.append(enr)

    TAB.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(prog_rows).to_csv(TAB / "operator_svd_programs.csv", index=False)
    pd.concat(enr_rows, ignore_index=True).to_csv(
        TAB / "operator_svd_enrichment.csv", index=False)

    FIG.mkdir(exist_ok=True, parents=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    ev = (S ** 2 / (S ** 2).sum())[:30]
    ax.plot(np.arange(1, len(ev) + 1), ev, "o-")
    ax.set_xlabel("component"); ax.set_ylabel("variance explained")
    ax.set_title("Operator matrix SVD — scree")
    fig.tight_layout(); fig.savefig(FIG / "32_operator_svd_scree.png", dpi=150)

    sig = pd.concat(enr_rows, ignore_index=True)
    sig = sig[(sig["rotation"] == "raw") & (sig["fdr"] < 0.05) & (sig["pc"] <= 5)]
    print(f"top-5 PCs with an FDR<0.05 label: {sorted(sig['pc'].unique())}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run it on real data**

Run: `python scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate`
Expected: writes the two CSVs and the scree figure; prints the list of top-5 PCs that earned an FDR<0.05 label.

- [ ] **Step 4: Verify the acceptance criterion**

Run:
```bash
python -c "import pandas as pd; e=pd.read_csv('docs/tables/operator_svd_enrichment.csv'); \
s=e[(e.rotation=='raw')&(e.fdr<0.05)&(e.pc<=5)]; \
print(s[['pc','tail','set_name','fdr']].sort_values('fdr').head(12).to_string(index=False))"
```
Expected (acceptance): **≥2 of the top-5 programs** carry a clean FDR<0.05 label — interferon (IFN_ISG) and a TCR-proximal program are the expected hits. If fewer than 2, do not proceed to interpretation; record it as a real negative (SVD programs not cleanly interpretable at this gene selection).

- [ ] **Step 5: Add Makefile target and commit**

Add `operator-svd` to `.PHONY` and:
```makefile
operator-svd:
	$(PY) scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate
```

```bash
git add scripts/decompose_operator_svd.py data/genesets/operator_genesets.gmt Makefile \
        docs/tables/operator_svd_programs.csv docs/tables/operator_svd_enrichment.csv \
        docs/figures/32_operator_svd_scree.png
git commit -m "feat(operator): Step 1 driver — gene programs from operator SVD (+varimax, offline enrichment)"
```

---

## Task 5: Kernel — masked CP fit + gauge fixing + degeneracy + split-half stability

**Files:**
- Modify: `scripts/_operator.py` (add `cp_fit_masked`, `fix_cp_gauge`, `cp_degeneracy`, `match_factors`, `split_half_stability`)
- Modify: `tests/test_operator.py`
- Test: `tests/test_operator.py`

**Interfaces:**
- Produces:
  - `cp_fit_masked(tensor, mask, rank, n_iter_max=400, n_init=10, random_state=0, weights=None) -> (lam, factors)` — best-of-`n_init` masked CP via `tensorly.decomposition.parafac(mask=...)`; NaNs in `tensor` are zeroed before the call (mask handles them). Optional `weights` is a `(R, C)` per-cell array folded into the mask. Returns CP weights `lam (rank,)` and `factors = [A (R,rank), B (G,rank), C (Cc,rank)]`.
  - `fix_cp_gauge(weights, factors) -> (lam, factors)` — push all scale into `lam[k]`, unit-normalize each mode's columns, fix sign so the largest-|.| entry of the **condition** factor column (`factors[2]`) is positive, propagating the paired flip to `factors[0]` so the rank-1 term is unchanged.
  - `cp_degeneracy(factors) -> array (rank,)` — per-factor max off-diagonal triple-cosine congruence (near 1 ⇒ degenerate).
  - `match_factors(B1, B2) -> (perm, cosines)` — Hungarian match of `B2` columns to `B1` by `|cosine|`; returns permutation and matched `|cosine|`.
  - `split_half_stability(tensor, mask, rank, n_splits=10, random_state=0) -> float` — mean matched `|cosine|` on the gene mode across random regulator half-splits.

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
    T, (A, B, Cc) = _synth_cp()
    mask = np.ones_like(T, dtype=bool)
    _, factors = op.cp_fit_masked(T, mask, rank=3, n_init=5, random_state=0)
    _, cos = op.match_factors(B, factors[1])
    assert cos.mean() > 0.95


def test_cp_masked_ignores_hidden_cells():
    T, _ = _synth_cp()
    mask = np.ones_like(T, dtype=bool)
    mask[::2, :, 2] = False           # hide half the Stim48hr slab
    _, factors = op.cp_fit_masked(T, mask, rank=3, n_init=5, random_state=0)
    assert factors[0].shape == (40, 3) and factors[2].shape == (3, 3)


def test_fix_cp_gauge_unit_norm_and_sign():
    T, _ = _synth_cp()
    mask = np.ones_like(T, dtype=bool)
    lam, factors = op.cp_fit_masked(T, mask, rank=3, n_init=3, random_state=1)
    _, f2 = op.fix_cp_gauge(lam, factors)
    for F in f2:
        assert np.allclose(np.linalg.norm(F, axis=0), 1.0, atol=1e-5)
    Cc = f2[2]
    for k in range(Cc.shape[1]):
        assert Cc[np.argmax(np.abs(Cc[:, k])), k] > 0


def test_split_half_stability_high_for_true_rank():
    T, _ = _synth_cp(R=80, noise=0.03)
    mask = np.ones_like(T, dtype=bool)
    stab3 = op.split_half_stability(T, mask, rank=3, n_splits=4, random_state=0)
    stab8 = op.split_half_stability(T, mask, rank=8, n_splits=4, random_state=0)
    assert stab3 > 0.7
    assert stab3 > stab8         # true rank more reproducible than an inflated one
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pip install tensorly && python -m pytest tests/test_operator.py -q`
Expected: FAIL — `AttributeError: ... 'cp_fit_masked'`.

- [ ] **Step 3: Implement the CP kernels**

Append to `scripts/_operator.py`:
```python
from scipy.optimize import linear_sum_assignment


def _unit_cols(X):
    n = np.linalg.norm(X, axis=0, keepdims=True)
    n[n == 0] = 1.0
    return X / n, n.ravel()


def cp_fit_masked(tensor, mask, rank, n_iter_max=400, n_init=10,
                  random_state=0, weights=None):
    import tensorly as tl
    from tensorly.decomposition import parafac
    T = np.nan_to_num(np.asarray(tensor, dtype=np.float64), nan=0.0)
    m = mask.astype(np.float64)
    if weights is not None:                      # (R, C) -> broadcast over genes
        w = weights[:, None, :]
        m = m * (w / (np.nanmax(w) + 1e-12))
    Tt, Mt = tl.tensor(T), tl.tensor(m)
    best = None
    for s in range(n_init):
        cp = parafac(Tt, rank=rank, mask=Mt, n_iter_max=n_iter_max,
                     init="random", random_state=random_state + s,
                     normalize_factors=False)
        rec = tl.cp_to_tensor(cp)
        err = float(np.sum(((rec - T) ** 2) * m))
        if best is None or err < best[0]:
            best = (err, cp)
    cp = best[1]
    lam = np.asarray(cp[0]) if cp[0] is not None else np.ones(rank)
    factors = [np.asarray(f) for f in cp[1]]
    return lam, factors


def fix_cp_gauge(weights, factors):
    factors = [f.copy() for f in factors]
    rank = factors[0].shape[1]
    lam = np.ones(rank)
    for mode, f in enumerate(factors):
        fn, norms = _unit_cols(f)
        factors[mode] = fn
        lam = lam * norms
    Cc = factors[2]
    for k in range(rank):
        if Cc[np.argmax(np.abs(Cc[:, k])), k] < 0:
            factors[2][:, k] *= -1
            factors[0][:, k] *= -1        # paired flip keeps the rank-1 term intact
    return lam, factors


def cp_degeneracy(factors):
    Aa, _ = _unit_cols(factors[0])
    Bb, _ = _unit_cols(factors[1])
    Cc, _ = _unit_cols(factors[2])
    cong = np.abs(Aa.T @ Aa) * np.abs(Bb.T @ Bb) * np.abs(Cc.T @ Cc)
    np.fill_diagonal(cong, 0.0)
    return cong.max(axis=1)


def match_factors(B1, B2):
    B1n, _ = _unit_cols(np.asarray(B1, float))
    B2n, _ = _unit_cols(np.asarray(B2, float))
    cost = -np.abs(B1n.T @ B2n)
    r, c = linear_sum_assignment(cost)
    return c, -cost[r, c]


def split_half_stability(tensor, mask, rank, n_splits=10, random_state=0):
    rng = np.random.default_rng(random_state)
    R = tensor.shape[0]
    scores = []
    for s in range(n_splits):
        perm = rng.permutation(R)
        h1, h2 = perm[: R // 2], perm[R // 2:]
        _, f1 = cp_fit_masked(tensor[h1], mask[h1], rank, n_init=3,
                              random_state=random_state + s)
        _, f2 = cp_fit_masked(tensor[h2], mask[h2], rank, n_init=3,
                              random_state=random_state + 100 + s)
        _, cos = match_factors(f1[1], f2[1])   # gene mode is shared across splits
        scores.append(float(cos.mean()))
    return float(np.mean(scores))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_operator.py -q`
Expected: PASS (9 passed). If `test_split_half_stability_high_for_true_rank` is flaky at the boundary, loosening `stab3 > 0.7` to `> 0.65` is acceptable — but `stab3 > stab8` must hold.

- [ ] **Step 5: Commit**

```bash
git add scripts/_operator.py tests/test_operator.py
git commit -m "feat(operator): masked CP + gauge fix + degeneracy + split-half stability kernels"
```

---

## Task 6: Step 2 driver — CP decomposition with scale control, stability rank, confound gate

**Files:**
- Create: `scripts/decompose_operator_cp.py`
- Modify: `Makefile` (add `operator-cp`)
- Test: run on real data; acceptance = a stability-selected rank + a gated factor that survives the scale control.

**Interfaces:**
- Consumes: `data/cache/operator_tensor.npz`; all Task 5 kernels + `frobenius_normalize_conditions`, `spearman_power`; `decompose_operator_svd.load_genesets`.
- Produces: `docs/tables/operator_cp_factors.csv` (`factor, lambda_, cond_Rest, cond_Stim8hr, cond_Stim48hr, gating_shape, power_rho, power_confounded, degeneracy, top_regulators, top_genes, program_label`), `docs/tables/operator_cp_stability.csv` (`rank, mean_matched_cosine`), `docs/tables/operator_cp_enrichment.csv`, `docs/figures/33_operator_cp_stability.png`, `docs/figures/34_operator_cp_condition_factors.png`.

**Nuisance controls this driver MUST implement (each maps to a spec failure mode):**
1. **Global condition scale** → `frobenius_normalize_conditions` before fitting; condition factors read as *relative* modulation. A `--scale-control none` counterfactual proves the control matters.
2. **Precision weighting** → `n_cells_target` proxy (no local `lfcSE`): `weights = n_cells / nanmax(n_cells)`. Print a WARNING that this is a power proxy, not inverse-variance.
3. **Rank by split-half stability** → sweep ranks 2..max, pick the largest rank with `mean_matched_cosine > threshold`; fall back to argmax if none clear the bar (and say so).
4. **Per-factor power test** → `spearman(a_k, n_cells_mean_per_regulator)`; `power_confounded = |rho| > 0.3`.
5. **Degeneracy** → `cp_degeneracy`; flag factors with congruence > 0.9.
6. **Init stability** → `n_init=10` inside `cp_fit_masked` (best-of).
7. **Selection-leakage cross-check** → `--sensitivity` prints the procedure to re-run on the `inclusion=any / gene-var=Rest` tensor and compare gating shapes.

- [ ] **Step 1: Write the driver**

```python
#!/usr/bin/env python3
"""Step 2 — CP decomposition of the operator tensor: regulator (x) gene (x) condition.

'TCR gated, chromatin constitutive' becomes the SHAPE of the condition factor c_k,
read off one factor — but only after the global condition-scale size effect is
removed, else every program is trivially 'gated'.

    python scripts/decompose_operator_cp.py --rank auto --scale-control frob
    python scripts/decompose_operator_cp.py --rank 4 --scale-control none
    python scripts/decompose_operator_cp.py --rank auto --sensitivity

Outputs: docs/tables/operator_cp_factors.csv, operator_cp_stability.csv,
         operator_cp_enrichment.csv, docs/figures/33_*, docs/figures/34_*
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scripts import _operator as op
from scripts.decompose_operator_svd import load_genesets

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def gating_shape(c):
    peak = COND_ORDER[int(np.argmax(c))]
    flat = (c.max() - c.min()) < 0.15 * (abs(c).max() + 1e-9)
    return "constitutive(flat)" if flat else f"gated(peak={peak})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rank", default="auto")           # "auto" or int
    ap.add_argument("--scale-control", choices=["frob", "none"], default="frob")
    ap.add_argument("--max-rank", type=int, default=8)
    ap.add_argument("--stab-threshold", type=float, default=0.7)
    ap.add_argument("--sensitivity", action="store_true")
    args = ap.parse_args()

    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    tensor, mask = d["tensor"].astype(np.float64), d["mask"]
    genes, regs = d["genes"].astype(str), d["regulators"].astype(str)
    n_cells = d["n_cells"].astype(np.float64)
    print(f"[confound] standing meter spearman(||slab||, n_cells)="
          f"{op.spearman_power(tensor, mask, n_cells):.3f}")

    if args.scale_control == "frob":
        tensor_fit, scales = op.frobenius_normalize_conditions(tensor, mask)
        print(f"[scale-control] per-condition Frobenius scales={scales.round(2)}")
    else:
        tensor_fit = tensor

    print("[WARN] precision weights use n_cells_target as a power PROXY, "
          "not inverse-variance (lfcSE not in local cache).")
    weights = n_cells / (np.nanmax(n_cells) + 1e-12)

    # --- rank selection by split-half stability ---
    stab_rows = []
    for r in range(2, args.max_rank + 1):
        s = op.split_half_stability(tensor_fit, mask, r, n_splits=6, random_state=0)
        stab_rows.append(dict(rank=r, mean_matched_cosine=s))
        print(f"[stability] rank={r} mean_matched_cosine={s:.3f}")
    stab = pd.DataFrame(stab_rows)
    if args.rank == "auto":
        ok = stab[stab["mean_matched_cosine"] > args.stab_threshold]
        rank = int(ok["rank"].max()) if len(ok) else int(
            stab.loc[stab["mean_matched_cosine"].idxmax(), "rank"])
        print(f"[rank] auto-selected rank={rank}"
              + ("" if len(ok) else " (NO rank cleared threshold; used argmax)"))
    else:
        rank = int(args.rank)

    lam, factors = op.cp_fit_masked(tensor_fit, mask, rank, n_init=10,
                                    random_state=0, weights=weights)
    lam, factors = op.fix_cp_gauge(lam, factors)
    A, B, C = factors
    degen = op.cp_degeneracy(factors)
    reg_ncell = np.nanmean(n_cells, axis=1)

    genesets, background = load_genesets(), list(genes)
    frows, erows = [], []
    for k in range(rank):
        rho = float(spearmanr(A[:, k], reg_ncell).statistic)
        top_reg = regs[np.argsort(-np.abs(A[:, k]))[:8]]
        gsel = np.argsort(-np.abs(B[:, k]))[:50]
        enr = op.hypergeometric_enrichment(list(genes[gsel]), genesets, background)
        enr.insert(0, "factor", k + 1); erows.append(enr)
        label = enr.iloc[0]["set_name"] if (len(enr) and enr.iloc[0]["fdr"] < 0.05) else "unlabeled"
        frows.append(dict(
            factor=k + 1, lambda_=float(lam[k]),
            cond_Rest=float(C[0, k]), cond_Stim8hr=float(C[1, k]),
            cond_Stim48hr=float(C[2, k]), gating_shape=gating_shape(C[:, k]),
            power_rho=rho, power_confounded=bool(abs(rho) > 0.3),
            degeneracy=float(degen[k]),
            top_regulators=";".join(map(str, top_reg)),
            top_genes=";".join(map(str, genes[gsel[:10]])),
            program_label=label))

    TAB.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(frows).to_csv(TAB / "operator_cp_factors.csv", index=False)
    stab.to_csv(TAB / "operator_cp_stability.csv", index=False)
    pd.concat(erows, ignore_index=True).to_csv(
        TAB / "operator_cp_enrichment.csv", index=False)

    FIG.mkdir(exist_ok=True, parents=True)
    f1, a1 = plt.subplots(figsize=(6, 4))
    a1.plot(stab["rank"], stab["mean_matched_cosine"], "o-")
    a1.axhline(args.stab_threshold, ls="--", c="grey")
    a1.set_xlabel("CP rank"); a1.set_ylabel("split-half matched cosine")
    a1.set_title("CP rank selection by split-half stability")
    f1.tight_layout(); f1.savefig(FIG / "33_operator_cp_stability.png", dpi=150)

    f2, a2 = plt.subplots(figsize=(7, 4))
    for k in range(rank):
        a2.plot(COND_ORDER, C[:, k], "o-",
                label=f"factor {k+1}: {frows[k]['program_label']}")
    a2.set_ylabel("condition modulation c_k"); a2.legend(fontsize=7)
    a2.set_title("CP condition factors (after Frobenius scale control)")
    f2.tight_layout(); f2.savefig(FIG / "34_operator_cp_condition_factors.png", dpi=150)

    gated = [r for r in frows if r["gating_shape"].startswith("gated")
             and not r["power_confounded"] and r["degeneracy"] < 0.9]
    flat = [r for r in frows if r["gating_shape"].startswith("constitutive")
            and not r["power_confounded"] and r["degeneracy"] < 0.9]
    print(f"[result] rank={rank}; clean gated factors={[r['factor'] for r in gated]}; "
          f"clean constitutive factors={[r['factor'] for r in flat]}")

    if args.sensitivity:
        print("[sensitivity] re-run build_operator_tensor.py with "
              "--inclusion any --gene-var-cond Rest, then this script, and compare "
              "gating_shape per matched factor. Flips => selection artifact.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python scripts/decompose_operator_cp.py --rank auto --scale-control frob`
Expected: prints the standing confound meter, per-condition Frobenius scales, a stability-vs-rank curve, an auto-selected rank (likely 3–6), and a `[result]` line listing clean gated vs constitutive factors. Writes 3 tables + 2 figures.

- [ ] **Step 3: Verify the flagship acceptance criterion**

Run:
```bash
python -c "import pandas as pd; f=pd.read_csv('docs/tables/operator_cp_factors.csv'); \
print(f[['factor','gating_shape','program_label','power_confounded','degeneracy']].to_string(index=False))"
```
Expected (acceptance): after the Frobenius scale control, **at least one clean factor** (`power_confounded=False`, `degeneracy<0.9`) has a **gated** shape (peaked on a Stim condition) with a TCR/immune `program_label`, **and** at least one clean factor is **constitutive(flat)** with a chromatin/interferon label. If the only gated factors are power-confounded or degenerate, the gating headline is an artifact — record that as the result (a real negative that saves shipping a false claim).

- [ ] **Step 4: Run the scale-control counterfactual (proves the control is load-bearing)**

Run: `python scripts/decompose_operator_cp.py --rank 4 --scale-control none`
Expected: with no scale control, condition factors lean toward the higher-DE Stim conditions for most/all factors (gating trivially "true" everywhere), showing why the `frob` control matters. Note the contrast in the writeup; **then re-run Step 2's default before committing** so the committed tables are the `frob` ones.

- [ ] **Step 5: Add Makefile target and commit**

Add `operator-cp` to `.PHONY` and:
```makefile
operator-cp:
	$(PY) scripts/decompose_operator_cp.py --rank auto --scale-control frob
```

```bash
git add scripts/decompose_operator_cp.py Makefile docs/tables/operator_cp_*.csv \
        docs/figures/33_operator_cp_stability.png docs/figures/34_operator_cp_condition_factors.png
git commit -m "feat(operator): Step 2 driver — CP with scale control, stability rank, confound gate"
```

---

## Task 7: Kernel — soft-impute completion + train-only standardization

**Files:**
- Modify: `scripts/_operator.py` (add `train_test_standardize`, `soft_impute`)
- Modify: `tests/test_operator.py`
- Test: `tests/test_operator.py`

**Interfaces:**
- Produces:
  - `train_test_standardize(M, train_mask) -> (Mc, mu)` — per-column mean computed on **train entries only** (`train_mask` True = train); returns centered matrix (test entries centered by the train mean) and `mu (G,)`. Never uses test entries to fit `mu`.
  - `soft_impute(M, observed_mask, rank, n_iter=100, tol=1e-4) -> M_hat` — iterative SVD hard-truncation-to-`rank` imputation: init unobserved to 0, repeat {SVD, truncate to `rank`, reset observed entries to true values} until Frobenius change < `tol`.

- [ ] **Step 1: Write the failing tests**

```python
def test_train_test_standardize_uses_train_only():
    M = np.array([[1.0, 10.0], [3.0, 30.0], [100.0, 100.0]])
    train = np.array([[True, True], [True, True], [False, False]])  # row 2 held out
    Mc, mu = op.train_test_standardize(M, train)
    assert np.allclose(mu, [2.0, 20.0])          # mean of train rows only
    assert np.allclose(Mc[2], [98.0, 80.0])      # test row centered by TRAIN mean


def test_soft_impute_recovers_low_rank():
    rng = np.random.default_rng(0)
    U = rng.normal(size=(30, 2)); V = rng.normal(size=(2, 20))
    M = U @ V
    obs = rng.random(M.shape) > 0.3              # 70% observed
    Mhat = op.soft_impute(M, obs, rank=2, n_iter=200)
    held = ~obs
    r2 = 1 - np.sum((M[held] - Mhat[held]) ** 2) / np.sum(
        (M[held] - M[obs].mean()) ** 2)
    assert r2 > 0.9
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_operator.py -q`
Expected: FAIL — `AttributeError: ... 'train_test_standardize'`.

- [ ] **Step 3: Implement**

Append to `scripts/_operator.py`:
```python
def train_test_standardize(M, train_mask):
    M = np.asarray(M, float)
    mu = np.zeros(M.shape[1])
    for j in range(M.shape[1]):
        col = M[train_mask[:, j], j]
        mu[j] = col.mean() if col.size else 0.0
    return M - mu[None, :], mu


def soft_impute(M, observed_mask, rank, n_iter=100, tol=1e-4):
    M = np.asarray(M, float)
    X = np.where(observed_mask, M, 0.0)
    Xr = X
    prev = np.inf
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

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_operator.py -q`
Expected: PASS (11 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/_operator.py tests/test_operator.py
git commit -m "feat(operator): soft-impute completion + train-only standardization kernels"
```

---

## Task 8: Step 3 driver — held-out prediction (entry-wise + condition extrapolation)

**Files:**
- Create: `scripts/operator_completion.py`
- Modify: `Makefile` (add `operator-completion`)
- Test: run on real data; acceptance = 3a beats per-gene-mean, 3b beats persistence.

**Interfaces:**
- Consumes: `data/cache/operator_tensor.npz`; `_operator.train_test_standardize`, `soft_impute`.
- Produces: `docs/tables/operator_completion_entrywise.csv` (`rank, r2, corr, mse_model, mse_gene_mean, beats_baseline`), `docs/tables/operator_completion_condition.csv` (`rank, r2_model, r2_persistence, beats_persistence, r2_lowSE, r2_highSE`), `docs/figures/35_operator_completion_curve.png`.

**Nuisance controls this driver MUST implement:**
- **Standardization leakage** → all centering via `train_test_standardize` on train entries only; gene selection inherited from the cached tensor (fixed, not re-selected on test).
- **Non-trivial baseline** → 3a vs per-gene-**mean** (not zero); 3b vs **persistence** (Stim48hr = Stim8hr).
- **Row-cold-start honesty** → the driver prints that pure low-rank cannot predict a regulator with zero observed entries; 3b holds out a condition *fiber*, never a whole regulator.
- **Target-noise stratification** → 3b reports r2 split by held-out `n_cells_target` median (low vs high power) as the `lfcSE` proxy.

- [ ] **Step 1: Write the driver**

```python
#!/usr/bin/env python3
"""Step 3 — is the operator recoverably low-rank? Held-out prediction.

3a entry-wise completion: mask 20% of observed entries, soft-impute at ranks
1..K, predict them; bar = beat the per-gene-mean baseline.
3b condition extrapolation: hold out entire (regulator, Stim48hr) fibers, predict
from that regulator's Rest+Stim8hr via the low-rank fit on OTHER regulators;
bar = beat persistence (Stim48hr := Stim8hr).

    python scripts/operator_completion.py --max-rank 12 --holdout 0.2

Outputs: docs/tables/operator_completion_entrywise.csv,
         operator_completion_condition.csv, docs/figures/35_*
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts import _operator as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def to_matrix(d):
    """Flatten fully-observed fingerprint rows to a (rows, G) matrix + index."""
    tensor, mask = d["tensor"], d["mask"]
    R, G, C = tensor.shape
    rows, idx = [], []
    for i in range(R):
        for c in range(C):
            if mask[i, :, c].all():
                rows.append(tensor[i, :, c]); idx.append((i, c))
    return np.asarray(rows, float), idx


def entrywise(M, max_rank, holdout, seed):
    rng = np.random.default_rng(seed)
    hide = rng.random(M.shape) < holdout
    train = ~hide
    Mc, mu = op.train_test_standardize(M, train)
    rows = []
    for r in range(1, max_rank + 1):
        Mhat = op.soft_impute(Mc, train, rank=r, n_iter=200)
        yt, yp = Mc[hide], Mhat[hide]
        mse_model = float(np.mean((yt - yp) ** 2))
        mse_base = float(np.mean(yt ** 2))       # gene-mean == 0 in centered space
        r2 = 1 - mse_model / (mse_base + 1e-12)
        corr = float(np.corrcoef(yt, yp)[0, 1])
        rows.append(dict(rank=r, r2=r2, corr=corr, mse_model=mse_model,
                         mse_gene_mean=mse_base, beats_baseline=bool(mse_model < mse_base)))
    return pd.DataFrame(rows)


def condition_extrap(d, max_rank, seed):
    tensor, mask, n_cells = d["tensor"], d["mask"], d["n_cells"]
    R, G, C = tensor.shape
    late, early = COND_ORDER.index("Stim48hr"), COND_ORDER.index("Stim8hr")
    full = np.array([i for i in range(R) if mask[i, :, :].all(axis=0).all()])
    rng = np.random.default_rng(seed)
    test = rng.permutation(full)[: max(10, len(full) // 5)]
    X = np.concatenate([tensor[full, :, 0], tensor[full, :, early],
                        tensor[full, :, late]], axis=1)         # Rest|Stim8|Stim48
    obs = np.ones_like(X, dtype=bool)
    test_pos = np.isin(full, test)
    obs[test_pos, 2 * G:3 * G] = False                          # hide test late fibers
    Xc, mu = op.train_test_standardize(X, obs)
    med = np.nanmedian(n_cells[test, late])
    rows = []
    for r in range(1, max_rank + 1):
        Xhat = op.soft_impute(Xc, obs, rank=r, n_iter=200)
        yt = Xc[np.ix_(test_pos, np.arange(2 * G, 3 * G))]
        yp = Xhat[np.ix_(test_pos, np.arange(2 * G, 3 * G))]
        pers = Xc[np.ix_(test_pos, np.arange(G, 2 * G))]        # persistence baseline
        ss = np.sum((yt - yt.mean()) ** 2) + 1e-12
        r2_model = float(1 - np.sum((yt - yp) ** 2) / ss)
        r2_pers = float(1 - np.sum((yt - pers) ** 2) / ss)
        lo = n_cells[test, late] <= med; hi = ~lo
        def _r2(sub):
            a, b = yt[sub], yp[sub]
            return float(1 - np.sum((a - b) ** 2) / (np.sum((a - a.mean()) ** 2) + 1e-12))
        rows.append(dict(rank=r, r2_model=r2_model, r2_persistence=r2_pers,
                         beats_persistence=bool(r2_model > r2_pers),
                         r2_lowSE=_r2(lo), r2_highSE=_r2(hi)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-rank", type=int, default=12)
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    print("[note] pure low-rank CANNOT predict a regulator with zero observed "
          "entries (no anchor). 3b holds out a condition FIBER, not a whole regulator.")
    M, _ = to_matrix(d)
    ew = entrywise(M, args.max_rank, args.holdout, args.seed)
    ce = condition_extrap(d, args.max_rank, args.seed)
    TAB.mkdir(exist_ok=True, parents=True)
    ew.to_csv(TAB / "operator_completion_entrywise.csv", index=False)
    ce.to_csv(TAB / "operator_completion_condition.csv", index=False)

    FIG.mkdir(exist_ok=True, parents=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ew["rank"], ew["r2"], "o-", label="3a entry-wise R² vs gene-mean")
    ax.plot(ce["rank"], ce["r2_model"], "s-", label="3b Stim48hr model R²")
    ax.plot(ce["rank"], ce["r2_persistence"], "--", c="crimson", label="3b persistence")
    ax.axhline(0, ls=":", c="grey")
    ax.set_xlabel("rank"); ax.set_ylabel("held-out R²"); ax.legend(fontsize=8)
    ax.set_title("Operator low-rank recoverability")
    fig.tight_layout(); fig.savefig(FIG / "35_operator_completion_curve.png", dpi=150)

    print(f"[3a] best rank={int(ew.loc[ew.r2.idxmax(),'rank'])} R²={ew.r2.max():.3f} "
          f"beats gene-mean at ranks={ew.loc[ew.beats_baseline,'rank'].tolist()}")
    print(f"[3b] beats persistence at ranks={ce.loc[ce.beats_persistence,'rank'].tolist()}; "
          f"max model R²={ce.r2_model.max():.3f} vs persistence={ce.r2_persistence.max():.3f}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it on real data**

Run: `python scripts/operator_completion.py --max-rank 12 --holdout 0.2`
Expected: writes 2 tables + 1 figure; prints the best entry-wise rank/R² and the ranks where 3a beats gene-mean and 3b beats persistence.

- [ ] **Step 3: Verify acceptance**

Run:
```bash
python -c "import pandas as pd; a=pd.read_csv('docs/tables/operator_completion_entrywise.csv'); \
b=pd.read_csv('docs/tables/operator_completion_condition.csv'); \
print('3a beats gene-mean:', bool(a.beats_baseline.any()), 'best R2', round(a.r2.max(),3)); \
print('3b beats persistence:', bool(b.beats_persistence.any()), 'best model R2', round(b.r2_model.max(),3))"
```
Expected (acceptance): **3a beats per-gene-mean** MSE at the selected rank with a visible elbow (→ operator effectively rank-k), **and** **3b beats persistence** on held-out Stim48hr. Either one cleanly converts the submission from descriptive to predictive; both is the strong result. If 3b does not beat persistence, report that late-stim is well-approximated by early-stim (still an honest finding).

- [ ] **Step 4: Add Makefile target and commit**

Add `operator-completion` to `.PHONY` and:
```makefile
operator-completion:
	$(PY) scripts/operator_completion.py --max-rank 12 --holdout 0.2
```

```bash
git add scripts/operator_completion.py Makefile docs/tables/operator_completion_*.csv \
        docs/figures/35_operator_completion_curve.png
git commit -m "feat(operator): Step 3 driver — held-out entry-wise + condition-mode prediction"
```

---

## Task 9: Kernel — principal angles + random-subspace null

**Files:**
- Modify: `scripts/_operator.py` (add `principal_angles`, `random_subspace_null`)
- Modify: `tests/test_operator.py`
- Test: `tests/test_operator.py`

**Interfaces:**
- Produces:
  - `principal_angles(Va, Vb) -> cos2` — for two `(G, k)` bases (not necessarily orthonormal), QR-orthonormalize each and return `cos²θ_i` sorted descending, length `k`. Uses `scipy.linalg.subspace_angles`.
  - `random_subspace_null(G, k, n=1000, random_state=0) -> (mean_cos2, p95_cos2)` — mean and 95th percentile of mean `cos²θ` between two random `k`-subspaces of `R^G`.

- [ ] **Step 1: Write the failing tests**

```python
def test_principal_angles_identical_subspace():
    rng = np.random.default_rng(0)
    V = rng.normal(size=(50, 4))
    cos2 = op.principal_angles(V, V.copy())
    assert np.allclose(cos2, 1.0, atol=1e-6)


def test_principal_angles_orthogonal_subspace():
    G = 20
    Va = np.zeros((G, 2)); Va[0, 0] = Va[1, 1] = 1.0
    Vb = np.zeros((G, 2)); Vb[2, 0] = Vb[3, 1] = 1.0
    cos2 = op.principal_angles(Va, Vb)
    assert np.allclose(cos2, 0.0, atol=1e-6)


def test_random_subspace_null_small_for_high_dim():
    mean_c2, p95 = op.random_subspace_null(G=200, k=4, n=200, random_state=0)
    assert mean_c2 < 0.1          # random k-subspaces in R^200 barely overlap
    assert p95 >= mean_c2
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `python -m pytest tests/test_operator.py -q`
Expected: FAIL — `AttributeError: ... 'principal_angles'`.

- [ ] **Step 3: Implement**

Append to `scripts/_operator.py`:
```python
from scipy.linalg import subspace_angles, qr


def principal_angles(Va, Vb):
    Qa = qr(np.asarray(Va, float), mode="economic")[0]
    Qb = qr(np.asarray(Vb, float), mode="economic")[0]
    ang = subspace_angles(Qa, Qb)
    cos2 = np.cos(ang) ** 2
    return np.sort(cos2)[::-1]


def random_subspace_null(G, k, n=1000, random_state=0):
    rng = np.random.default_rng(random_state)
    means = np.empty(n)
    for t in range(n):
        A = rng.normal(size=(G, k)); B = rng.normal(size=(G, k))
        means[t] = principal_angles(A, B).mean()
    return float(means.mean()), float(np.percentile(means, 95))
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `python -m pytest tests/test_operator.py -q`
Expected: PASS (14 passed).

- [ ] **Step 5: Commit**

```bash
git add scripts/_operator.py tests/test_operator.py
git commit -m "feat(operator): principal angles + random-subspace null kernels"
```

---

## Task 10: Step 4 driver — donor-subspace stability across disjoint donor pairs

**Files:**
- Create: `scripts/operator_donor_angles.py`
- Modify: `Makefile` (add `operator-donors`)
- Test: run on real data **if** per-donor matrices are available; else the driver reports the missing-data condition and exits 0.

**Interfaces:**
- Consumes: per-donor-pair log-FC matrices. **These are NOT in the local cache** (`data/cache/` has only `donor_obs.csv`, a per-regulator summary — verified). The driver supports, in priority order:
  1. `data/cache/by_donors_index.csv` (columns `regulator,gene`) + `data/cache/by_donors_pair_<ab>.npy` (`R x G`) if a fetch step (analogous to `scripts/fetch_fingerprint_matrix.py`) has produced them;
  2. absent → print a clear NEEDS-DATA message naming the expected files + fetch pattern, write an empty-but-headed table, and exit 0 (do not crash the pipeline).
- Produces: `docs/tables/operator_donor_angles.csv` (`pair_a, pair_b, k, mean_cos2, null_mean_cos2, null_p95, above_null`), `docs/figures/36_operator_donor_angles.png` (only when data present).

**Nuisance controls this driver MUST implement:**
- **Donor-pair non-independence** → compare **disjoint** pairs only. With donors {1,2,3,4}, the three disjoint comparisons are `(1,2)vs(3,4)`, `(1,3)vs(2,4)`, `(1,4)vs(2,3)`. Hardcode these three; never emit the 15 inflated overlapping comparisons.
- **Identical gene/regulator axis** → assume the shared `by_donors_index.csv` defines a common axis for every pair matrix (the fetch step guarantees it).
- **Subspace, not per-vector** → use `principal_angles` only; never correlate `v_k` pairwise.

- [ ] **Step 1: Write the driver**

```python
#!/usr/bin/env python3
"""Step 4 — are the gene programs donor-reproducible AS SUBSPACES?

Compare the top-k gene-program subspace of one donor-pair modality against another
via principal angles, but ONLY across DISJOINT donor pairs (the 6 modalities share
donors, so overlapping pairs inflate overlap by construction). With donors {1,2,3,4}
the honest comparisons are (1,2)vs(3,4), (1,3)vs(2,4), (1,4)vs(2,3).

Requires per-donor-pair matrices (NOT in the default local cache). If absent, prints
the expected files + fetch pattern and exits cleanly.

    python scripts/operator_donor_angles.py --k 5

Outputs: docs/tables/operator_donor_angles.csv, docs/figures/36_*
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts import _operator as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
DISJOINT = [((1, 2), (3, 4)), ((1, 3), (2, 4)), ((1, 4), (2, 3))]


def load_pair_matrices():
    """Return {pair_tuple: M (R,G)} using the shared axis, or None if unavailable."""
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
    _, _, Vt = np.linalg.svd(Mc, full_matrices=False)
    return Vt[:k].T                    # (G, k)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--null-n", type=int, default=500)
    args = ap.parse_args()

    TAB.mkdir(exist_ok=True, parents=True)
    cols = ["pair_a", "pair_b", "k", "mean_cos2", "null_mean_cos2", "null_p95", "above_null"]
    mats = load_pair_matrices()
    if mats is None:
        print("[NEEDS-DATA] per-donor-pair matrices absent. Expected:\n"
              "  data/cache/by_donors_index.csv (columns: regulator,gene) AND\n"
              "  data/cache/by_donors_pair_12.npy ... by_donors_pair_34.npy  (R x G)\n"
              "Produce them with a fetch analogous to scripts/fetch_fingerprint_matrix.py\n"
              "reading the per-donor layers of the remote h5ad. Writing empty table, exit 0.")
        pd.DataFrame(columns=cols).to_csv(TAB / "operator_donor_angles.csv", index=False)
        return

    G = next(iter(mats.values())).shape[1]
    null_mean, null_p95 = op.random_subspace_null(G, args.k, n=args.null_n)
    rows = []
    for (a, b) in DISJOINT:
        if a not in mats or b not in mats:
            continue
        Va = top_subspace(mats[a], args.k)
        Vb = top_subspace(mats[b], args.k)
        c2 = op.principal_angles(Va, Vb)
        rows.append(dict(pair_a=f"{a[0]}{a[1]}", pair_b=f"{b[0]}{b[1]}", k=args.k,
                         mean_cos2=float(c2.mean()), null_mean_cos2=null_mean,
                         null_p95=null_p95, above_null=bool(c2.mean() > null_p95)))
    df = pd.DataFrame(rows, columns=cols)
    df.to_csv(TAB / "operator_donor_angles.csv", index=False)

    FIG.mkdir(exist_ok=True, parents=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(df["pair_a"] + " vs " + df["pair_b"], df["mean_cos2"])
    ax.axhline(null_p95, ls="--", c="crimson", label="random-subspace 95th pct")
    ax.set_ylabel(f"mean cos²θ (top-{args.k})"); ax.legend()
    ax.set_title("Program-subspace overlap across disjoint donor pairs")
    fig.tight_layout(); fig.savefig(FIG / "36_operator_donor_angles.png", dpi=150)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

Run: `python scripts/operator_donor_angles.py --k 5`
Expected: since per-donor matrices are not in the local cache, prints the `[NEEDS-DATA]` message, writes an empty-headed `operator_donor_angles.csv`, and exits 0. (Correct behavior for the current cache — do not fake data to force a number.)

- [ ] **Step 3: (Conditional) verify the real result when data is present**

Only if `data/cache/by_donors_pair_*.npy` have been produced: re-run and confirm `above_null=True` for all three disjoint comparisons. Acceptance (when data available): mean `cos²θ` well above the random-subspace 95th percentile → programs are donor-reproducible as subspaces.

- [ ] **Step 4: Add Makefile target and commit**

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

## Task 11: Integration — umbrella Makefile target, writeup, README pointer

**Files:**
- Modify: `Makefile` (umbrella `operator` target + help line)
- Create: `docs/OPERATOR_ANALYSIS.md`
- Modify: `README.md` (one pointer line)
- Test: `make operator` runs Steps 0–3 end to end; the writeup reflects the actual committed numbers.

**Interfaces:**
- Consumes: all tables/figures produced by Tasks 2, 4, 6, 8, 10.

- [ ] **Step 1: Add the umbrella target and help line**

Add `operator` to `.PHONY`; add to the `help:` block:
```makefile
	@echo "  make operator - empirical regulatory operator: tensor + SVD + CP + completion [+donors if fetched]"
```
Add the target (Steps 0–3 fully local; donors best-effort):
```makefile
operator: operator-tensor operator-svd operator-cp operator-completion operator-donors
```

- [ ] **Step 2: Run the full local pipeline**

Run: `make operator`
Expected: Steps 0→1→2→3 complete and write their tables/figures; Step 4 prints `[NEEDS-DATA]` and exits 0 (does not fail the target).

- [ ] **Step 3: Run the kernel test suite once more (full green before writeup)**

Run: `python -m pytest tests/test_operator.py -q`
Expected: PASS (14 passed).

- [ ] **Step 4: Write `docs/OPERATOR_ANALYSIS.md` from the actual numbers**

Read the committed tables and write the doc. Required sections, each stating the **real** value produced (no placeholders):
```markdown
# The Empirical Regulatory Operator

One matrix, four questions. L[reg, gene] = log-FC of every gene under every
regulator knockdown — the object the fingerprint PCA computed and half-discarded.

## The object and the tensor (Step 0)
- Tensor shape, inclusion rule, gene selection, and the standing confound meter
  spearman(||slab||, n_cells) = <value from operator_tensor_summary.json>.

## Gene programs (Step 1)
- Top-5 programs and their FDR-significant labels (operator_svd_enrichment.csv),
  raw vs varimax. State how many of the top-5 earned a clean label (acceptance ≥2).

## Condition gating (Step 2) — the flagship
- Stability-selected rank; condition-factor shapes AFTER the Frobenius scale
  control; which factors are clean (not power-confounded, not degenerate); the
  gated-vs-constitutive result. Include the scale-control counterfactual: without
  it, gating is trivially true for all factors.

## Prediction (Step 3)
- 3a: does low-rank beat per-gene-mean, and the effective-rank elbow.
- 3b: does low-rank beat persistence on held-out Stim48hr; note the row-cold-start
  limitation (cannot impute a never-paneled regulator from factorization alone).

## Donor reproducibility of programs (Step 4)
- Either the disjoint-pair principal-angle result, or the explicit NEEDS-DATA note
  (per-donor matrices not in the local cache) + how to fetch them.

## What the nuisance controls bought us
- One line per control (power, condition-scale, selection leakage, standardization
  leakage, gauge, rank stability, donor-pair independence, baselines) — for each,
  say whether a naive version would have produced a different (wrong) headline.
```

- [ ] **Step 5: Add a README pointer and commit**

Add one line under the analyses section of `README.md`:
```markdown
- **Empirical regulatory operator** (`make operator`) — tensor/CP/completion/donor-subspace analysis of the log-FC operator; see `docs/OPERATOR_ANALYSIS.md`.
```

```bash
git add Makefile docs/OPERATOR_ANALYSIS.md README.md
git commit -m "docs(operator): umbrella make target, analysis writeup, README pointer"
```

---

## Task 12 (STRETCH — only if Tasks 6 and 8 land cleanly): network deconvolution + asymmetric subsumption

**Files:**
- Create: `scripts/operator_deconvolution.py`
- Modify: `Makefile` (add `operator-deconv`)
- Test: run on real data; acceptance = recovered asymmetry presented as hypothesis-generating, with explicit condition-number reporting.

**Interfaces:**
- Consumes: `data/cache/operator_tensor.npz`; `docs/tables/top_robust_regulators.csv` (existing) for the donor-robust shortlist.
- Produces: `docs/tables/operator_deconvolution_edges.csv` (`source, target, weight, block_condition, cond_number`), `docs/tables/operator_subsumption.csv` (`source, target, condition, subsumption_frac`).

**Nuisance controls (spec's stretch caveats, mandatory):**
- Deconvolution `A ≈ I − L⁻¹` is meaningful only on the **square block** where perturbed regulators are themselves readout genes; restrict to that intersection, report its size, and **regularize** (ridge). Report the block **condition number** per condition — an ill-conditioned inverse invalidates the edges.
- It is a **linear approximation to a nonlinear system** — asymmetry is hypothesis-generating, never ground truth. Filename and doc must say so.
- Subsumption `s(i→j)` = fraction of `j`'s significant set contained in `i`'s; directional, expresses hierarchy cosine cannot.

- [ ] **Step 1: Write the driver**

```python
#!/usr/bin/env python3
"""Step 5 (stretch) — square-block network deconvolution + asymmetric subsumption.

HYPOTHESIS-GENERATING ONLY. A ≈ I - L^{-1} is valid only on the square block where
perturbed regulators are also readout genes, is a linear approx to a nonlinear
system, and is ill-conditioned — every edge ships with the block condition number.
Subsumption expresses hierarchy that cosine cannot.

    python scripts/operator_deconvolution.py --n-robust 50 --ridge 1e-2

Outputs: docs/tables/operator_deconvolution_edges.csv, operator_subsumption.csv
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def deconvolve(L, ridge):
    n = L.shape[0]
    cond = float(np.linalg.cond(L + ridge * np.eye(n)))
    A = np.eye(n) - np.linalg.inv(L + ridge * np.eye(n))
    return A, cond


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-robust", type=int, default=50)
    ap.add_argument("--ridge", type=float, default=1e-2)
    ap.add_argument("--sig-thr", type=float, default=0.5)
    args = ap.parse_args()

    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    tensor = d["tensor"]
    regs, genes = d["regulators"].astype(str), d["genes"].astype(str)
    robust_path = TAB / "top_robust_regulators.csv"
    shortlist = (pd.read_csv(robust_path).iloc[:, 0].astype(str).tolist()
                 if robust_path.exists() else list(regs))
    reg_set, gene_set = set(regs), set(genes)
    block = [g for g in shortlist if g in reg_set and g in gene_set][: args.n_robust]
    if len(block) < 10:
        print(f"[deconv] square block too small ({len(block)}); most regulators are "
              "not in the readout gene panel. Deconvolution NOT supported by the data.")
    ri = [list(regs).index(g) for g in block]
    gi = [list(genes).index(g) for g in block]
    print(f"[deconv] square block size = {len(block)} regulators")

    edges = []
    for c, cond in enumerate(COND_ORDER):
        L = np.nan_to_num(tensor[np.ix_(ri, gi, [c])][:, :, 0], nan=0.0)
        A, condnum = deconvolve(L, args.ridge)
        for i, s in enumerate(block):
            for j, t in enumerate(block):
                if i != j and A[i, j] != 0:
                    edges.append(dict(source=s, target=t, weight=float(A[i, j]),
                                      block_condition=cond, cond_number=condnum))
    pd.DataFrame(edges).to_csv(TAB / "operator_deconvolution_edges.csv", index=False)

    sub = []
    for c in range(3):
        sets = {g: set(np.where(np.abs(np.nan_to_num(tensor[i, :, c])) > args.sig_thr)[0])
                for i, g in enumerate(regs)}
        for i in block:
            for j in block:
                if i != j and len(sets[j]):
                    frac = len(sets[i] & sets[j]) / len(sets[j])
                    sub.append(dict(source=i, target=j, condition=COND_ORDER[c],
                                    subsumption_frac=float(frac)))
    pd.DataFrame(sub).to_csv(TAB / "operator_subsumption.csv", index=False)
    print(f"[deconv] wrote {len(edges)} edges; condition numbers reported per condition.")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it**

Run: `python scripts/operator_deconvolution.py --n-robust 50 --ridge 1e-2`
Expected: prints the square-block size and writes the two tables. If the block is < ~10, it reports that deconvolution is not supported by the data rather than forcing edges.

- [ ] **Step 3: Verify the condition-number caveat is honored**

Run:
```bash
python -c "import pandas as pd; e=pd.read_csv('docs/tables/operator_deconvolution_edges.csv'); \
print(e.groupby('block_condition')['cond_number'].first() if len(e) else 'no edges (block too small)')"
```
Expected: every edge carries the block condition number; if condition numbers are huge (≫1e6), the writeup must flag the inverse as unreliable and present edges as speculative only.

- [ ] **Step 4: Add Makefile target and commit**

Add `operator-deconv` to `.PHONY` and:
```makefile
operator-deconv:
	$(PY) scripts/operator_deconvolution.py --n-robust 50 --ridge 1e-2
```

```bash
git add scripts/operator_deconvolution.py Makefile \
        docs/tables/operator_deconvolution_edges.csv docs/tables/operator_subsumption.csv
git commit -m "feat(operator): Step 5 stretch — square-block deconvolution + asymmetric subsumption"
```

---

## Global nuisance checklist → where each is enforced

| # | Confound | Enforced in |
|---|----------|-------------|
| 1 | Power / cell-count | Step 0 standing meter (Task 2); per-factor `power_rho`/`power_confounded` gate + precision-proxy weights (Task 6) |
| 2 | Global condition scale vs gating | `frobenius_normalize_conditions` before CP + `--scale-control none` counterfactual (Task 6) |
| 3 | Selection leakage | `--inclusion any --gene-var-cond Rest` alt build (Task 2) + `--sensitivity` flip check (Task 6) |
| 4 | Standardization leakage | `train_test_standardize` fit on train entries only (Tasks 7, 8) |
| 5 | Sign/scale/rotation gauge | `fix_cp_gauge` (Task 5); `varimax` for interpretable programs (Tasks 3, 4); principal angles for subspaces (Tasks 9, 10) |
| 6 | Rank selection | `split_half_stability` sweep, no scree-eyeballing (Tasks 5, 6); `n_init=10` init stability |
| 7 | Donor-pair non-independence | disjoint pairs only (Task 10) |
| 8 | Baselines | per-gene-mean (3a) + persistence (3b), never zero (Task 8) |
| 9 | Multiple testing | BH-FDR within each enrichment family (Task 3 kernel) |
| 10 | Orthogonality artifact | raw vs varimax reported; orthogonality never called biology (Tasks 3, 4, writeup) |

---

## Self-Review (run against the spec)

**Spec coverage:** Step 0 → Tasks 1–2. Step 1 → Tasks 3–4. Step 2 (flagship CP) → Tasks 5–6. Step 3a/3b → Tasks 7–8. Step 4 → Tasks 9–10. Step 5 stretch → Task 12. Global nuisance checklist → mapped above and enforced per-task. Integration/writeup → Task 11. No spec section is unmapped.

**Known data-reality departures from the spec (called out, not hidden):**
- The spec assumes `lfcSE` for precision weighting; it is **not in the local cache**. Task 6 uses `n_cells_target` as an explicit power *proxy* and prints a WARNING. Upgrading to true inverse-variance weights requires fetching per-donor SE (a future step).
- The spec's Step 4 references `by_donors.h5mu`; only `donor_obs.csv` (a summary) is local. Task 10 degrades gracefully (NEEDS-DATA path) and documents the exact files + fetch pattern required to produce the real result.

**Placeholder scan:** every code step contains complete, runnable code; every run step states the exact command and expected output; every acceptance step is a concrete, checkable criterion. No TBD/TODO/"handle edge cases"/"similar to Task N".

**Type consistency:** kernel signatures in the Interfaces blocks match their call sites — `cp_fit_masked` returns `(lam, factors)` with `factors = [A, B, C]`; `fix_cp_gauge` preserves shapes; `match_factors`/`split_half_stability` use the gene mode `factors[1]`; `principal_angles`/`random_subspace_null`/`train_test_standardize`/`soft_impute`/`varimax`/`hypergeometric_enrichment` signatures match Tasks 4/6/8/10 usage. `COND_ORDER` is identical everywhere. The npz keys written in Task 2 are exactly those read in Tasks 4, 6, 8, 10, 12.
