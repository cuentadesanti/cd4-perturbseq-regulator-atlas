# Modeling — from the DE matrix to regulators with uncertainty

Two small, **dependency-light** models (only `scipy` + `statsmodels`) that separate signal
from noise with uncertainty instead of ranking by raw counts and `adj_p_value < 0.1`.

> **Honest naming:** both are **empirical-Bayes / pseudo-Bayesian**. There is no PPL,
> no formal random effects, no jointly MCMC-sampled posterior. Where we say
> "posterior" we mean the normal EB approximation with prior parameters estimated from the data.

---

## Model 2 — regulator ranking (core, runs locally)

**Script:** `scripts/model_hubs.py` · **Runs on:** `DE_stats.suppl_table.csv` (local, no downloads).
**Grain:** 1 row = perturbed gene × condition.

### Specification

1. **Fixed effects (conditional mean).** GLM on `n_downstream`:

   ```
   n_downstream ~ C(culture_condition) + ontarget_significant + offtarget_flag
   ```

   Poisson and NB share the same mean model; since we only use the fitted mean `μᵢ`
   (no inference on coefficients) we fit by stable IRLS: Poisson GLM → NB `α` by
   method of moments (`Var = μ + α·μ²`) → NB GLM with fixed `α`. This avoids the
   convergence problems of full NB MLE.

2. **Empirical-Bayes shrinkage of the per-gene effect.** Log-rate deviation from the baseline:

   ```
   workᵢ = log(yᵢ + 0.5) − log(μᵢ + 0.5)
   ```

   Per gene g:  `d_g = mean(work)`,  `s²_g = σ²_e / n_g`.
   Prior `u_g ~ Normal(0, τ²)` with `τ²` by method of moments (`Var(d_g) − mean(s²_g)`).
   Approximate posterior:

   ```
   u_g | data ~ Normal( shrink·d_g ,  shrink·s²_g ),   shrink = τ²/(τ²+s²_g)
   ```

   Genes with few conditions / little signal are shrunk toward 0.

### Outputs

`docs/tables/hub_ranking_bayes.csv` (all genes) and `docs/tables/top_regulators_for_review.csv`
(top 30, judge-facing). Key columns: `regpower_eb_mean/sd` (log-rate regulatory power),
`p_top_1pct` (EB probability of exceeding the empirical top-1% threshold, not "P of being in the top 1%"),
`expected_downstream`. Figure `07_hub_posterior_ranking.png`.

### Reading the result

The robust ranking surfaces **chromatin/transcription** machinery consistent across conditions
— SAGA complex (TADA1/TADA2B/SGF29/SUPT20H/TAF6L), Mediator (MED12/CCNC), KDM1A, SETD2, CTBP1 —
above the raw TCR-signaling hubs that were Stim8hr-specific. In other words: the shrinkage
rewards regulators with a large **and** stable effect.

### Caveats

- `xcond_reproducibility` is an **exploratory feature** (cross-condition stability). It does **not**
  replace cross-donor / cross-guide reproducibility, which requires `DE_stats.h5ad`.
- The fixed-effects baseline is treated as known (plug-in) → pseudo-Bayesian, not full Bayes.
- `single_guide_estimate` and `n_guides` are NOT in the CSV; in the core review table they appear as
  `NA (requires DE_stats.h5ad)` — they are present in the sensitivity audit below.

### Guide/donor-aware sensitivity audit (optional)

When `de_obs_reproducibility_metadata.csv` exists (extracted from the `.obs` of `DE_stats.h5ad`,
without `.layers`), `model_hubs.py` runs an audit: it **reweights** the EB score with real
reproducibility (`reweighted_score = regpower_eb_mean · repro_weight`) and reports which regulators
survive (`reproducibility_audit.csv`, fig 19).

- **It is a sensitivity analysis, not a new posterior**: the EB model is NOT re-estimated.
- **Partial coverage**: `guide_correlation_all` ~78% of contrasts, `donor_correlation_hits_mean`
  only ~19% → in practice more *guide-aware* than *donor-aware*. Where the metric is missing, a
  **neutral weight** (0.75) is used, so **a gene is not penalized just for lacking donor metadata**.
- The **core** ranking does **not depend** on this file (`make all` runs without it).

---

## Model 1 — uncertainty-aware effect network (STRICTLY OPTIONAL)

**Scripts:** `scripts/model_edges_spike.py` (validation) and `scripts/model_edges.py` (scaling).
**Rule:** if the remote spike fails or is slow, the official deliverable is Model 2 + docs.

### Idea

Exact normal-normal EB on `log_fc` / `lfcSE` from the h5ad `.layers`:

```
yᵢ | θᵢ ~ Normal(θᵢ, seᵢ²)          # observed
θᵢ     ~ Normal(0, τ²)              # shrinkage prior
θᵢ|yᵢ  ~ Normal(mᵢ, vᵢ),  vᵢ = 1/(1/τ² + 1/seᵢ²),  mᵢ = vᵢ·yᵢ/seᵢ²
```

Per-edge outputs: `theta_post_mean/sd`, `p_effect_positive`, `p_abs_effect_gt_1p5x`.
**Decision rule** (more interpretable than FDR): `p_abs_effect_gt_1p5x > 0.8 AND ontarget_significant`.

### Memory/compute-aware strategy

Disk: **9.8 GB free < 17 GB** for the h5ad → not downloaded. Instead:
- Only the edges of the **candidate regulators** (top of Model 2) are needed, not the ~350M.
- `model_edges_spike.py` **measures** (does not assume) the layout/chunking and the real cost of
  reading one row per slice from S3 (anonymous `fsspec` + `h5py`). If viable, `model_edges.py`
  fetches only those rows and runs the vectorized EB (seconds, ~15 MB RAM).
- `τ²` is estimated from a sample of rows, not the whole matrix (documented approximation).

See the real spike verdict in `docs/report.md` (Model 1 section).

---

## How to run

```bash
make model          # Model 2 (core)
make spike          # Model 1 spike (optional, requires: pip install h5py s3fs fsspec)
```

## Next steps (not included)

- Strengthen `regpower` with real cross-donor/cross-guide reproducibility from `DE_stats.h5ad`.
- A condition-specific term `γ_{p,c,g}` and a spike-and-slab prior (`z ~ Bernoulli(π)`) for the network.
- Full Bayes (NumPyro/PyMC) if EB stops being sufficient.
