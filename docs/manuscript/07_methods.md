# 7. Methods

## Dataset and access

We analyze the genome-scale CRISPRi Perturb-seq screen in primary human CD4⁺ T cells (four donors;
three activation states: resting, Stim8hr, Stim48hr). The full atlas is ~1.8 TB; all analyses in
this paper run without loading the cell-level data. Instead we use the published supplementary
differential-expression tables (~15 MB), a remote per-condition z-score slice fetched on demand
(~0.1 GB, cached locally), and, for the cross-cell-type comparison, the Replogle et al. (2022)
pre-aggregated K562 bulk product (~375 MB). The entire pipeline is therefore reproducible on a
laptop; no step requires the cell-level `.h5ad`/`.h5mu` files.

## Effect estimates and the power confound

Per-regulator effect profiles are the calibrated differential-expression estimates from the
screen's harmonized pipeline; we do not re-derive them. Because raw log-fold-change magnitude is
confounded with statistical power — regulators assayed in more cells show systematically larger
apparent effects (Spearman ρ between per-regulator effect-vector norm and cell count ≈ −0.69 in
raw space) — all operator analyses are performed in a pooled z-score space, in which each entry is
standardized at the layer level across the full regulator population rather than per condition.

## Operator tensor construction

The operator is assembled as a regulator × gene × condition tensor (`build_operator_tensor.py`).
Genes are the top-variance panel (2,000 genes by default). Regulators are selected by a
per-condition cell-count floor (retaining regulators above the median), which yields 3,106
regulators for the escalated analysis and 800 for the validation pilot. After the z-score fetch, a
confound guard re-measures the rank correlation between per-regulator slab norm and cell count on
the panel genes and refuses to cache the tensor if |ρ| ≥ 0.15 (fail-closed); the escalated tensor
passes at ρ = −0.006, and a power gate on the leading singular subspace holds at max |ρ| = 0.049.

## CP decomposition

The tensor is decomposed with a masked CANDECOMP/PARAFAC factorization
(`decompose_operator_cp.py`; `tensorly` `parafac` with a missingness mask, `n_iter_max = 400`,
`n_init = 10`, `random_state = 0`). Before fitting, condition slabs are RMS scale-controlled so
that a single high-variance condition cannot dominate the fit. Rank is selected by split-half
stability: for each candidate rank we compute the matched-factor cosine between decompositions of
two random halves (subsample 400), and take the largest rank whose stability exceeds 0.70. Each
factor carries a regulator loading, a gene loading, and a three-entry condition-modulation profile;
condition gating is assessed by bootstrap (100 resamples over conditions), and a factor is called
gated only if the bootstrap confidence interval on its modulation profile excludes a flat profile.
Degenerate factors — those with a maximum cross-factor cofactor cosine near 1 (≈0.97) — are
excluded before interpretation. Factor loadings are varimax-rotated (`gamma = 1.0`, convergence on
rotation change, `tol = 1e-6`).

## Out-of-panel completion

Predictive structure is assessed by low-rank matrix completion of a held-out condition fiber
(`operator_completion.py`, kernels in `scripts/_opkernels.py`). We hold out the late-stimulation
(Stim48hr) entries of regulators outside the original variance-selected panel (a seed-0 20%
sample; `--holdout 0.2 --seed 0`), fit a low-rank model to the remaining entries, and predict the
held-out fiber. Standardization is computed on the training entries only
(`train_test_standardize`, centering each gene on its train-mask values) to prevent leakage.
Completion uses soft-impute with a truncated SVD (`soft_impute`, `n_iter = 100`, `tol = 1e-4`); the
truncated solver (`scipy.sparse.linalg.svds`) is numerically identical to a full-SVD soft-impute to
machine precision (a regression test asserts agreement to ~1e-13) and is roughly an order of
magnitude faster. The baseline is persistence — predicting the late-stimulation response equals the
early-stimulation one. We report the model−persistence margin per rank, and stratify the held-out
regulators by knockdown strength (`operator_completion_stratified.py`), with the stratum thresholds
hardcoded and documented in that script: a median split on effect strength, and a panel-matched
stratum (regulators at least as strong as the median in-panel regulator).

## K562 cross-cell-type concordance

Cell-type specificity is assessed against the Replogle et al. (2022) genome-wide K562 Perturb-seq
screen, using its pre-aggregated, Z-normalized perturbation × gene bulk product
(`K562_gwps_normalized_bulk_01.h5ad`). Regulators are collapsed to gene level with an unweighted
mean over guides/promoters, non-targeting controls are dropped (602 removed), and genes are matched
by Ensembl identifier, yielding 6,407 shared regulators. For each shared regulator we compute a
Pearson concordance between its Z-scored CD4⁺ and K562 effect profiles over co-measured genes (no
imputation of non-measured entries), converted to a per-regulator Z-score. Three controls are
applied: a label-permutation null (permuting the cross-dataset pairing; 95th percentile 0.038); a
coverage check (Spearman between concordance and number of co-measured genes, −0.010); and a power
control (Spearman between concordance and the minimum per-side knockdown rank, 0.200 via the
`min_kd_rank` column). Because concordance depends on power, we lead with the well-powered subset
(both screens knocking the regulator down strongly) and report the donor-reproducible core as a
conservative floor. Regulators are classified universal / T-cell-specific / intermediate on the
sign and magnitude of their concordance Z-score.

## Reproducibility

The full operator pipeline runs from `make operator` (tensor → SVD → CP → completion → donor
angles). Numerical kernels are isolated in `scripts/_opkernels.py` and covered by a unit-test suite
(17 tests, `tests/test_opkernels.py`), including the svds/full-SVD equivalence regression test.
Random seeds are fixed (`random_state = 0`, `--seed 0`) throughout. All result tables cited in the
text are written to `docs/tables/` with a `_3106` suffix for the escalated analysis, and the
pilot-scale tables are retained alongside them so the pilot↔escalation comparison is itself
reproducible from disk.
