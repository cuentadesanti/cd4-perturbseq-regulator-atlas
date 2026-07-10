# Tier-A prepared prose (NOT integrated — user decides manuscript integration)

Prepared paragraphs from the Block-A / tier-A hardening pass. **These are drafts held here on
purpose; they are not wired into `docs/manuscript/*`.** Numbers are sourced against the shipped
tables (`operator_bcv_rank_3106.csv`, `operator_community_null_3106.csv`).

---

## 1b — Evaluation reframe: state the mean-predictor floor explicitly

> We evaluate out-of-panel prediction against **two** baselines, and report both because they
> answer different questions. The first is the **mean predictor** — assigning every held-out
> gene its across-regulator mean, i.e. R² = 0 by construction. The second is **persistence** —
> carrying the early-stimulation response forward (Stim48hr := Stim8hr), which scores
> R² = −0.31 on the held-out late-stimulation fibers: late stimulation genuinely diverges from
> early, so persistence is *worse than the mean*. The low-rank operator beats **both** — a
> positive out-of-panel R² ≈ 0.07 that clears the mean-predictor floor and exceeds persistence
> by a margin of ≈ +0.38 (`operator_completion_multiseed_3106.csv`,
> `operator_completion_stratified_3106.csv`).
>
> The mean-predictor floor is not a formality: beating it out-of-panel is precisely the bar that
> current deep perturbation-effect predictors do **not** reliably clear — genome-scale
> foundation models fail to outperform simple linear/mean baselines
> ([Ahlmann-Eltze et al. 2025](https://doi.org/10.1038/s41592-025-02772-6)), and generalizable
> single-cell perturbation-response prediction remains an open benchmark
> ([Wei et al. 2025](https://doi.org/10.1038/s41592-025-02980-0)). A small, positive, *honestly
> cross-validated* R² that clears the mean-predictor on regulators never seen at training is
> therefore a real and falsifiable result, not an underachievement relative to an inflated
> baseline.

---

## 1a — Predictive rank vs signal rank (the honest ceiling)

> The correlation structure of the operator carries **≈ 92 signal directions** above the
> empirical Marchenko–Pastur noise edge (2.95; `operator_community_null_3106.csv`), but only a
> handful of these **generalize across the condition split**. Cross-validating the predictive
> rank — bi-cross-validation of the condition-extrapolation task over regulator folds and
> held-out conditions — puts the optimal predictive rank at **≈ 7**
> (`operator_bcv_rank_3106.csv`). The gap is the point: of ≈ 92 reproducible signal directions,
> only ≈ 7 predict the response of unseen regulators; the rest is
> regulator-idiosyncratic. This honestly explains the modest R² ≈ 0.07 ceiling as **mostly
> irreducible** — the predictable subspace is genuinely low-dimensional — rather than as a
> failure of the model to fit.

*Config note (report honestly):* the predictive rank is **target-condition-specific**. For the
flagship late-stimulation (Stim48hr) extrapolation it is **7** [95% CI 7–7], robust to holdout
fraction (rank 7 at 10/20/40%). For the earlier conditions it is higher — **11 (Rest), 12
(Stim8hr)** — and those extrapolations are also *easier* (R² ≈ 0.13 vs ≈ 0.06 for late-stim):
late stimulation is the most divergent, hence the hardest and lowest-rank to predict. Either way
the predictive rank (7–12) is an order of magnitude below the ≈ 92 signal directions.
