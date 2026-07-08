# Response to the senior-researcher review

This document answers the six-point critique with tested numbers, not caveats. Each item states
what we changed, what we found, and — where the reviewer was right and where the fix was more
subtle than proposed — the honest conclusion. Reproduce with the scripts named below.

## Summary table

| # | Critique | What we did | Result | Verdict |
|---|---|---|---|---|
| 1 | Ranks a DE-count, not an effect size | Built a continuous magnitude over the **FDR-significant** DE set; audited vs power | Proper magnitude ≈ `n_downstream` (ρ=0.90); naive all-gene norm is *worse* | Metric-robust; ranking stands |
| 2 | Modeled the summary table, not the data | Reframed: magnitude + fingerprints are the spine; EB-count is fast triage | (reorg) | Adopted |
| 3 | Condition-confounded permutation | Added a within-condition null; report both | On the reported `zscore` space TCR z 11.2→11.2 (confound negligible) | Claim survives |
| 4 | "Donor-aware" is 81% inert | Relabeled honestly now; real donor axis is Tier 1 | (label fixed) | Adjective dropped |
| 5 | Scope sprawl / caveat armor | Cut nested caveats; disease bridge → one hypothesis line | (edit) | Trimmed |
| 6 | No external ground truth | Concordance vs an independent screen + the paper's own coefficients | Top-100 overlap 25–30/100, p<1e-26 | Externally validated |

## #1 — The effect metric  (`scripts/rank_effect_size.py`)

The reviewer is right that `n_downstream` is a thresholded count, not an effect size. But the
proposed fix (L2 / Σ|log2FC|) has a trap, and testing it properly changes the conclusion.

We streamed `layers/{zscore, adj_p_value}` from the DE h5ad and, per contrast, computed several
magnitudes, then audited each against `n_downstream` and against `n_cells_target` (the power
proxy). Power-decoupling audit (`docs/tables/effect_size_metric_audit.csv`):

| metric | ρ vs n_downstream | ρ vs n_cells |
|---|---|---|
| `n_downstream` | 1.00 | −0.22 |
| `mag_sig` = Σ\|log2FC\| over FDR-significant genes | **0.90** | **−0.21** |
| `mean_sig` = mean\|log2FC\| over FDR-significant genes (intensity) | −0.76 | +0.17 |
| `l2_z` = ‖zscore‖₂ | 0.74 | −0.29 |
| `l2_all` = ‖log2FC‖₂ over **all** genes | 0.44 | **−0.68** |

Three findings:

1. **The naive continuous surrogate is a trap.** L2 / Σ|log2FC| over *all* genes correlates
   **−0.68** with cell count — fewer cells → noisier per-gene log-FC → larger norm. It measures
   imprecision, and it is *more* power-confounded than the count it was meant to replace. It
   discards the real machinery (TADA1/SGF29/MED12/ITK) for low-cell noise (TOMM20/SEC61B).
2. **The fix done properly barely moves the ranking.** Restricting the magnitude to the genes
   that are actually FDR-significant (`mag_sig`) gives ρ=0.90 vs `n_downstream` (per-gene 0.93,
   top-30 overlap 21/30) with the *same* mild power dependence. Among the significant set, "more
   genes" and "more total displacement" are ~the same thing here — so `n_downstream` was a
   defensible proxy, not a category error.
3. **The one genuinely new axis is intensity** — mean|log2FC| among affected genes — which *is*
   power-decoupled (ρ≈0) but is a **different biology** (focal potency), anti-correlated with
   breadth. It is worth reporting as a second dimension, not as the ranking key (see #6).

**Adopted:** rank on `mag_sig` (a proper continuous effect size) as primary, report
`n_downstream` as the equivalent count and `mean_sig` as an orthogonal intensity axis.

## #3 — The permutation null  (`scripts/analyze_fingerprints.py`)

We added a **within-condition** null (each draw matches the tested complex's condition mix)
alongside the original cross-condition null (`docs/tables/fingerprint_complex_validation.csv`):

| complex | cross-condition z | within-condition z |
|---|---|---|
| SAGA | 9.4 | 11.1 |
| Mediator | 3.2 | 4.0 |
| TCR | 11.2 | **11.2** |

The reviewer's mechanism is real: on the **raw log_fc** matrix, where the Stim/Rest axis
dominates variance, TCR's z deflates **19.7 → 14.9** under the within-condition null. But the
reported result is computed on the **standardized zscore** space, where that condition signal is
muted — and there TCR's cohesion is **unchanged** (11.2 → 11.2; its within-condition null mean is
actually lower). So the published z≈11 was **not** materially inflated by condition. We now report
both nulls so the reader can see the confound is representation-dependent and negligible for the
claim as made.

## #6 — External ground truth  (`scripts/external_concordance.py`)

We correlated our per-gene ranking against two references bundled by the original analysis repo
(`docs/tables/external_concordance.csv`):

- **Schmidt & Steinhart 2022** — an *independent* CRISPRi screen (cytokine phenotypes) in CD4⁺ T cells.
- **Polarization regulator coefficients** — the *paper's own* regulator-importance ranking from a
  different modelling task on the same cells.

Results (top-100 overlap by hypergeometric test):

| our metric | vs Schmidt2022 (independent) | vs paper coefficients |
|---|---|---|
| `n_downstream` | ρ=0.17, 25/100 (p<1e-25) | ρ=0.13, 23/100 (p<1e-14) |
| `mag_sig` | ρ=0.15, **30/100** (p<1e-33) | ρ=0.13, 22/100 |
| `mean_sig` (intensity) | ρ=−0.14, 0/100 (n.s.) | ρ=−0.09, 1/100 (n.s.) |

Two conclusions: (a) the ranking is **externally real** — our top regulators are strongly enriched
among an independent screen's functional hits and the paper's own top regulators; it is not
internally-consistent noise. (b) As an **arbiter for #1**, the breadth metrics
(`n_downstream`/`mag_sig`/`l2_z`) are statistically indistinguishable, while **intensity is
negatively concordant** — confirming that focal potency alone does not identify functional
regulators, and retroactively justifying a breadth-based ranking.

## #2 / #5 — Scope and the caveat budget

- **#2:** the spine is now effect-magnitude + fingerprints (the parts that touch the multivariate
  structure); the empirical-Bayes model on the count is presented as fast CSV triage, not the flagship.
- **#5:** the "donor-aware" audit is relabeled to what the data supports — a guide-level
  sensitivity check with **partial (19%) donor coverage** — pending the real donor axis (Tier 1).
  The interferon/disease bridge is reduced to a single hypothesis line; its nested requalifications
  are removed (if a claim needs three caveats to survive, we cut the claim).

## What we deliberately did NOT do (and why)

- **True cell-level E-distance (Tier 2).** It is the field-standard metric and needs a ~140 GB
  cell file. We gate it behind the evidence: the summary-level metrics are mutually indistinguishable
  *and* externally validated (#1, #6), so cell-level E-distance is very unlikely to overturn the
  ranking. We do not stage 140 GB to defend a point already made — that would be the scope sprawl of #5.
- **The real donor axis (#4, Tier 1)** is a bounded next step: stream `GWCD4i.DE_stats.by_donors.h5mu`
  (15.7 GB) for per-regulator cross-donor concordance and a hard "concordant in ≥3/4 donors" flag.
