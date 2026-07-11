# Manuscript outline — the empirical regulatory operator of CD4⁺ T cells

## Thesis (the sentence everything hangs from)

> Perturbation effects in CD4⁺ T cells are not an independent hit list but a **shared regulatory
> operator**: (i) its **predictive subspace is low-dimensional (~7 components)** and predicts the
> unobserved response of **out-of-panel** regulators, (ii) it decomposes into **gene programs gated
> by activation state** (some peaking in resting cells, others in early stimulation), and (iii) it
> separates **universal from T-cell-specific regulation** — all recovered with explicit statistical
> rigor from a 1.8 TB atlas, on a laptop.

Guardrails inherited from the repo (do not re-open in the manuscript):
- **Low-rank applies to the predictive subspace (~7), not the operator** (full effective rank ~86;
  `OPERATOR_ANALYSIS.md:92`). Never write "the operator is low-dimensional."
- **3b holds out a condition fiber, never a whole regulator** (`OPERATOR_ANALYSIS.md:37`). "Fuera del
  panel / respuesta no observada," never "nunca medidos."
- **Absolute precision softens at scale** (0.154→0.089); the *margin* generalizes, not absolute R².
- **Naming is thinner than gating**: gating is bootstrap-CI proven; program *labels* rest on
  single clean anchors (SUPT20H, DOCK2). State both at the same evidentiary tier.

---

## Section 1 — Introduction

- **Frame:** genome-scale Perturb-seq in *primary* CD4⁺ T cells (Marson 2025), 4 donors × 3
  activation states. The field ranks perturbations or extracts programs; **none builds the
  regulator×gene operator and interrogates its geometry** — that is the gap this work fills.
- **Position against methods literature** (from `literature_positioning.md` §Methodology): front-end
  = calibrated/harmonized DE inherited (Barry SCEPTRE `10.1186/s13059-024-03254-2`; Peidli scPerturb
  `10.1038/s41592-023-02144-y`); contribution = the operator, its condition factors, its
  out-of-panel predictivity. GSFA (`10.1038/s41592-023-02017-4`) = probabilistic analogue, a
  comparison not a scoop.
- **The rigor posture as a stated methodological choice**, not decoration: fewer, better,
  power-controlled claims — against a frontier where DL predictors don't beat baselines
  (Ahlmann-Eltze `10.1038/s41592-025-02772-6`; Wei `10.1038/s41592-025-02980-0`).
- **Closing paragraph — the substrate bridge (load-bearing for structure).** The thesis has three
  legs (operator / programs / K562), but the results open with the regulator ranking (§2, Figure 1),
  which is *not* one of the three. Reconcile this explicitly so the lead figure earns its place: state
  that we **first establish the regulator ranking is trustworthy — metric-robust and
  power-controlled (§2) — then build the operator on that audited substrate (§3–5)**. The ranking is
  the foundation the operator stands on, not a fourth claim; promoting it into the thesis would
  dilute the operator through-line that is the paper's actual novelty. One sentence in the intro's
  close makes Figure 1 motivated rather than orphaned.

## Section 2 — The regulator ranking is metric-robust and power-controlled (foundation)

*Establishes the substrate before the operator. This is the "audited, trustworthy" layer.*
- `n_downstream` ranking; the power confound (ρ≈−0.68 for naive magnitude, wrong sign for
  `n_downstream` at −0.22; 0.90 concordance with FDR-restricted magnitude).
- Knockdown gate (62% sig carry 85% of trans-effects); donor-robustness as a *column*, not a re-sort.
- cis-off-target audit (Step 1): `offtarget_flag` = cis at 99.4%; exclusion legitimate; SAGA
  survives on cis-clean subunits (SUPT20H #4→#3), corrected subunit list.
- **Figures:** `07_hub_posterior_ranking.png`, `08_kd_gate_changes_ranking.png`,
  `11_ranking_stability.png`, `19_reproducibility_aware_ranking_shift.png`.
- **Tables:** `hub_ranking_bayes.csv`, `hub_ranking_bayes_reproducibility_aware.csv`.

## Section 3 — The operator and its predictive subspace [THESIS (i)]

*The flagship. Descriptive → predictive.*
- Build: tensor regulator×gene×condition in pooled z-score space; confound guard (ρ=−0.006 at 3106);
  escalation to full above-floor axis (4×), power-gate holds (0.049).
- **3b out-of-panel prediction:** beats persistence at **every rank in every stratum** (621
  held-out out-of-panel regulators; aggregate margin **+0.379 @ rank 7**). Strength-stratified
  (source: `operator_completion_stratified_3106.csv`): **weak +0.501, strong +0.277 (median split),
  strong panel-matched +0.205** — the margin *decreases* with regulator strength because persistence
  is a progressively fairer baseline for stronger regulators (−0.21 median-split, −0.15
  panel-matched) and collapses for weak ones (−0.42). The
  anti-skeptic point is therefore direct and sourced: even the strongest, panel-matched stratum with
  the fairest baseline still beats persistence (+0.205, min +0.151 across ranks) — the win is not an
  artifact of weak-regulator collapse. **⚠ The old "+0.332 ≈ pilot +0.330 apples-to-apples" line is
  RETIRED** — it never reproduced from disk and its direction was inverted (strong stratum is below
  pilot, not equal to it). Effective predictive rank ~7 (curve peaks at 7, declines — genuine
  turnover, not right-censored).
- **The honest caveat, in-text:** absolute R² softens (0.154→0.089); relative margin generalizes,
  absolute precision does not fully. Keep this next to the margin, **not** in the discussion —
  margin-generalizes/precision-softens is one honest claim; splitting them is how it becomes an
  overclaim in results and a buried caveat later.
- **⚠ Drafting guardrail — do NOT lift these verbatim from `OPERATOR_ANALYSIS.md`:** line 37 says
  the operator predicts "regulators **never characterized**" and line 96's net headline says the
  operator "**is low-rank (~7)**". Both predate the thesis-tightening fixes and reintroduce the exact
  overclaims the thesis guards against. The correct scope is: **out-of-panel** regulators, a held-out
  **condition fiber** (not a whole/never-measured regulator); and low-rank applies to the
  **predictive subspace (~7)**, not the operator (whose effective rank is ~86). Draft §3 from the
  tightened thesis, not from those lines.
- **Figures (anchor):** `35_operator_completion_curve_3106.png` (the margin-vs-rank, the money
  figure), `32_operator_svd_scree_3106.png` (spectrum / the ~7 vs ~86 distinction made visually).
- **Tables:** `operator_completion_condition_3106.csv`, `operator_svd_power_3106.csv`.

## Section 4 — Programs modulated by activation state [THESIS (ii)]

*The CP decomposition; condition identity is the proven part, naming the thinner part.*
- CP gives factors with a condition-modulation profile read off one factor; **gating is bootstrap-CI
  proven** (factor 1 peak-Rest, factor 6 peak-Stim8hr, CIs exclude flat). Pooling re-confirmed
  *positively* at 3106.
- **Naming, non-circular** (assignment by `a_k` blind to annotation — state it): factor 6 = early
  activation coupled with SAGA, anchored on cis-clean **SUPT20H alone** (TADA2B/TAF6L cis-inflated,
  excluded); factor 1 = Rest-gated proteostasis/RNA-processing (downstream-gene name; no curated
  gene-set, so identity + cross-decomposition, *not* a p-value) with a multi-module regulator side
  (AHR–ARNT + immune Rac→WAVE→actin, suggestive-not-decisive split at unstable higher ranks).
- **Figures (anchor):** `34_operator_cp_condition_factors_3106.png` (the gating, the money figure for
  this section), `33_operator_cp_stability_3106.png` (rank selection).
- **Tables:** `operator_cp_factors_3106.csv`, `operator_svd_programs_3106.csv`.

## Section 5 — Universal vs T-cell-specific regulation [THESIS (iii)]

*Cross-cell-type validation against an external dataset; independently recovers the naming.*
- Concordance vs K562 (Replogle 2022 pre-aggregated bulk); per-regulator z-scored Pearson;
  label-permutation null (q95 0.038); **power control (ρ≈0.20 → lead with the well-powered split**,
  508 universal / 227 T-specific; core donor-robust 201/49, stated as a **conservative floor**).
- **The inversion (the strong form):** SAGA→universal (anchor SUPT20H +0.186); immune
  Rac→WAVE→actin→T-specific (clean anchor **DOCK2 alone**; NCKAP1L/AHR/ARNT donor-robust but not
  well-powered — floor-argument, same evidentiary tier as SUPT20H). Two independent analyses
  (§4 naming + §5 K562) converge on immune-cytoskeleton T-specificity.
- **Figure (anchor):** `36_k562_universal_vs_specific.png` (in PR #15 branch).
- **Table:** `k562_concordance.csv` — includes the `min_kd_rank` column (added in fcf8aad), so the
  power diagnostic `spearman(pearson_z, min_kd_rank) = 0.200` reproduces directly from disk. Lives on
  the PR #15 branch; lands on main when #15 merges. (Note: recomputing from the rounded `kd_*` columns
  gives 0.122 because K562 ranks carry NaNs — the `min_kd_rank` column is the authoritative one.)

## Section 6 — Discussion

- What the operator view buys that ranking/program-extraction alone does not: prediction beyond the
  panel, condition-resolved programs as one object, universal/specific as a data-built contrast.
- Honest limits, gathered: predictive-subspace low-rank ≠ operator low-rank; naming rests on single
  clean anchors; absolute precision softens at scale; the one *conceptual* gap (count-model effect
  estimation with guide-assignment uncertainty — SCEPTRE/GSFA territory) named, with the compute
  trade that motivates deferring it.
- Future directions: donor subspace stability (principal angles, `operator_donor_angles.csv`);
  network deconvolution benchmarked against Chevalley 2025 (`10.1038/s42003-025-07764-y`); Option B
  (lfcSE sub-floor) only if the low-power tail is the target.

## Section 7 — Methods

- Dataset & access (1.8 TB atlas; ~15 MB tables + 375 MB K562 bulk + ~0.1 GB z-score fetch; never
  loads cell-level). Tensor build + confound guard. CP (rms scale-control, split-half rank,
  bootstrap CIs, degeneracy exclusion). Completion (train-only standardization, svds soft-impute,
  persistence baseline). K562 concordance (Ensembl join, gene-level collapse, controls dropped,
  z-scored Pearson, permutation null, power/coverage/donor diagnostics). Reproducibility: `make
  operator`, `scripts/_opkernels.py` + 16 tests.

---

## Figure plan (what anchors what)

| Fig | File | Anchors |
|-----|------|---------|
| 1 | `07_hub_posterior_ranking.png` + `08_kd_gate_changes_ranking.png` | §2 ranking foundation |
| 2 | `35_operator_completion_curve_3106.png` | §3 thesis (i) — predictivity **[main money fig]** |
| 3 | `34_operator_cp_condition_factors_3106.png` | §4 thesis (ii) — condition gating |
| 4 | `36_k562_universal_vs_specific.png` | §5 thesis (iii) — universal/specific |
| S1–Sn | scree, stability, fingerprint PCA, donor | supplement |

## What NOT to put in (carried from the rigor pass)

- No GWAS/disease convergence as a result (overclaim magnet; specificity control already walked it
  back — `disease_and_specificity.md`).
- No "predicts the genome" / "the operator is rank-k" / "novel master regulators."
- No fifth analytical axis — four legs (ranking, operator, programs, K562) is a coherent paper.