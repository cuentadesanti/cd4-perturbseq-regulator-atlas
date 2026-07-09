# The Empirical Regulatory Operator (z-score, Option-1 pilot)

One matrix, four questions. The **operator** is the effect of every regulator knockdown on every measured gene — the object the fingerprint PCA computed and half-discarded. We build it as a 3-way tensor `T[regulator, gene, condition]` in **precision-decoupled z-score space** over an expanded **~800-regulator** panel, then ask: what are its gene programs, is condition-gating real, is it predictive, and are its programs donor-reproducible.

Reproduce: `make operator` (Step 0 fetch is one-time then cached; Steps 1–4 offline). Kernels are unit-tested (`tests/test_opkernels.py`, 16 tests).

## Representation, and why z-score (Step 0)

Raw log-FC row-norm is **negatively** correlated with statistical power: `spearman(‖row‖, n_cells_target) = −0.683` on significant rows. Norm tracks noise, not signal — so an SVD/CP on raw log-FC would fit heteroscedastic noise first, exactly where the leading factors (the programs) live. The remote `layers/zscore` is a per-(perturbation, gene) z-score (effect/SE), precision-decoupled and pooled across conditions.

**Tensor:** 800 regulators × 2000 genes × 3 conditions (Rest, Stim8hr, Stim48hr). The 800 are the regulators with the broadest downstream footprints (top by `n_downstream`) **above a median cell-count floor** — 733 of them outside the original 200-regulator fingerprint panel.

- **Selection matters.** Ranking the expansion by `ontarget_effect_size` (an earlier draft) reintroduces the power confound (that ranking preferentially selects low-cell / high-z perturbations → z-space confound ρ = −0.27). Breadth (`n_downstream`) above a power floor decorrelates it: the built tensor has **confound ρ = +0.080** (guard requires |ρ| < 0.15). The floor is essential — breadth-ranking *without* it still leaves ρ ≈ −0.36.
- **One fail-closed guard.** Only the confound guard (|ρ| < 0.15) is asserted (the driver refuses to cache otherwise). Two earlier magnitude-based guards (row-norm CV, per-anchor cross-condition spread) were **retired after verification**: on a per-cell-z, breadth-homogeneous panel, *any* magnitude statistic is nearly invariant to the pooled-vs-within-condition distinction. Row-norm CV (0.229 here) is driven by selection homogeneity, not within-condition z — proven in raw space, where no z exists (this selection's raw row-norm CV 0.28 < a random panel's 0.34). Pooling instead rests on a documented evidence chain: (1) it is a **layer-level** property established on the same `layers/zscore` slice (row-norm CV 0.36 on the heterogeneous 200-panel, before any selection); (2) the confound guard passes; (3) within-condition z would show up downstream as **all-constitutive** CP condition factors — which the Step-2 gating test would expose. That downstream test, not a Step-0 magnitude proxy, is the authoritative pooling check.

## Gene programs (Step 1)

Right singular vectors of the pooled fingerprint matrix, oriented by an ISG anchor, both raw and varimax-rotated, enriched offline (hypergeometric + BH-FDR).

- **The power fix worked.** Every one of the top-10 program's left factors is power-clean: `max |spearman(u_k, n_cells)| = 0.121`, none flagged confounded. In raw-log-FC space the leading programs would risk being power axes; in z-space they are not. This is the clean positive that validates the whole representation choice.
- **Top programs are the expected immune ones.** PC1 = interferon (IFN_ISG, FDR = 1.6e-4, clean). PC2 = TCR-proximal (FDR = 0.063, near-significant).
- **Honest limit:** only **1 of the top-5** PCs clears the strict "FDR < 0.05 label AND power-clean" bar. This is a **coverage limit of the deliberately-minimal 4-set offline enrichment reference** (~11–17 gene sets, 2% tail), not a failure of the SVD programs — the programs recovered are exactly interferon and TCR-proximal as predicted. A richer offline `.gmt` (MSigDB Hallmark/Reactome dropped into `data/genesets/`) would give the enrichment fair power; `load_genesets()` already merges any present.

## Condition gating (Step 2) — flagship of the descriptive half, and the pooling proof

Masked CP (`regulator ⊗ gene ⊗ condition`) with RMS per-condition scale control, split-half stability rank selection, per-factor power/degeneracy/collinearity gates, and a **bootstrap CI on each condition factor**. Gating is `gated_ci=True` only if the CI excludes flat.

- **Pilot config:** rank sweep to 6, stability-subsample 250, 50 bootstrap resamples (the full defaults — sweep to 8, 100 resamples — are the Makefile default / an Option-3 cluster run). Confound meter on the tensor = 0.080.
- **The pooling proof HOLDS.** Three **clean** factors (not power-confounded, not degenerate, `max_cofactor_cosine` < 0.7) are `gated_ci=True`, peaked on Stim48hr (factor 2), Stim8hr (factor 3), and Rest (factor 4). Because within-condition-normalized data cannot produce gating (every slab would be independently scaled → all-constitutive), **≥1 clean gated factor is positive proof the representation is genuinely pooled** — this is the load-bearing check that replaced Step 0's retired magnitude guards, and it passes decisively.
- **Reconciliation with the prior permutation result.** An earlier within-condition permutation found TCR's *complex-cohesion z* was not condition-inflated in z-space (11.2 → 11.2). That is a **different quantity** from CP *program-magnitude* gating across conditions: cohesion of a complex's fingerprints vs the magnitude of a program's transcriptional response by condition. Both can hold at once — they do not contradict.
- **Honest limits (three).** (1) **All factors are `unlabeled`** — the 4-set enrichment panel cannot name them (same coverage limit as Step 1), so we have *proven* gating but *unnamed* programs; the "gated TCR vs constitutive chromatin" biological headline is not cleanly supported by naming here. (2) **No clean constitutive factor** — the constitutive factors (1/5/6) form a collinear cluster (`gene_mode_cosine` ≈ 0.8), so they fail the cleanliness filter. (3) Stability stays > 0.7 across all ranks and *peaks* at rank 3 (0.937), so the "largest rank clearing threshold" rule lands at the ceiling (rank 6); the operator's effective rank is more honestly ~3 by the stability peak.

## Prediction (Step 3) — flagship of the predictive half

3b (out-of-panel condition extrapolation) is the headline; 3a (entry-wise) is a sanity check. Train-only standardization throughout (no leakage).

- **3b beats persistence, cleanly, out-of-panel.** Holding out entire (regulator, Stim48hr) fibers for **160 regulators, 100% outside the original panel**, and predicting them from Rest + Stim8hr via the low-rank fit on other regulators: the low-rank model **beats persistence at every rank 1–12** — the result does not hinge on a lucky rank pick. Model R² rises from 0.028 (rank 1) to **0.154 (rank 11)**, versus a **persistence** baseline (Stim48hr := Stim8hr) of **R² = −0.176**. Persistence is negative — late stim genuinely diverges from early stim — so low-rank structure is really predicting unseen late-stim responses for regulators never characterized. **This converts the operator from descriptive to predictive.** No genome-scale-imputation overclaim: pure low-rank cannot predict a regulator with zero observed entries; 3b holds out a condition *fiber*, never a whole regulator.
- **3a (sanity):** low-rank beats a rank-1 baseline at every rank ≥ 2; R² is still rising at rank 12 (elbow right-censored — a wider sweep would locate it). Beating the per-gene-mean baseline is near-trivial and is *not* claimed.
- **One disclosure:** the 2000-gene axis is selected once (top variance) on the *full* tensor, so 3b's held-out fibers live on a gene set chosen with those rows present. This is a diffuse feature-selection-on-full-data effect, **not** target leakage — standardization and the low-rank fit never read held-out *values* (verified: corrupting held-out entries leaves predictions bit-identical). It does not invalidate the 3b result, but a fully airtight version would re-select genes on train rows only.

## Donor reproducibility of programs (Step 4)

Principal angles between top-k gene-program subspaces across **disjoint** donor pairs only — `(1,2)vs(3,4)`, `(1,3)vs(2,4)`, `(1,4)vs(2,3)` — never the 15 inflated overlapping comparisons (the 6 modalities share donors). **Data not local:** the per-donor-pair matrices are not in the cache (only a `donor_obs.csv` summary), so the driver correctly prints `[NEEDS-DATA]`, writes an empty-headed table, and exits 0. Producing them needs a fetch of the per-donor layers analogous to `scripts/fetch_fingerprint_matrix.py`; the honest cross-donor result is then a one-command run.

## Escalation trigger

**Step 3b beat persistence cleanly on out-of-panel regulators.** Per the pilot's pre-registered trigger, that is the evidence to escalate to **Option 3**: the full 6209-regulator axis, completion as the flagship, and the CP stability sweep + bootstrap on cluster compute (the pilot ran a reduced config to fit local runtime). Option 1 was the pilot that earned Option 3 — decide on this evidence, not by relitigating scope.

**Two curve-driven requirements for the 6209 run** (the pilot R² had **not** plateaued at rank 11 — it was still rising at the sweep ceiling): (1) sweep rank **well past 12** to actually locate saturation, otherwise any "effective rank" number is right-censored and unsupported; (2) frame the headline as **predictivity, not compressibility** — the honest claim is "low-rank structure predicts unseen late-stim out-of-panel," not "the operator is rank-k." The curve does not yet support a low-dimensionality claim, and configuring the cluster sweep without (1) risks quietly overclaiming it.

### 6209 runbook

**(a) "Full axis" is a size to be found empirically, not 3106.** The floor (n_cells ≥ 567) is necessary but **not sufficient**: above it, the breadth-ranked confound still grows monotonically with axis size. On the panel's top-2000 genes (the guard's own pooled metric), raw `ρ(‖slab‖, n_cells)` is **−0.05 at N=800, −0.16 at 1600, −0.22 at 2500, and −0.24 across all 3106 above-floor** (the literal full 6209 is −0.53 on panel genes, −0.66 on all 10282). Only **N=800 is verified guard-clean** — the sole size with a *measured z-space* confound (**+0.08**); for every larger size only the raw number exists, and the raw→z map is unreliable (at 800, raw −0.05 flipped to z **+0.08**), so the raw trend marks larger axes as *riskier*, not safe. The floor keeps exactly **3106 regulators (50.0% of 6209)**, but **3106 is an upper bound on regulator count, not a verified-clean axis size.** The escalation must sweep axis size upward, re-fetch z-scores, and re-assert the confound guard (**|ρ|<0.15**) on each candidate, taking the largest that passes — or add per-cell precision weighting with real `lfcSE` (not local; needs the per-cell layers) to include the low-power tail directly. Decide sweep-to-max-clean-size vs lfcSE-weighting *before* the fetch: it sets both the fetch cost and the headline ("N-regulator operator, N found empirically" vs "6209 with precision weighting").

**(b) Starting config (drivers are escalation-ready — CLI knobs only, no code change):**
- `python scripts/build_operator_tensor.py --n-total 1200 --top-genes 2000` (rebuild the tensor at the larger panel; re-fetches z-score once to `data/cache/panel_zscore_<hash>.npy`, writes `data/cache/operator_tensor.npz` after the confound guard)
- `python scripts/operator_completion.py --max-rank 30 --holdout 0.2` (push past 12 to find the plateau; raise further if R² still climbs at 30 → writes `docs/tables/operator_completion_condition.csv`)
- `python scripts/decompose_operator_cp.py --rank auto --scale-control rms --max-rank 12 --stab-subsample 800 --boot-n 100` (dial the reduced-pilot knobs — pilot ran subsample 250 / 50 boots — up to full config)
- `python scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate` (re-run Step 1 on the new tensor; cheap)

**(c) The plateau check gates any effective-rank claim.** An effective-rank number is only reportable if the 3b R²-vs-rank curve actually saturates on the escalated (~1200–1600) axis. No plateau ⇒ report "predictive, not low-dimensional" and do **not** cite an effective rank. This is the one place the cluster run can silently overclaim, so it is the gate on the headline.

## What the nuisance controls bought us

| Control | Naive version would have… |
|---|---|
| **z-score representation** | fit the −0.68 power axis first; leading "programs" would be noise magnitude. (Power gate: max \|ρ\|=0.121, clean.) |
| **breadth + power-floor selection** | (effect-size ranking) re-imported the confound at z-ρ = −0.27; the guard caught it. |
| **confound guard, retired magnitude guards** | either a false block (row-norm CV on a homogeneous panel) or a magnitude proxy that can't see the real failure mode. |
| **RMS condition scale control** | (no control) made *every* factor trivially "gated" by Stim's larger magnitude. |
| **bootstrap-CI gating** | called gating off 3 numbers with no error bars — one resample from flipping. |
| **train-only standardization** | leaked test statistics into the fit → inflated, invalid R². |
| **persistence / rank-1 baselines** | "beating zero" — meaningless; the honest bars are persistence and rank-1. |
| **disjoint donor pairs** | 15 inflated overlapping comparisons instead of 3 honest ones. |
| **gauge fixing + varimax (true, not quartimax)** | sign/scale/rotation ambiguity read as biology. |

## Status

Pilot complete. Descriptive flagship (Step 2): pooling proven via clean gated factors; program *naming* limited by the offline enrichment panel. Predictive flagship (Step 3b): a clean out-of-panel win → escalate. Donor step: awaiting the per-donor fetch.


## Escalation to 3106 (Option A, executed)

Reproduce: rebuild `operator-tensor` at `--n-total 3106`, then the pipeline. Escalation outputs are suffixed `_3106` (the `.csv`/`.png` without a suffix remain the 800-pilot record — the escalation is additive, never destructive, so the pilot↔3106 comparison that *is* the result stays verifiable).

**The axis scales 4× clean, and the predictive advantage generalizes — measured, not projected.** Sweeping `--n-total` with the confound guard as arbiter, the guard-clean ceiling turned out to be the **entire above-floor set — 3106 regulators** (confound ρ = **−0.006**; the z-space confound *improved* as the panel grew — +0.08 at 800 → +0.030 at 2000 → +0.012 at 2500 → −0.006 at 3106 — refuting the raw-space expectation that 3106 would fail). The power gate stayed clean at 4× (max |ρ| = **0.049**, no PC confounded): the larger, weaker-tail panel did not reintroduce power.

**Prediction (3b), the flagship, holds against the pre-registered bar — apples-to-apples.** On 621 held-out out-of-panel regulators, the low-rank model beats persistence at all 50 ranks, **though the model's own R² peaks at rank 7 and declines thereafter** (overfitting) — a robust *advantage* over baseline, on genuinely low-rank structure. Restricted to strong regulators (matched to the pilot's type), the **margin (model − persistence) is +0.332 ≈ the pilot's +0.330** — the relative predictive advantage generalizes to 4× the axis and is not an artifact of test-set composition. (The weak stratum shows a *larger* margin, +0.470, but only because persistence collapses there — −0.427 — not because absolute prediction is better.)

**But the absolute precision softens at scale, and we report it.** The model's absolute R² fell — pilot 0.154 → 0.089 **even in the strength-matched strong stratum**. Restricting to strong regulators recovers only *part* of the drop (0.074 → 0.089), so the residual is a genuine **scale effect**: a larger, more heterogeneous training population makes the shared low-rank structure a slightly worse per-regulator predictor, even at equal strength. **The relative margin generalizes; the absolute precision does not fully.** Reporting both halves is what keeps this research-worthy rather than an inflated headline.

**The out-of-panel predictive structure is effectively low-rank (~7) at this scale.** The 3b curve peaks at rank 7 and declines through rank 50 — a genuine turnover, not the pilot's right-censored "still climbing" (rank-11, sweep cut at 12). This licenses a low-dimensionality claim **for the ~7 components that predict unseen regulators — not the full operator's rank** — and it is a property of this 3106/holdout configuration, **not an intrinsic constant** (a 6209 run may move it).

**CP re-confirms pooling at 3106.** Even at reduced config (`--stab-subsample 250 --boot-n 50`), **2 clean factors are gated with bootstrap CIs excluding flat** (peaks Rest and Stim8hr) and **2 clean factors are constitutive** — positive proof the representation stays genuinely pooled at 4× the axis (cleaner than the pilot, which had no clean constitutive factor). The full-config bootstrap is deferred to cluster but was not needed for the re-confirmation; a reduced-config CP has less bootstrap power, so a null there would not be evidence against pooling — here it is not null, it re-confirms.

**Net (the honest headline):** *the predictive advantage over baseline generalizes to 4× the regulator axis and is low-rank (~7)* — with representation clean (ρ=−0.006), power-gate holding (0.049), pooling re-confirmed, and the explicit caveat that absolute predictive precision softens at scale. **Not** "predicts the genome." Per the escalation trigger, this clean out-of-panel result is the evidence that would justify Option B (the `lfcSE` sub-floor build) — *if and only if* the low-power tail is specifically the target.

*Tooling: `soft_impute` was switched to a truncated-SVD (`svds` top-k) kernel — identical to full-SVD truncation to machine precision (pinned by a 1e-13 equivalence regression test, `test_soft_impute_matches_full_svd_reference`), ~10× faster — which made the rank-50 sweep and the 6209 escalation tractable.*


## Cis-off-target gate audit (Step 1 insurance)

Before building further descriptive cards on top of the ranking and the operator, we confirmed — from the code and the data — what the two paper-feeding analytics actually do with off-target contrasts, and measured the sensitivity to excluding them. Reproduce: `python scripts/audit_offtarget_gate.py` → `docs/tables/offtarget_gate_sensitivity.csv`. Additive throughout: the shipped ranking and the merged-3106 tensor are untouched; the leave-out tensor is a separate gitignored `operator_tensor_cisclean.npz`.

**The flag is cis, and "exclude" is therefore legitimate — proven, not assumed.** `DE_stats.offtarget_flag` matches the tensor obs's `neighboring_gene_KD` at **99.4%**; of flagged rows **92% are cis** (guide knocked down a *neighbouring* gene), 1% distal. Excluding it removes exactly the "cis read as trans" artifact. It covers **11.6%** of the significant rows.

**The ranking was never undefended — off-target is handled three ways.** `model_hubs.py` enters it as a **covariate** (`off_i` in `n_downstream ~ C(culture_condition) + sig_i + off_i`), **annotates** it (`any_offtarget` / "possible off-target"), and **penalises** it in `robust_score` (×0.6 if flagged). The row filter is `ontarget_significant` only, so flagged rows are not *excluded* — but they are covariate-adjusted, flagged, and down-weighted. A hard exclusion leaves the global structure intact: Spearman rank correlation **0.9897**, top-30 overlap **25/30**. What moves is a handful of *individuals* whose height came partly from a neighbouring-gene knockdown (TADA2B #2→#93, TAF6L #9→#130, ELOF1/ATAD5/SUPT7L out of the top-30).

**The operator tensor was the one genuine exposure — and it is robust.** `build_operator_tensor.py` gates on `ontarget_significant_bool` only; it does not reference the cis/distal/low-target flags the obs cache carries. An offline regulator leave-out (drop the **404/3106** regulators with any cis-flagged significant row — no new z-scores fetched — and re-run flagship 3b) leaves the predictive result essentially unchanged: the margin over persistence holds (base **+0.379** @rank 7 → cis-clean **+0.390** @rank 6), beats persistence at every rank 1–12 in both, the low-rank elbow stays **~6–7**, and the slab-matrix effective rank is unchanged (**86.5 → 87.1**). **The leave-out changes the out-of-panel test population (621 → 540 held-out regulators), so the small margin delta is *not* interpretable as an improvement** — the defensible claim is robustness (no degradation), not that cleaning helps.

**SAGA survives, but its leading-subunit list is corrected — and that is the audit working, not failing.** The cis gate discriminates *within* the complex: the cis-clean subunits rise (SUPT20H #4→#3, TADA1 #5→#4, USP22 #27→#22) while the cis-driven ones collapse (TADA2B #2→#93, TAF6L #9→#130, SUPT7L #20→out). In the operator panel, 3/6 subunits (SUPT20H, TAF6L, USP22) anchor the program under the clean gate. The SAGA *identity* is robust to cis exclusion; the specific leading subunits are not the same ones — TADA2B's previous prominence (rank 2) was partly neighbouring-gene off-target. We report the correction rather than bury it.
