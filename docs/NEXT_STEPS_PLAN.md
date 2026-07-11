# Next steps — the remaining research-worthy moves

**Status going in.** The operator axis is closed in its honest form (see
[`docs/OPERATOR_ANALYSIS.md`](OPERATOR_ANALYSIS.md), "Escalation to 3106"): the predictive
advantage generalizes to 4× the regulator axis (matched margin +0.332 ≈ pilot +0.330), is
effectively low-rank (~7), with a clean representation (confound ρ = −0.006) and pooling
re-confirmed positively by CP. **Do not spend more effort on the predictive axis** — it is in its
defensible form and pushing it further is chasing decimals. The value now is in the descriptive
cards that raise the ceiling, in dependency order below.

This plan is theoretical/explicit — it names the exact files, columns, and the reasoning for each
move. It is not timeboxed. Execute in order; steps 1→2→3 are ordered by dependency (an artifact in
step 1 poisons steps 2–3 if skipped), step 4 is independent and can run anytime.

---

## Step 1 — Confirm the cis-off-target gate (insurance before you build)

**Why first.** This is the one "obvious but wrong" trap whose failure is *asymmetric*: if
off-target contrasts leaked into the regulator ranking or the operator, every downstream result
(K562, programs) inherits the artifact. It is cheap to confirm and must be confirmed *before*
building anything on top.

**Correction to the earlier framing.** A prior note assumed columns `neighboring_gene_KD` /
`distal_offtarget_flag`. Those do **not** exist in this repo's local data. The real column is
**`offtarget_flag`** in [`data/suppl_tables/DE_stats.suppl_table.csv`](../data/suppl_tables/DE_stats.suppl_table.csv)
— boolean, **2,837 True / 31,146 False**. That is the flag to check.

**What to do, theoretically.**
1. Establish where "effective, on-target perturbation" is defined. The gating predicate lives in the
   scripts that consume the DE table — grep confirms `offtarget_flag` / `ontarget_significant` are
   referenced in `scripts/_opkernels.py`, `scripts/build_operator_tensor.py`,
   `scripts/model_hubs.py`, `scripts/audit_ranking.py`, `scripts/analyze_fingerprints.py`,
   `scripts/model_edges.py`, `scripts/eda.py`, `scripts/donor_concordance.py`.
2. For each of the two *built* analyses that feed the paper — the EB regulator ranking
   (`scripts/model_hubs.py` → `audit_ranking.py`) and the operator tensor
   (`scripts/build_operator_tensor.py`) — confirm the selection filters `offtarget_flag == True`
   *out*, at the same point it applies the `ontarget_significant` gate. The gate should be: keep rows
   where `ontarget_significant == True` **and** `offtarget_flag == False`.
3. **The outcome is binary and both branches are cheap:**
   - **If the flag is already dropped** → write one sentence into
     [`docs/OPERATOR_ANALYSIS.md`](OPERATOR_ANALYSIS.md) and
     [`docs/report.md`](report.md) documenting that effective-perturbation selection excludes
     `offtarget_flag`. That sentence is a reviewer shield.
   - **If it is NOT dropped** → this is a real finding. Re-run the affected selection with the flag
     excluded and report what changed in the ranking / operator. Do this before K562 or programs.

**Acceptance.** A documented, code-traceable statement that both built analyses used the
`offtarget_flag == False ∧ ontarget_significant == True` set — or a corrected re-run if they didn't.

---

## Step 2 — K562 universal-vs-T-cell-specific wiring (highest-value unspent card)

**Why.** This is the cheapest genuinely *comparative* result available, and comparison is what
moves the work from the middle of "research-worthy" toward the top. Partition regulators by whether
their effects agree between a cancer line (K562) and primary CD4 T cells: concordant = universal
wiring, divergent = **cell-type-specific regulation**. The divergent set is the interesting result,
and it survives review because the comparison is built into the data rather than asserted.

**Correction to the earlier framing — read this before starting.** A prior note said "you already
hold the `K562_comparison` table, no fetch needed." **That is false for this repo's local state.**
`data/suppl_tables/` contains only three tables — `DE_stats`, `sample_metadata`,
`sgrna_library_metadata` — and no K562 comparison. K562 appears only in prose in
[`docs/literature_positioning.md`](literature_positioning.md). **So this move requires acquiring the
cross-cell-type comparison data first.** Two honest options:
- **(2a) Fetch the paper's supplementary cross-cell-type table** (the K562↔primary comparison, with
  its `logfc_pearson_r` concordance and `random_r1/r2/r3` negative-control columns) from the
  study's supplement / S3 bucket, into `data/suppl_tables/`. This is the intended path — but it is a
  network fetch, and the exact object must be located first. Do **not** assume the column names
  until the file is in hand; verify them the way we verified `DE_stats` (`.columns` on the first
  rows).
- **(2b) If that table cannot be obtained,** the comparison is not doable as the cold-take
  described. Note the limitation honestly rather than fabricating a substitute — do not invent a
  concordance from data that isn't the real cross-cell-type measurement.

**What to do, theoretically (assuming 2a succeeds).**
1. Load the fetched comparison table; confirm it carries a per-regulator concordance statistic
   (`logfc_pearson_r` or equivalent) and the negative-control columns (`random_r1/r2/r3`).
2. **Set the null from the negative controls, not from zero.** The random columns define what
   concordance looks like by chance in *this* data. A regulator is "universal" only if its
   concordance exceeds the random-control distribution; "T-cell-specific" if it is at/below that null
   (or T-cell-only). This is the step that makes the split defensible — the comparison against an
   internal null, not an asserted threshold.
3. Partition regulators into universal (concordant above null) vs T-cell-specific (divergent), and
   cross-reference the divergent set against the donor-robust regulators already established
   (`docs/tables/program_label_evidence.csv` and the ranking's `donor_robust` column) — a
   T-cell-specific regulator that is *also* donor-robust is the strongest claim.
4. Deliverables: a table `docs/tables/k562_concordance.csv` (regulator, concordance, null quantile,
   class) and a figure `docs/figures/NN_k562_universal_vs_specific.png` (concordance vs null, colored
   by class).

**Acceptance.** A clean universal/specific split with the random-control null explicit, the
divergent set named, and its overlap with donor-robust regulators reported.

---

## Step 3 — Name the programs from the regulator side (highest leverage for the paper)

**Why.** The condition-gated programs pass the pooling test (they're real) but are `program_label ==
unlabeled` in the CP output — see [`docs/tables/operator_cp_factors_3106.csv`](../docs/tables/operator_cp_factors_3106.csv),
where all 6 factors carry `program_label = unlabeled`. "Gated TCR vs constitutive chromatin" is
therefore *not yet earned*. Naming even one program cleanly is worth more to the paper than another
predictive decimal.

**The clean input is already computed.** `operator_cp_factors_3106.csv` has, per factor: the
condition profile (`cond_Rest/cond_Stim8hr/cond_Stim48hr`), the gating verdict (`gating_shape`,
`gated_ci`), and — critically — **`top_regulators`** and **`top_genes`**. The two clean gated
factors are:
- **factor 1** — `gated(peak=Rest)`, `gated_ci=True`; top_regulators `UBXN1;NCKAP1L;AHR;ZMYM2;
  ARHGAP30;ARNT;WASF2;DOCK2`.
- **factor 6** — `gated(peak=Stim8hr)`, `gated_ci=True`; top_regulators `SLC3A2;TADA2B;SENP5;TAF6L;
  SUPT20H;DOLPP1;APPL2;CD28`.
- (factors 2/5 are degenerate — `max_cofactor_cosine ≈ 0.97`, `degeneracy ≈ 0.92` — exclude them;
  factors 3/4 are clean `constitutive(flat)`.)

**The non-circularity requirement (state this explicitly in the writeup).** Naming a program by its
regulators is circular *only if* the regulator set was used to define the program. Here it was not:
the CP assigns regulators to a factor via the loading vector `a_k`, computed from the z-score tensor
**blind to any biological annotation**. So checking whether a factor's high-`a_k` regulators are
coherent by *known* function is a genuine test, not a tautology — but the writeup must say the
assignment was annotation-blind for the claim to hold.

**What to do, theoretically.**
1. **Structural path (lead with this).** For each clean gated factor, take `top_regulators` and ask
   whether they cohere by known function using the annotation already in the repo:
   [`docs/tables/program_label_evidence.csv`](../docs/tables/program_label_evidence.csv) maps
   regulators to known complexes (SAGA/chromatin: `TADA2B;SUPT20H;TADA1;TAF6L;SUPT7L;USP22`;
   Mediator; TCR: `ZAP70;LCK;LAT;CD3E;PLCG1`). Cross-reference the factor's regulators against these
   sets. Note the honest reading of the current data: factor 6 (Stim8hr-gated) contains **CD28**
   (costimulation) and **SLC3A2** (activation-associated transport) *plus* SAGA subunits — a mixed
   membership, not a pure "TCR" factor; factor 1 (Rest-gated) carries **AHR/ARNT** (aryl-hydrocarbon
   pathway) and **DOCK2/NCKAP1L/WASF2** (Rac/actin, immune motility). Report what the regulators
   actually say — do not force a clean "TCR vs chromatin" label the membership doesn't support.
2. **Downstream corroboration (second, not primary).** Take the factor's `top_genes` and test for
   enrichment against curated T-cell signatures, using the marker sets already in the repo:
   [`docs/tables/fingerprint_program_markers.csv`](../docs/tables/fingerprint_program_markers.csv)
   (columns `program, direction, gene, mean_z, consistency`) and the paper's Th1/Th2/activation
   supplementary signatures. Hypergeometric test with BH-FDR; a program the generic offline GO panel
   left anonymous often labels cleanly against the *right* directed gene sets.
3. **Cross-check with the SVD programs** in
   [`docs/tables/operator_svd_programs_3106.csv`](../docs/tables/operator_svd_programs_3106.csv)
   (columns `rotation, pc, gene, loading, tail`) — if a CP gated factor's gene program matches an
   SVD PC's high-loading tail, that is independent corroboration of the same axis from two
   decompositions.
4. Deliverable: fill the `program_label` column in a new
   `docs/tables/operator_cp_factors_3106_labeled.csv` with the earned label (or an honest
   "mixed/unresolved" where the regulators don't cohere), plus one paragraph in
   `docs/OPERATOR_ANALYSIS.md` stating the annotation-blind assignment and what each clean factor's
   regulators support.

**Acceptance.** At least one clean gated factor carries an earned, non-circular label backed by
regulator-side coherence (primary) and downstream enrichment (corroboration); factors whose
membership is mixed are labeled honestly as such, not forced.

---

## Step 4 — (independent) Fold the loose literature memo into the repo

`docs_literature_methods_positioning.md` (the single-cell-CRISPR *statistics* positioning memo:
Squair 2021 pseudoreplication, Peidli 2024 scPerturb/E-distance, Barry SCEPTRE) currently lives as a
**loose artifact in the artifact store, not a file in the repo working tree** (it will not appear in
`ls docs/` — retrieve it from the saved artifacts). Fold its verified citations into
[`docs/literature_positioning.md`](literature_positioning.md) as a "Methodology positioning"
section. **Before folding, re-resolve the one unverified DOI:** Chevalley 2025
(`10.1038/s42003-025-07764-y`) — its only CrossRef check failed and it was never re-verified. Drop it
or confirm it; do not carry an unverified DOI into the repo.

---

## What NOT to do (carried from the operator work)

- **Do not re-open the GWAS/disease convergence** until everything else is locked and the specificity
  control is formalized — the repo already walked back its own interferon/lupus enrichment under a
  specificity control (`docs/disease_and_specificity.md`). It is the overclaim magnet.
- **Do not push the predictive result past its honest ceiling.** The rank-7 plateau and the absolute
  R² softening (0.154→0.089) are the honest form; that is the finished state.
- **Do not add a fifth axis.** Two built analyses (reproducibility-weighted ranking, activation
  rewiring), one comparative win (K562), one labeling pass (programs) is a coherent paper. Breadth
  past that dilutes.
- **Option B (full 6209 via lfcSE per-cell weighting)** is not needed — Option A reached the entire
  above-floor axis. Take it only if the sub-floor low-power tail is specifically the object of
  interest, and only as its own planned project.

---

## Order of execution

1. **Step 1 (cis gate)** — one grep + one sentence, or a corrected re-run. Blocks everything else.
2. **Step 2 (K562)** — highest-value, but confirm the comparison data can be obtained first
   (it is **not** local). If it can't, note the limitation and go to Step 3.
3. **Step 3 (name programs)** — highest leverage for the paper; inputs are all local and computed.
4. **Step 4 (literature fold-in)** — independent, do anytime; re-resolve Chevalley first.