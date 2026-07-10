# 4. The operator decomposes into gene programs gated by activation state

Section 3 established that the operator has shared, low-dimensional predictive structure. This
section asks what that structure *is*: when the operator is decomposed into factors, do those
factors correspond to coherent gene programs, and are those programs organized by the biological
variable the screen was designed around — the cell's activation state? We decompose the operator
with a CANDECOMP/PARAFAC (CP) factorization, which represents it as a sum of rank-one components,
each carrying a regulator loading, a gene loading, and a condition-modulation profile. The
condition profile is what makes CP the right tool here: it reads off, for each program, how its
transcriptional magnitude varies across resting, early-stimulated (Stim8hr), and late-stimulated
(Stim48hr) cells.

## Condition gating is the proven result

Two factors are gated by activation state with bootstrap confidence intervals that exclude a flat
profile (Figure&nbsp;3 (`docs/figures/34_operator_cp_condition_factors_3106.png`); source `operator_cp_factors_3106.csv`). One peaks in
resting cells (factor 1: condition weights 0.91 / 0.39 / 0.15 across Rest / Stim8hr / Stim48hr)
and one peaks in early stimulation (factor 6: 0.08 / 0.98 / 0.17). Two further factors are cleanly
constitutive — flat across conditions. That both a resting-peaked and a stimulation-peaked program
emerge, each with a confidence interval excluding flatness, is the load-bearing claim of this
section: the operator's shared structure is not a single undifferentiated block but resolves into
programs whose activity is *modulated by* activation state, in both directions. This is a proven,
bootstrap-supported statement about condition identity, and it is what clause (ii) of the thesis
claims — nothing more.

We are deliberate about the word "modulated" rather than "activated": the two gated programs move
in opposite directions with stimulation. Reading these as "the activation programs" would miss
that one of them is a resting-cell program that stimulation turns *down*. The decomposition
recovers bidirectional gating, and the prose should not flatten it into a one-directional
activation story.

## Naming the programs — non-circular, and honestly thin where it is thin

Assigning a biological name to a factor is a separate, weaker step than proving it is gated, and
we hold it to two independent filters. First, non-circularity: a regulator is assigned to a factor
by its loading in the CP fit, which is computed with no access to any biological annotation — so
recognizing that a factor's top regulators form a known complex is a genuine test, not a
tautology. Second, cis-cleanliness: a regulator can anchor a name only if it survives the
cis-off-target gate of Section 2, because a regulator whose apparent effect came from knocking
down a neighboring gene is not evidence for anything.

Under both filters, **factor 6 reads as an early-activation program coupled to SAGA chromatin
machinery.** Its top regulators include three SAGA subunits — TADA2B, TAF6L, and SUPT20H — but
here the cis gate does real work: TADA2B and TAF6L are cis-inflated, collapsing from ranks 2 and 9
to 93 and 130 under a hard cis-off-target exclusion, so the SAGA anchor rests on the one cis-clean
subunit present, **SUPT20H**. This makes factor 6's SAGA signal thinner than its raw subunit list
suggests — cis-clean, it is essentially one subunit — and we state it that way rather than citing
the full list. Its gene loadings are enriched for early T-cell-activation effector transcripts
(leukocyte actin- and TCR-signaling effectors — CORO1A, ARPC1B, GPSM3, LIMD2, MYL12A; empirical
"TCR signaling|up" markers, hypergeometric p = 1.9e-5), consistent with the Stim8hr peak. The
curated TCR-proximal kinase set (ZAP70/LCK/LAT) shows no overlap with factor 6, as expected —
these are downstream readout genes, not proximal signaling components.

**Factor 1 reads as a resting-gated program spanning proteostasis/RNA-processing on the gene side
and a multi-module regulator side.** Its downstream genes co-load chaperone and RNA-processing
transcripts (TCP1, HSPD1, HSPE1, HNRNP family); because the repository carries no curated
proteostasis gene set, this is an identity-and-cross-decomposition call, not a hypergeometric
p-value, and we label it as such. On the regulator side, factor 1 collects two recognizable
modules: the AHR–ARNT xenobiotic-sensing heterodimer, and the immune Rac→WAVE→actin cytoskeletal
axis (NCKAP1L/Hem1, WASF2, DOCK2, ARHGAP30). Whether these two modules are one program or two that
the rank-6 truncation fuses is a question we can pose but not settle: they co-load at the stable
rank and separate only at ranks 8–10, which sit below the decomposition's stability floor, so the
split is suggestive, not decisive — we cannot distinguish two genuine programs emerging from
ordinary factor rotation in an unstable regime. We report the co-loading and stop there.

The honest summary of this section is a two-tier one: condition gating is proven with bootstrap
CIs; the program *names* are a weaker, filtered overlay, strongest for factor 6's SUPT20H-anchored
SAGA coupling and explicitly tentative for factor 1's internal structure. The next section tests
one of these names against an entirely independent axis of evidence.
