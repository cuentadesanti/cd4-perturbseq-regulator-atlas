# 6. Generalization: specificity across cell types, and the predictive boundary

The programs named in Section 5 are a within-dataset result: they come from decomposing the CD4⁺
T-cell operator alone. This section tests one of the strongest naming claims — that the immune
Rac→WAVE→actin cytoskeletal module is T-cell biology rather than generic cellular machinery —
against an entirely independent axis: how the same regulators behave in a different cell type. If
a regulator's perturbation signature is conserved between CD4⁺ T cells and an unrelated lineage,
its regulation is likely universal cellular housekeeping; if the signature diverges, the
regulation is cell-type-specific. The comparator is the genome-wide Perturb-seq screen of Replogle
et al. (2022) in K562, a chronic-myelogenous-leukemia line with no T-cell identity.

## Construction and controls

We use Replogle's pre-aggregated, Z-normalized bulk product — a perturbation × gene signature
matrix, the correct unit for a per-regulator concordance, rather than the 8.8 GB cell-level file.
Regulators are collapsed to gene level and non-targeting controls dropped (602 removed), genes are
matched by Ensembl identifier, and 6,407 regulators are shared between the two screens. For each
shared regulator we compute a per-regulator Pearson concordance between its Z-scored CD4⁺ and K562
effect profiles over co-measured genes, and convert it to a Z-score against a per-regulator null.

Three controls guard the split. A label-permutation null — permuting the cross-dataset pairing so
each CD4⁺ regulator is compared against a random K562 regulator — is centered at zero with a 95th
percentile of 0.038, so the observed concordances are not a baseline artifact. A coverage check
finds essentially no dependence of concordance on the number of co-measured genes (Spearman
−0.010; the data are near-complete), so missingness drives nothing. The decisive control is for
statistical power: concordance correlates with knockdown strength (Spearman 0.200 against the
minimum per-side knockdown rank; source `k562_concordance.csv`, `min_kd_rank` column), meaning a
regulator can look "T-specific" simply because it was weakly perturbed in one screen. We therefore
do not report the raw split; we lead with the well-powered subset, where both screens knocked the
regulator down strongly.

## The result, and the inversion

Among well-powered regulators, 508 are universal and 227 are T-cell-specific; restricting further
to the donor-reproducible core gives 201 universal and 49 T-specific. Both are conservative
counts, and deliberately so — for immune regulators that are lowly expressed in K562, a weak K562
knockdown *is itself* a form of cell-type specificity, so the power filter discards some genuine
T-specific signal along with the noise. The T-specific count is a floor, not an estimate.

The biological pattern inverts the naive expectation. One might guess that the activation-coupled
program (factor 6) would be the T-cell-specific one and the constitutive machinery universal; the
data give the opposite, and it is more coherent for being so. **SAGA chromatin machinery is
universal**: SUPT20H — the cis-clean SAGA anchor from Section 5 — is concordant between the two
cell types (Z = +0.186, well-powered and donor-reproducible). Chromatin/transcription machinery
does the same job in any cell, so its conservation is expected. **The immune Rac→WAVE→actin module
is T-cell-specific**: its clean anchor is DOCK2, a lymphocyte Rac guanine-nucleotide-exchange
factor, which is well-powered and donor-reproducible and shows no cross-cell-type concordance.
This is the same module Section 5's decomposition assigned to factor 1 — so two independent
analyses, the within-dataset CP decomposition and the cross-cell-type comparison, converge on the
identity of the immune cytoskeletal module as T-cell-specific.

We anchor each arm on a single clean regulator, symmetrically, because that is what the evidence
supports. On the universal side, SUPT20H alone survives every filter — TADA2B and TAF6L also score
universal but are cis-inflated (Section 5) and not well-powered. On the T-specific side, DOCK2
alone is well-powered; the module's other members — NCKAP1L (Hem1), AHR, and ARNT — are
donor-reproducible and classified T-specific but fall below the well-powered threshold, so their
T-specific label rests on the conservative-floor argument rather than on independent power. The
strong form of the claim is therefore carried by one clean anchor per arm, SUPT20H and DOCK2, with
the surrounding modules consistent but not independently powered — the same evidentiary discipline
applied on both sides.

The value of this section is that it converts a within-dataset naming call into a cross-dataset,
externally validated distinction: the operator does not merely factor into programs, it separates
regulation that is generic to cells from regulation that is specific to the T-cell lineage, and it
does so in agreement with an independent screen the decomposition never saw.

## Specificity is sustained across a third cell type (RPE1)

The universal-versus-specific distinction drawn above rests, so far, on two cell types — primary
CD4⁺ T cells and leukemic K562 — both of them extremes. A single pairwise comparison cannot separate
a genuine lineage-specificity axis from an idiosyncrasy of one contrast. We add a third, genuinely
distinct context: RPE1, a non-cancerous, hTERT-immortalized, near-euploid, p53-competent epithelial
line, using Replogle's essential-scale screen and the identical concordance pipeline. Because RPE1 is
essential-scale rather than genome-wide, we gate on regulator overlap before proceeding — 962
regulators are shared, past the pre-registered threshold, with clean power and coverage controls.

The distinction is sustained, asymmetrically and honestly. SAGA reappears universal: of the three
SAGA subunits present in the essential-scale set, two clear the null (TADA2B and TAF6L,
z = +0.163 and +0.150 against q95 = 0.121) and the third sits at threshold (SUPT20H, z = +0.117) —
the K562 pattern, now in a third lineage. The T-cell-specific side is carried at the class level:
of the 142 regulators classified T-specific in the CD4↔K562 comparison and also present in RPE1,
none become universal (mean concordance z = −0.013). We flag the honest limit rather than paper over
it: DOCK2, the clean T-specific anchor, is not in RPE1's essential-scale screen at all — it is a
hematopoietic-specific exchange factor, not a core-essential gene — so its arm is carried by the
class-level result, not a direct measurement, and its very absence from a core-essential screen is
itself weakly consistent with an immune-restricted role. Across three diverse backgrounds, SAGA
regulation reads as universal and the immune cytoskeletal module as T-cell-specific — a distinction
built from data in three lineages rather than asserted from one.

## Operator structure does not imply predictability for unseen regulators

That the operator has compressible structure does not mean that structure can be predicted, for a
new regulator, from external covariates. We tested this directly. Holding out whole regulators — all
their conditions at once — we asked whether their trans-response could be predicted from
side-features alone: transcription-factor and DNA-binding-domain annotation, protein–protein network
connectivity, and functional-category composition, covering 98–99% of the panel.

Real and permuted regulator features yielded indistinguishable performance, indicating that the
available annotations contain no detectable signal for leave-regulator-out prediction. A linear
model and a non-linear gradient-boosted learner agree: out-of-sample R² is statistically
indistinguishable from zero, and — decisively — from the same model trained on features permuted
against the responses. When real and scrambled inputs perform identically, there is no learnable
map to recover. Per a pre-registered criterion, a model that fails to beat its shuffled control ends
the inquiry rather than motivating a larger one, so we did not pursue further model classes.

This is a conceptual boundary, not a tuning failure, and it is consistent with the field's finding
that current perturbation foundation models do not beat a train-mean predictor on unseen
perturbations (Ahlmann-Eltze et al. 2025). Low-dimensional organization is not equivalent to
out-of-sample predictability: the operator's structure is rich enough to group known regulators into
reproducible modules and to extrapolate across conditions for regulators it has seen, yet the map
from a regulator's identity to its trans-effects is not recoverable, on this data, from the
covariates available. We state the boundary rather than obscure it.
