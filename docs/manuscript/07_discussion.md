# 7. Discussion

## What the operator view buys

The main result of this work is that a quality-controlled perturbational operator contains
recoverable modular organization beyond noise and hub structure. After spectral denoising isolates
roughly ninety genuine signal directions from an atlas of thousands of regulators, the surviving
regulator–regulator structure is strongly non-random and resolves into discrete co-regulated
modules — organization that is invisible to a ranked hit list and that emerges only when the screen
is treated as a single object.

Complex I provides the cleanest validation, precisely because it was recovered without annotations.
A stability-gated community, defined purely by co-movement of perturbation signatures, matches a
curated molecular complex at FDR 1.4 × 10⁻⁷; the clustering never saw the subunit assignments. When
an unsupervised partition reconstructs a known complex, the denoising, the fingerprint
representation, and the clustering are validated simultaneously — there is no interpretive degree of
freedom to absorb the result.

SAGA illustrates a different evidentiary mode, and the contrast is worth stating plainly because it
reflects methodological maturity rather than a weakness to hide. The SAGA-centered module is not a
stably recovered community — its Leiden partition stability sits below our gate — but three
independent sources converge on the same grouping: database identity, factor concentration, and
cohesion in a second cell type. We distinguish these two classes of discovery — a *recovered
community* versus a *convergent regulatory module* — rather than merge them under one label, because
the difference is real and a careful reader would otherwise find it themselves.

The same object is what makes the rest of the analysis possible: it predicts the unobserved
late-stimulation response of out-of-panel regulators from shared low-rank structure (Section 3),
reads condition-gating off directly as bootstrap-tested factors (Section 5), and turns each
regulator's effect profile into the unit of a cross-cell-type comparison sustained across three
lineages (Section 6). But recoverable structure is not predictive reach. Holding out whole
regulators, real and scrambled side-features perform identically, so the map from a regulator's
identity to its trans-effects is not recoverable from available annotations. The operator organizes
what it has seen; it does not, on this data, generalize to what it has not — a boundary we make
explicit rather than eliding.

## Honest limits

We gather the caveats here rather than leave them distributed, because together they define the
scope of the claims. First, low-rank is a property of the *predictive subspace*, not the operator:
the completion optimum at ~7 components describes the transferable structure, while the operator's
own effective rank is roughly an order of magnitude higher. A reader should not take "low-rank
operator" as a description of the object — only of the part of it that generalizes. Second, the
program names rest on single clean anchors. After the cis-off-target and power filters, the SAGA
program is anchored on SUPT20H alone and the T-cell-specific cytoskeletal module on DOCK2 alone;
the surrounding subunits and modules are consistent but not independently powered, and we have been
explicit about which regulators carry the claim and which ride along on a conservative-floor
argument. Third, absolute predictive precision is modest and softens at scale: the out-of-panel
model R² is 0.06–0.08 across strata, below the 0.154 the smaller pilot reached, so the robust
finding is the *advantage* over baseline, not a high absolute accuracy.

Beyond these bounded caveats there is one genuine conceptual gap. We inherit the effect estimates
from the screen's harmonized differential-expression pipeline rather than jointly modeling the raw
counts with guide-assignment uncertainty — the territory of calibrated single-cell-CRISPR testing
and guided factor models. A fully generative treatment that propagated assignment uncertainty into
the operator would be more principled, and we defer it deliberately: it requires the cell-level
data (1.8 TB) and the compute that entails, against a design goal of a laptop-reproducible
pipeline. Naming this gap is part of the honest accounting, not a claim that our shortcut is free.

## Future directions

Three extensions follow naturally. The donor axis can be examined as a subspace-stability question
— whether the operator's leading structure is preserved across donors, measurable as principal
angles between per-donor operators (`operator_donor_angles.csv`) — turning donor reproducibility
from a per-regulator annotation into a statement about the shared structure itself. The operator's
off-diagonal structure invites a network-deconvolution treatment, which should be benchmarked
against the recent single-cell perturbation network-inference benchmark of Chevalley et al. (2025,
`10.1038/s42003-025-07764-y`) rather than asserted, given that field's finding that inference
methods often fail to beat simple baselines. And the low-power regulator tail, excluded here by the
cell-count floor, could be recovered with a precision-weighted (lfcSE) sub-floor build if those
specific regulators — rather than the well-powered structure — became the target of a follow-up.

The contribution is a way of reading genome-scale perturbation data as an auditable operator whose
geometry contains recoverable modular and context-dependent biological organization, while placing a
measured boundary on what that organization can predict — recovered with explicit statistical rigor
from a 1.8 TB atlas, on a laptop.
