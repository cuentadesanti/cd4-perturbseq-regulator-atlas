# 6. Discussion

## What the operator view buys

The three results share a common source: treating the screen as one object rather than a list.
Each is something the ranking or program-extraction views cannot produce on their own. The
predictive result (Section 3) is the clearest case — a ranked list says nothing about a regulator
outside the panel, whereas the operator's shared structure predicts the unobserved
late-stimulation response of out-of-panel regulators, beating an honest persistence baseline at
every rank. The condition-gated programs (Section 4) exist in program-extraction analyses too, but
here they emerge as factors of a single object whose condition-modulation profile is read off
directly and tested with bootstrap confidence intervals, rather than as separate per-condition
decompositions that must then be aligned. And the universal-versus-specific contrast (Section 5) is
not a property of the CD4⁺ screen at all until it is set against a second cell type; the operator
makes each regulator's effect profile the unit of comparison, so the contrast becomes a
data-built distinction rather than an annotation. The through-line is that shared structure, once
made explicit, is usable — for prediction, for condition resolution, and for cross-cell-type
comparison — in ways a hit list is not.

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

The contribution of this work is not a new estimator or a new master regulator, but a way of
looking: the perturbation screen as a single operator whose geometry is predictive, condition-
resolved, and cell-type-discriminating — recovered with explicit statistical rigor from a 1.8 TB
atlas, on a laptop.
