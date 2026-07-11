# 1. Introduction

Genome-scale Perturb-seq — pooled CRISPR perturbation with single-cell transcriptomic readout —
has made it possible to measure the transcriptional consequence of knocking down essentially every
expressed gene in a cell [Dixit 2016 `10.1016/j.cell.2016.11.038`; Replogle 2022
`10.1016/j.cell.2022.05.013`]. In primary human CD4⁺ T cells, a recent
genome-scale CRISPRi screen [Marson Lab 2025 `10.64898/2025.12.23.696273`] extends this to a
setting that matters for immunology and for therapy: four donors, three activation states
(resting, early- and late-stimulated), and genome-scale coverage of candidate regulators. That
screen — the dataset we analyze throughout — is the primary data source for every result below. The
standard way such a screen is read is as a ranked list — which perturbations move the most genes,
which are the "hub" regulators — and that reading is useful but incomplete. It treats each
regulator's effect as an independent entry and discards the possibility that the effects share
structure: that the response to one perturbation is informative about the response to another.

Here we take the opposite starting point. We treat the screen as a single object — the
regulator × gene × condition tensor of perturbation effects, which we call the **regulatory
operator** — and ask what its geometry reveals. This is a different question from ranking
perturbations or extracting gene programs, and to our knowledge it has not been asked of a
primary-cell Perturb-seq screen: existing analyses either rank regulators by effect breadth or
factor the expression matrix into programs, but none builds the regulator × gene operator and
interrogates it as an operator — its predictive structure, its low-dimensional subspace, its
decomposition by condition. That is the gap this work addresses, and the results that follow are
properties of the operator: it predicts within its panel, it is gated by condition, its denoised
structure resolves into reproducible modules — one of them a known molecular complex recovered blind
to annotation — and it separates universal from cell-type-specific regulation, while failing to
predict the effects of regulators it never saw.

Positioning this against the single-cell-CRISPR methodology literature clarifies what is inherited
and what is contributed. The front end — turning raw counts into calibrated, harmonized
per-regulator effect estimates — is a solved problem we build on rather than redo: calibrated
single-cell-CRISPR testing [Barry 2024 `10.1186/s13059-024-03254-2`] and the scPerturb
effect-quantification framework [Peidli 2024 `10.1038/s41592-023-02144-y`] define that layer. Our
contribution sits above it: the operator, its condition factors, and its out-of-panel
predictivity. The closest methodological neighbor is guided sparse factor analysis
[Zhou 2023 `10.1038/s41592-023-02017-4`], which recovers gene programs and their driving
perturbations probabilistically; our CP/SVD operator decomposition is a deterministic analogue of
that object, and we treat it as a comparison rather than a claim of priority. The predictive
result is framed against a frontier where deep perturbation-response predictors do not yet
reliably beat simple baselines [Ahlmann-Eltze 2025 `10.1038/s41592-025-02772-6`; Wei 2025
`10.1038/s41592-025-02980-0`]: we therefore ask not whether a flexible model can fit, but whether
the operator's shared structure beats an honest persistence baseline out-of-panel — a lower, more
falsifiable bar.

A methodological posture runs through all results, and we state it as a choice rather than
leave it implicit: fewer claims, each controlled for the confound that most threatens it. Every
result in this paper is gated by an explicit control — for the power–magnitude confound that
inflates naive effect sizes, for cis-off-target artifacts, for donor reproducibility, and for
decomposition stability — and where a control costs signal, we report the reduced number. In a
field where predictive claims often do not survive a fair baseline, this discipline is part of the
argument, not decoration around it.

A note on structure. The results open (Section 2, Figure 1) not with the operator but with the
regulator ranking: the operator is built on top of the ranking, and the ranking is only trustworthy
if its effect metric is not an artifact of statistical power — Section 2 establishes exactly that,
so everything built on it stands on an audited foundation rather than an assumed one. On that
substrate we (i) construct a quality-audited regulatory operator, filtering effects for
reproducibility and context rather than trusting raw hubs or nominal significance; (ii) denoise it
spectrally and show its community structure is strongly non-random; (iii) recover Mitochondrial
Complex I by annotation-blind clustering, the cleanest evidence the operator contains real
biological organization; (iv) identify a SAGA-centered coactivator module and condition-dependent
programs by convergent evidence; (v) show the universal-versus-specific distinction is sustained
across three cell types; and (vi) demonstrate, with a scrambled-feature control, that unseen
regulators are not inductively predictable from the annotations available. The organizing conclusion
is that a perturbational operator can contain recoverable biological structure while that structure
remains predictively inaccessible from external covariates — a distinction we make explicit rather
than eliding.
