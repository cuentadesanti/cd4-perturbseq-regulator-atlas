# Abstract

Genome-scale Perturb-seq screens resolve thousands of regulators at once, but the resulting effects
mix genuine signal with noise, connectivity hubs, and context dependence, and it is unclear whether
they assemble into a reproducible regulatory architecture — or whether any such architecture
generalizes beyond the perturbations observed. We address both questions in a genome-scale CRISPRi
Perturb-seq atlas of primary human CD4⁺ T cells by building a quality-audited regulatory operator —
regulators filtered for reproducibility and context, effects taken as a single object rather than a
hit list — and denoising it spectrally. After separating a global co-variation mode and a
Marchenko–Pastur noise bulk, an empirically calibrated **92 signal directions** (not the 336 a
closed-form edge would admit) carry a strongly non-random community structure (consensus modularity
z = 259). Blind, annotation-free community detection recovers a known molecular complex:
Mitochondrial Respiratory Chain Complex I, at BH-FDR 1.4 × 10⁻⁷ — a stable, unsupervised community
matching a curated complex, which validates the denoising, representation, and clustering at once.
Convergent evidence — CORUM identity, CP-factor concentration, and cross-cell-type cohesion —
supports a second, SAGA-centered coactivator module, and the module organization nominates
autoimmune risk regulators against a length-matched null. The universal-versus-T-specific
distinction the operator draws is sustained across three diverse cell types. Yet the same structure
does not translate into prediction: real and permuted regulator features perform identically for
leave-regulator-out inference, so unseen regulators are not predictable from available annotations.
The perturbational operator contains recoverable biological structure — but recoverable structure is
not equivalent to inductive predictability.
