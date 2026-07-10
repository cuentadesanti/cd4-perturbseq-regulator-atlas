# Abstract (draft)

Genome-scale Perturb-seq screens are typically read as ranked lists of individual regulator
hits. We argue that, in primary human CD4⁺ T cells, the perturbation effects instead form a
shared **regulatory operator** that can be estimated, interrogated geometrically, and used to
predict — and that treating it as one object, rather than a hit list, recovers biology the ranked
view cannot. Working from a genome-scale CRISPRi Perturb-seq atlas
(four donors × three activation states; ~1.8 TB raw) entirely on a laptop, we assemble a
regulator × gene × condition effect tensor in a power-decoupled z-score space (confound
ρ = −0.006) and establish three results. First, the operator's **predictive subspace is
low-dimensional (~7 components)**: a low-rank fit predicts the unobserved late-stimulation
response of out-of-panel regulators, beating a persistence baseline at every rank across a 4×
expansion of the regulator axis (3,106 regulators; aggregate margin Δ R² = +0.379). The advantage
survives strength stratification — even the strongest, panel-matched regulators, for which
persistence is a fair baseline, stay positive (+0.21) — so it is not an artifact of weak-regulator
baseline collapse. This converts the operator from descriptive to predictive, with the honest
caveat that the absolute precision is modest (model R² ≈ 0.06–0.08) and softens at scale. Second, a CANDECOMP/PARAFAC decomposition resolves **gene programs gated by activation
state** (bootstrap-CI–supported factors peaking in resting versus early-stimulated cells),
assigned to regulators blind to annotation. Third, cross–cell-type concordance against an
external K562 screen **separates universal from T-cell-specific regulation**: chromatin/SAGA
machinery is universal (anchor SUPT20H), whereas the immune Rac→WAVE→actin cytoskeletal module
(anchor DOCK2) is T-cell-specific and donor-reproducible — an independent screen recovering the
same axis the decomposition named. Every claim is gated by an explicit control for the
power–magnitude confound, cis-off-target artifacts, donor reproducibility, and decomposition
stability; caveats are reported in the direction that costs signal. The contribution is not a new
estimator or a new biological hit, but the operator view itself — inheriting calibrated
differential-expression front-ends and adding low-rank predictivity, condition-resolved programs,
and cell-type specificity on top.
