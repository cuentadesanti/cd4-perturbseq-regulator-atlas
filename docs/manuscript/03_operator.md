# 3. A shared regulatory operator with a low-dimensional predictive subspace

Having established that the regulator ranking is trustworthy (Section 2), we now treat the screen
not as a ranked list but as a single object: the regulator × gene × condition tensor of
perturbation effects, which we call the regulatory operator. The question this section answers is
whether that object has *shared* structure — structure that lets the measured perturbations
inform the unmeasured response of others — or whether each regulator's effect profile is
independent, in which case the operator is merely a stack of unrelated rows and the ranked-list
view loses nothing.

## Construction

We assemble the operator as a regulator × gene × condition tensor in a pooled z-score space,
where each entry is a knockdown effect standardized at the layer level rather than per condition.
This choice is load-bearing: raw log-fold-change magnitude is confounded with statistical power
(regulators assayed in more cells show larger apparent effects; Section 2), and a tensor built in
raw space would encode that confound as structure. In the pooled z-score space the confound is
removed — the rank correlation between per-regulator slab norm and cell count is ρ = −0.006 on
the panel genes the decomposition uses, and a power gate on the singular subspace holds at a
maximum |ρ| of 0.049. We build on the set of regulators above a per-condition cell-count floor,
which yields 3,106 regulators — a fourfold expansion over the 800-regulator pilot on which the
approach was first validated. The confound guard is re-asserted after the expanded fetch and
passes, so the fourfold-larger operator is estimated in the same clean representation as the
pilot.

## The operator predicts the unobserved response of out-of-panel regulators

The central claim is predictive. We hold out an entire condition fiber — the late-stimulation
(Stim48hr) response — for regulators that lie outside the original variance-selected panel, and
ask whether a low-rank fit trained on the other regulators, reading only the resting and
early-stimulation entries of the held-out regulators, recovers their held-out late-stimulation
response. The baseline is persistence: predicting that the late-stimulation response equals the
early-stimulation one. Persistence is not a straw man — if late stimulation simply extended early
stimulation, it would be hard to beat; its negative R² (−0.305 on the full held-out set) shows
that late stimulation genuinely diverges from early, so any method that beats it is using
cross-regulator structure to anticipate that divergence.

On 621 held-out out-of-panel regulators, the low-rank model beats persistence **at every rank**,
with an aggregate margin of Δ R² = +0.379 at rank 7 (Figure&nbsp;{{fig:completion_curve}}). Beating
the baseline at every rank, rather than at a single tuned rank, is what makes the result a
property of the operator rather than of a lucky hyperparameter.

The size of that margin, however, must be read carefully, and stratifying by knockdown strength
is what keeps it honest (Table&nbsp;{{tab:stratified}}; source
`operator_completion_stratified_3106.csv`). The margin is *largest* for weak regulators
(+0.501, n=310) and *smallest* for strong ones (+0.277 on a median split, n=311; +0.205 on the
subset matched to the original panel's strength, n=203). This ordering is mechanistic, not
mysterious: persistence is a fair baseline for strongly knocked-down regulators (its R² is only
−0.146 in the panel-matched stratum) and collapses for weakly knocked-down ones (−0.421), so the
aggregate margin is, if anything, concentrated where the baseline is weakest. The conservative
reading is therefore the one we adopt: even for the strongest, panel-matched regulators — where
persistence is a genuinely fair baseline — the model still beats it, by +0.205 at rank 7 and by
at least +0.151 at every rank. The out-of-panel predictive advantage is not an artifact of
weak-regulator baseline collapse; it holds where the baseline is hardest to beat.

Two boundaries on this claim are worth stating explicitly, because both are places a reader could
otherwise over-read the result. First, the held-out unit is a condition *fiber*, not a whole
regulator: each held-out regulator is observed in the resting and early-stimulation conditions,
and only its late-stimulation response is predicted. A purely low-rank model cannot predict a
regulator with no observed entries at all, and we make no such claim — the result is extrapolation
across a regulator's conditions, not imputation of never-measured regulators. Second, the absolute
predictive precision is modest: the model's held-out R² is in the range 0.06–0.08 across all
strata, well below the 0.154 the pilot reached on its smaller 160-regulator holdout. The larger,
more heterogeneous 3,106-regulator population makes the shared low-rank structure a modestly worse
per-regulator predictor in absolute terms, even at matched strength. The beats-persistence
*advantage* is robust and reproducible across strata and ranks; the absolute precision is modest
and does not match the pilot. We report both halves, because reporting only the first would be an
inflated headline.

## The predictive subspace is low-dimensional

The rank at which prediction is best is itself informative. The out-of-panel completion R² peaks
at rank 7 and declines thereafter through rank 50 (Figure&nbsp;{{fig:completion_curve}}) — a genuine
turnover, not a curve still climbing at the edge of the swept range. This licenses a
low-dimensionality statement, but a carefully scoped one: it is the **predictive subspace** that
is low-dimensional (~7 components), not the operator as a whole. The operator's own effective rank,
measured on its singular spectrum, is roughly an order of magnitude higher
(Figure&nbsp;{{fig:svd_scree}}); the ~7-dimensional structure is the part of it that generalizes to
unseen regulators. The distinction matters — a reader who took "low-rank operator" at face value
would expect a spectrum that collapses after seven components, and the spectrum does not. What is
low-rank is the shared, transferable structure, and it is a property of this 3,106-regulator
holdout configuration rather than an intrinsic constant of T-cell regulation.

## The representation stays pooled at scale

A fourfold expansion of the regulator axis could in principle break the pooling that makes the
operator one object rather than many: if the added regulators occupied their own disjoint
subspace, the "shared" structure would be an average over unrelated blocks. A CANDECOMP/PARAFAC
decomposition of the expanded operator rules this out. Two factors are gated by condition with
bootstrap confidence intervals excluding a flat profile (peaking in resting and in
early-stimulated cells respectively; Section 4), and two further factors are cleanly constitutive
— a cleaner decomposition than the pilot, which recovered no clean constitutive factor. That the
decomposition finds coherent, condition-structured, partly-shared factors at 3,106 regulators is
positive evidence that the representation remains genuinely pooled at the expanded scale, not a
mosaic of independent per-regulator effects.

Taken together, these results convert the operator from a descriptive summary into a predictive
object: a shared, low-dimensional structure that anticipates the unobserved late-stimulation
response of regulators outside the panel, robustly beats a fair baseline, and remains one pooled
object at four times the scale at which it was first validated — with an absolute precision we
report honestly as modest.
