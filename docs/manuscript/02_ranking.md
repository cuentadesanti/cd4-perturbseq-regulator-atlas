# 2. The regulator ranking is metric-robust and power-controlled

Before treating the screen as an operator, we establish that its most basic summary — which
regulators have the broadest transcriptional effect — is trustworthy. This matters because the
operator is built on the same effect estimates the ranking uses; if those estimates were dominated
by a statistical-power artifact, the operator would inherit it. This section audits the ranking and
shows it survives the three challenges most likely to undermine it: the choice of effect metric,
the confound between apparent effect and assay power, and cis-off-target contamination. It is the
foundation the rest of the paper stands on, not a fourth scientific claim.

## The ranking metric is not an artifact of the metric choice

A recurring criticism of Perturb-seq hub rankings is that they rank a *count* — how many genes pass
a differential-expression threshold — rather than an effect size, so the ranking could be an
artifact of the counting. We tested this directly by building a continuous effect magnitude (the
summed absolute log-fold-change over the FDR-significant genes for each regulator) and comparing it
to the count-based rank. The two agree closely: Spearman ρ = 0.90 between the proper
FDR-restricted magnitude and `n_downstream` (per-gene ρ = 0.93). A naive
all-gene magnitude — summing over every gene rather than the significant ones — is *worse*, not
better, because it sums noise. The ranking is therefore robust to the metric: doing the effect-size
calculation properly reproduces the count-based order, and the order does not depend on which of
the two defensible metrics is used.

## The ranking is controlled for statistical power

The most serious confound is power. Regulators assayed in more cells can show larger apparent
effects for purely statistical reasons, so a ranking that simply correlated with cell count would
be measuring assay depth, not biology. In raw log-fold-change space this confound is real and
large: the per-regulator effect-vector norm correlates with cell count at Spearman ρ = −0.68
(the sign reflects the particular normalization, but the magnitude is the concern). The ranking
metric we use does not carry it: the FDR-restricted magnitude correlates with cell count at only
ρ = −0.21, and `n_downstream` even less. The power confound is a property of the naive magnitude,
not of the ranking, and it is precisely why every downstream operator analysis is performed in the
pooled z-score space (Section 7) where the confound is removed.

Two further properties make the ranked signal interpretable. Effective knockdown gates it: 62% of
contrasts have a statistically significant on-target knockdown, and those concentrate 85% of all
trans-effects — so the broad-effect regulators are overwhelmingly the ones whose intended
perturbation actually took hold. And the ranking is externally corroborated: the top regulators
overlap an independent T-cell CRISPR screen at 25–30 of the top 100 (p < 1e-26), a concordance far
above chance for two screens with different readouts.

## Donor reproducibility is reported, not used to re-sort

The screen spans four donors, and cross-donor reproducibility is informative about which
high-ranking regulators are robust. We carry donor-robustness as an annotation column on the
ranking (`hub_ranking_bayes.csv`), not as a re-sorting criterion. This is deliberate: re-ranking on
donor concordance would silently hide the regulators that rank high but fail reproducibility, which
is exactly the population a reader needs to see. As a column, it lets a reader see "ranked highly,
but fails cross-donor concordance": SMG1, for example, stays at rank 24 in the primary ranking
flagged `donor_robust = False`, rather than being silently demoted. That this matters is shown by
the alternative — a reproducibility-aware re-ranking would push SMG1 from rank 24 to 288
(Figure&nbsp;S2, `docs/figures/19_reproducibility_aware_ranking_shift.png`) — so re-sorting on
reproducibility would bury exactly the high-ranked-but-unreproducible regulators a reader most
needs to see. We therefore keep the primary ranking (`hub_ranking_bayes.csv`) intact and carry
reproducibility as the flag.

## The ranking is robust to cis-off-target exclusion

CRISPRi can knock down not only the targeted gene but its chromosomal neighbors, so an apparent
regulator could be acting through a cis neighbor rather than the annotated target. The screen's
`offtarget_flag` captures this, and we confirmed it is overwhelmingly a cis signal: it agrees with
an independent neighboring-gene-knockdown call at 99.4%, of which the flagged cases are 92% cis.
Excluding flagged regulators is therefore a legitimate specificity control. It does reshape parts
of the ranking — the SAGA chromatin complex is the clearest case, where two subunits (TADA2B,
TAF6L) prove cis-inflated and collapse under exclusion while the cis-clean subunit SUPT20H rises
(rank 4 → 3). The complex's identity survives on its cis-clean members, but the leading-subunit
list is corrected; we carry this correction forward into the program naming of Section 4, where
SUPT20H alone anchors the SAGA program.

Taken together, these audits establish the ranking as an audited substrate: metric-robust,
power-controlled, externally corroborated, transparent about donor reproducibility, and cleaned of
cis-off-target artifacts. The operator built on these effect estimates in the following sections
therefore inherits a controlled foundation rather than an assumed one.
