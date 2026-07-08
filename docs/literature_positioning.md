# How the Robust-Regulator Analysis Relates to the Current Literature

This positions the hackathon submission, a robust-regulator ranking of a genome-scale CD4⁺ T
cell CRISPRi Perturb-seq screen, built with empirical-Bayes shrinkage and a knockdown-gated
signal filter, against three literatures it sits at the intersection of: the Perturb-seq
method lineage, functional-genomic screens in primary human T cells, and empirical-Bayes
statistics for high-throughput expression data.

## Position in the Perturb-seq lineage

Perturb-seq, pooled CRISPR perturbation read out by single-cell RNA-seq, was introduced to
turn genetic screens from single-phenotype readouts into transcriptome-wide ones, so that a
perturbation's effect is a full expression profile rather than a growth score
([Dixit et al. 2016](https://doi.org/10.1016/j.cell.2016.11.038); the companion platform paper
[Adamson et al. 2016](https://doi.org/10.1016/j.cell.2016.11.048)). The scale problem the
submission confronts, a heavy-tailed matrix dominated by a few hubs and a mass of
near-null perturbations, is exactly the regime that emerged once the method went genome-scale in cell lines ([Replogle et al. 2022](https://doi.org/10.1016/j.cell.2022.05.013)),
which mapped ~10,000 perturbations against single-cell transcriptomes and established that most
perturbations produce small effects while a minority reshape global state. The submission's
finding that the median perturbation yields 2 downstream DEGs while ~1.5% are >1000-DEG hubs is
the same effect-size structure, now observed in *primary* human cells rather than a K562/RPE1
line, the harder and less-characterized setting. Two methodological pillars the submission
relies on also trace to this lineage: direct guide-capture sequencing, which makes reliable
guide→cell assignment at genome scale feasible
([Replogle et al. 2020](https://doi.org/10.1038/s41587-020-0470-y)), and the use of the full
expression profile (not a single reporter) as the phenotype, which
[Norman et al. 2019](https://doi.org/10.1126/science.aax4438) pushed furthest by resolving
genetic-interaction manifolds from rich single-cell phenotypes. Positioning the submission
honestly: it is a downstream *analysis* of a primary-cell genome-scale Perturb-seq atlas, and
its novelty is in the ranking methodology and the memory-aware engineering, not in the assay.

## Knockdown-gating and effect structure

The submission's strongest methodological claim, that contrasts with a significant on-target
knockdown (62%) concentrate ~85% of all trans-effects, so `ontarget_significant` is the
mandatory first filter, is well aligned with how the primary-T-cell screening literature treats
perturbation quality. Genome-wide screens in primary human T cells were established as feasible
and informative by [Shifrut et al. 2018](https://doi.org/10.1016/j.cell.2018.10.024), which
recovered known and novel regulators of proliferation and made the case that primary cells,
despite their noise, capture bona fide immune-function biology. The submission's EDA
recapitulation of this, top hubs enriched for canonical TCR-proximal signaling (CD3 chains,
ZAP70, LCK, LAT, VAV1, PLCG1), is a sanity check the field would recognize as the screen
"working." Where the submission adds discipline is in refusing to rank on raw counts: filtering
on effective knockdown before trusting a trans-effect is the analytic analogue of the
guide-efficacy and reproducibility controls that
[Schmidt et al. 2022](https://doi.org/10.1126/science.abj4008) built into their paired CRISPRa/CRISPRi
screens of stimulation responses in primary human T cells, where concordance between orthogonal
perturbation directions was used to separate real regulators from artifacts. The submission does
not have orthogonal-direction data, but its knockdown gate plus cross-condition/cross-guide
reproducibility audit serves the same purpose, it is a within-modality substitute for the
cross-modality control the field considers best practice.

## The chromatin-machinery result, corroborated

The submission's headline biological finding, that robust ranking surfaces
chromatin/transcription machinery (the SAGA complex: TADA1/TADA2B/SGF29/SUPT20H/TAF6L; Mediator:
MED12/CCNC; plus KDM1A/LSD1, SETD2) as large-and-stable regulators, above the
Stim8hr-specific TCR-signaling hubs, is consistent with an independent and prominent line of
T-cell functional genomics. Genome-wide screens of T-cell *exhaustion* converged on epigenetic
and chromatin-remodeling factors as dominant regulators of T-cell state and persistence
([Belk et al. 2022](https://doi.org/10.1016/j.ccell.2022.06.001)), and the broader review
literature frames chromatin regulators as recurrent hits across immune CRISPR screens
([Shi et al. 2022](https://doi.org/10.1038/s41577-022-00802-4)). That an EB ranking built to
reward *stable* effects promotes exactly this class, machinery active across all three culture
conditions, while demoting condition-specific signaling nodes is therefore the expected result
if the method is doing what it claims. It is corroboration of the ranking's construct validity:
the "robust" axis recovers the general-purpose gene-expression machinery, and the
"context-specific" axis (ZAP70, LCK) recovers the stimulation-gated signaling that only fires
under TCR engagement. The distinction the submission draws between global and context-specific
regulators is real biology, not a statistical artifact of the shrinkage.

The one caveat worth stating plainly: recovering chromatin/transcription machinery as top hits
is partly *expected* on general grounds, perturbing core transcriptional coactivators
broadly dysregulates expression, so these genes will have many downstream DEGs in almost any
Perturb-seq screen. The submission's contribution is not the discovery that SAGA matters, but
the demonstration that a principled robustness filter ranks these above the raw-count hubs; the
biological novelty for *T-cell programs specifically* would require the downstream target-gene
analysis the submission lists as a next step.

## Closing the "expected" gap: what the classes actually converge on

The convergence analysis (`make class-programs`, `docs/FINGERPRINT_ANALYSIS.md` for the
fingerprint side) takes that next step and directly addresses the "SAGA is expected" caveat. If
chromatin machinery merely dysregulated expression broadly, its convergent targets would look like
everyone else's; instead a **balanced 30-regulator panel chosen by class** shows the classes
converge on **distinct** target sets (median off-diagonal Jaccard ≈ 0.05). More specifically, the
SAGA-family regulators converge on a **163-gene module that is 24× enriched for interferon-stimulated
genes** (hypergeometric P≈7.5e-21), and the direction is uniform: knockdown *de-represses* ISGs, i.e.
SAGA-family chromatin machinery normally *restrains* the interferon program in CD4⁺ T cells. That an
interferon/ISG axis is a coherent, directional readout of chromatin-regulator perturbation is
consistent with the exhaustion/epigenetics literature that flags chromatin remodelers as state
regulators ([Belk et al. 2022](https://doi.org/10.1016/j.ccell.2022.06.001)); here it is expressed as
a *specific convergent program* rather than a generic hub signature. The condition comparison adds the
expected axis: TCR-signaling programs are stimulation-gated while the chromatin/interferon program is
largely constitutive — the same global-vs-context-specific split the ranking draws, now visible at the
level of downstream targets. This is candidate-program evidence (ISG-flagged convergent targets, not
causal pathway proof), but it moves the biological claim from "SAGA has many DEGs" to "the SAGA class
specifically converges on de-repression of an interferon module," which is the part that was *not*
guaranteed a priori.

## The empirical-Bayes ranking, in context

Shrinking noisy per-gene effect estimates toward a data-estimated prior is the standard remedy
for exactly the problem the submission faces, ranking thousands of genes each measured with
different precision, and it originates with
[Smyth 2004](https://doi.org/10.2202/1544-6115.1027) (limma), whose empirical-Bayes moderation of
gene-level variances became the default for high-throughput differential expression. The modern
single-cell/bulk-DE tool the field reaches for, DESeq2, applies the same philosophy as
empirical-Bayes shrinkage of log-fold-changes and dispersions
([Love et al. 2014](https://doi.org/10.1186/s13059-014-0550-8)), and the DE_stats table the
submission consumes is itself a DESeq2-style output. So the submission's "Model 1" (normal-normal
EB shrinkage of `log_fc`/`lfcSE`) is not an exotic choice, it is the same shrinkage logic the
upstream pipeline already used, applied one level up, at the level of *which regulator* rather
than *which gene*. The submission's honesty about nomenclature, calling this empirical-Bayes /
pseudo-Bayesian rather than a full hierarchical posterior, and flagging that `p_top_1pct` is
"P(exceeds the top-1% threshold)" not "P(in the top 1%)", is the correct level of claim for a
plug-in EB estimate and avoids the over-statement that referees penalize.

What the submission does *not* do, relative to the specialized single-cell-CRISPR statistics
literature, is model the guide- and cell-level count structure directly. Purpose-built methods
for single-cell CRISPR screens estimate perturbation effects from the raw count model with
explicit handling of assignment uncertainty and multiple testing rather than from a pre-summarized
DE table. The submission's design deliberately trades that for a laptop-scale, CSV-only pipeline —
a reasonable engineering choice for a hackathon, and one it states as a limitation. The listed
next steps (folding real cross-guide/cross-donor reproducibility into the ranking; a
spike-and-slab prior for the edge network; full-Bayes via a PPL) are precisely the moves that
would close the gap to the specialized methods.

## Bottom line

The submission is well-grounded. Its effect-size structure, its knockdown-gating principle, and
its chromatin/transcription-machinery result each independently match established findings in the
Perturb-seq, primary-T-cell-screening, and immune-CRISPR literatures, and the empirical-Bayes
ranking is a textbook-appropriate tool applied at a sensible level of aggregation, with unusually
honest scoping of its own claims. The work's originality is not a new biological discovery or a
new statistical method; it is the *packaging*, a reproducible, memory-aware pipeline that turns a
1.8 TB primary-cell atlas into an uncertainty-aware regulator ranking on a laptop, with an
explorable atlas on top. Framed that way, it holds up against the literature it draws on. The
nearest contemporary comparators to cite for context are the atlas-guided TF-discovery efforts in
T cells ([Chung et al. 2025](https://doi.org/10.1101/2023.01.03.522354)) and the paired
activation/interference screens of primary-T-cell stimulation
([Schmidt et al. 2022](https://doi.org/10.1126/science.abj4008)), both of which pursue the same
goal, principled regulator discovery in T cells, with complementary designs.
