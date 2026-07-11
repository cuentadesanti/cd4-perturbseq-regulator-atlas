# 4. Denoised operator structure reveals regulatory modules

The preceding sections establish which effects are trustworthy: regulators identified above a
reproducibility- and context-audited bar, not raw hubs or nominal significance. This section asks
what those trustworthy effects build when taken together. The operator induces a
regulator–regulator similarity, and that similarity is not structureless — it carries discrete,
co-regulated modules. This is where the analysis changes scale, from individual regulators to the
organization of the system, and it is the strongest use of the shared structure in this work: we
recover the modules blind to any annotation, and only then ask what they are
(Figure&nbsp;3 (`docs/figures/38_operator_insignia_3106.png`)).

## 4.1 Spectral denoising isolates low-dimensional regulatory structure

Each regulator is represented by its 2,000-gene z-score fingerprint with the three conditions
concatenated (3,106 × 6,000, aspect ratio q = N/T = 0.52), and the regulator–regulator correlation
of these fingerprints is the object we cluster. Raw, it is dominated by two nuisance structures that
random-matrix theory isolates cleanly. A single global co-variation mode (leading eigenvalue
λ₀ = 167, an order of magnitude above the rest) is the shared perturbation-response axis loading
every regulator; it is deflated. The Marchenko–Pastur noise bulk (fitted edge λ₊ = 1.49 after
deflation) is discarded. What remains — the correlation reconstructed from the signal eigenvalues
between the bulk edge and the global mode — is the denoised operator.

The count of genuine signal directions must be set honestly, because the closed-form MP edge is a
*lower bound*: the 6,000 feature columns are not independent (the same gene appears in all three
conditions), so the effective sample size is below 6,000 and noise eigenvalues leak upward. An
empirical null fixes this — the 95th percentile of the top eigenvalue across 100 column-permuted
correlations gives an edge of λ₊ = 2.95, leaving **92 signal eigenvalues, not the 336 the closed
form would admit**. We use 92 throughout. The point of this subsection is that the clustering below
runs on genuine low-dimensional structure, not on noise the denoising failed to remove.

## 4.2 The denoised operator contains highly non-random community structure

Consensus community detection (Leiden with an RBConfiguration objective, swept over resolution and
50 random seeds, 500 partitions consensus-aggregated) yields eight communities, three of which pass
a partition-stability gate of s_c ≥ 0.8. The structure is overwhelmingly non-random: the consensus
modularity Q = 0.506 sits **z = 259** above a label-permutation null. This is a property of the
denoised operator, established before any community is named — the community structure is strongly
supported relative to the matched null, and the
question of what the modules correspond to is therefore worth asking.

## 4.3 Blind community detection recovers mitochondrial Complex I

The cleanest result in this work is what the clustering produces with no annotation in the loop.
Of the three stability-gated communities, one (n = 87, 45% donor-reproducible) matches
**Mitochondrial Respiratory Chain Complex I** (NADH dehydrogenase) under CORUM enrichment at
BH-FDR 1.4 × 10⁻⁷ — 8 of the 13 CORUM holoenzyme members, with 10 NDUF-family genes present in the
community overall, corroborated across seven Complex-I assembly and intermediate entries. **An annotation-blind partition recovered a known molecular complex.** The
community was defined purely by co-movement of perturbation signatures across the transcriptome and
only afterward matched to a database; nothing about its construction knew the subunits belonged
together. This single result validates the denoising, the fingerprint representation, and the
clustering at once — a stable, unsupervised community corresponding exactly to a curated complex is
not something noise or a flexible interpretation can produce.

That a mitochondrial electron-transport module surfaces coherently under CD4⁺ perturbation is
expected biology, not a housekeeping artifact: ETC respiration is induced and required for T-cell
activation, proliferation, and memory formation (Tarasenko et al., *Nat. Commun.* 2025), so knocking
down its subunits produces one coordinated, activation-relevant transcriptional response — which is
why the subunits co-cluster into a stable, donor-robust community. The other two stability-gated
communities (n = 177 and 132) match no annotated CORUM complex at FDR < 0.05: real, reproducible
regulator modules the atlas has not otherwise surfaced, reported as such rather than named by force.

## 4.4 Convergent evidence identifies a SAGA-centered regulatory module

The second module is a richer inference resting on a different, and honestly weaker, evidentiary
mode. Six SAGA chromatin-coactivator subunits (SUPT20H, SUPT7L, TADA2B, TAF6L, USP22, SGF29) fall in
one community, with only KAT2B elsewhere — but this community's partition stability is 0.56, *below*
the 0.8 gate. We therefore do not claim that stable community detection recovered SAGA. Instead the
SAGA-centered module rests on three independent lines converging on the same grouping. CORUM
enrichment names it: GCN5-linked SAGA at BH-FDR 3 × 10⁻⁴ (8/8), STAGA at 10⁻³, core SAGA and Mediator
below 0.05. The CP factors carrying SAGA loadings (factors 2, 3, 5) map to it at 100% concentration.
And rebuilding the regulator–regulator correlation *in K562* from Replogle's independent screen, the
same subunits re-form a coherent cluster — mean intra-correlation 0.244 against 0.010 ± 0.014 for
random same-size sets (z = 16.1, p < 10⁻⁴), with five of the six shared subunits scoring universal
in the cross-cell-type concordance of Section 6.

The distinction from Complex I is deliberate and we keep it visible throughout: Complex I is a
**recovered community** — isolated directly by stable, annotation-blind clustering — while SAGA is a
**convergent regulatory module**, supported by agreement across CORUM identity, CP-factor
concentration, and cross-cell-type cohesion, but not by a sufficiently stable Leiden partition.
Stating the 0.56 openly is what licenses the convergent-evidence claim; conflating the two under one
"community" label would blur a real methodological difference.

## 4.5 Module structure nominates disease-linked regulators

Known biology validated the structure; the module organization can, in turn, generate nominations.
Because the communities are groups of T-expressed regulators, we ask whether any carries more
autoimmune genetic-risk signal than chance — as a nomination, not a causal claim, against a null
built to avoid the obvious trap. Immune genes are long, SNP-dense, and clustered near the MHC, so a
genome-wide null is trivially "enriched"; we instead resample 10,000 length-decile-matched sets from
the 3,106 panel regulators themselves, reported with and without the chr6 MHC block (Open Targets
genetic_association, five autoimmune diseases). The controls validate the null exactly: the recovered
Complex-I community is not enriched (fold 0.45, p = 0.98) — a leaky null would falsely light it up —
while the SAGA-centered module is enriched (fold 1.35, p = 0.0025, surviving MHC exclusion) broadly
across four diseases (RA/SLE/MS/T1D loads 15.1/8.4/8.3/9.7), including canonical T-cell risk
regulators GATA3, CBLB, UBASH3A, IL12RB2, and RASGRP3. The two stable unnamed communities carry no
signal above the matched null (folds 0.97 and 1.31; p = 0.55 and 0.10). The coactivator module that
organizes T-cell activation touches autoimmune risk genes; the metabolic module does not — a
nomination a follow-up can act on, illustrating what the validated structure enables rather than
promising a clinical result.
