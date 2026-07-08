# Tier-0 rigor pass — design

**Date:** 2026-07-08
**Status:** approved (design), Tier 0 in implementation
**Trigger:** senior-researcher critique of the CD4+ Perturb-seq submission (six judgment-call findings, not sloppiness).

## The critique, distilled

1. **Wrong effect metric.** Everything ranks on `n_downstream` — a count of FDR-crossing DE
   calls. That conflates biological magnitude with statistical power (dispersion, baseMean,
   cell number). Rank on a continuous quantity instead.
2. **Modeled the summary table, not the data.** The flagship model never touches the
   multivariate fold-change structure it claims to characterize. The fingerprint analysis (which
   does) should be the spine.
3. **Condition-confounded permutation.** TCR members are all stimulation-specific; the panel
   mixes conditions; the null permutes across the whole panel. TCR's z=11 is inflated by "same
   stimulation state," not only complex cohesion.
4. **Donor robustness 81% absent.** Donor metadata covers only 19% of contrasts, so the
   "donor-aware" reweight is inert for 81% of genes. The adjective oversells.
5. **Scope sprawl.** Six explorable objects, each with matching caveat paragraphs. The disease
   bridge is confirmatory by construction.
6. **No external ground truth.** Everything is internal consistency; no comparison to the
   paper's own hits or a held-out donor.

Meta-tell: caveat-to-result ratio. Fix = fewer, better-supported claims, delete the rest.

## Data landscape (verified 2026-07-08)

Public bucket `s3://genome-scale-tcell-perturb-seq/marson2025_data/`:

| object | size | grain | use |
|---|---|---|---|
| `suppl_tables/*.csv` | ~15 MB | contrast/guide | current CSV core |
| `GWCD4i.DE_stats.h5ad` | 15.6 GB | **contrast × gene** log-FC | cached locally as `data/cache/log_fc.f32.npy` **(33983, 10282)** |
| `GWCD4i.DE_stats.by_donors.h5mu` | 15.7 GB | **donor** DE | real donor axis (#4) |
| `GWCD4i.pseudobulk_merged.h5ad` | 41.5 GB | pseudobulk | — |
| `D*_*.assigned_guide.h5ad` ×12 | 110–161 GB ea | **cell** | only source of individual cells → true E-distance |

Key consequence: the cached `log_fc.f32.npy` is **contrast-level**, one summary vector per
perturbation×condition. True energy-distance (scperturb `edist`) is a **cell-level** quantity and
is NOT computable from this matrix — it needs a cell file. Labeling an L2 norm of this matrix
"E-distance" would be a mislabel. We therefore split the metric fix in two:

- **`#1a` magnitude** = `L2 / Σ|log2FC|` of the contrast vector — a continuous, power-decoupled
  effect size, fully local. Honest name: effect *magnitude*, not E-distance.
- **`#1b` true E-distance** = cell-level, deferred to Tier 2, probe-gated.

## Plan — three tiers by credibility-per-byte

### Tier 0 — zero download (this pass)
- **#1a** `scripts/rank_effect_size.py`: KD-gated magnitude ranking from the cached matrix.
  Deliverable = **Spearman(magnitude, n_downstream)** + **top-30 reshuffle**. If it reshuffles,
  that reshuffle is the result; if ρ≈0.95, the count was defensible and we say so.
- **#3** within-condition permutation null in `analyze_fingerprints.py`; report cross-condition z
  (old) and within-condition z (new) side by side.
- **#6a** rank concordance vs the paper's reported hits
  (`emdann/GWT_perturbseq_analysis_2025`).
- **#5** relabel "donor-aware" → "partial guide-level sensitivity check"; one-caveat-per-claim
  pass; disease bridge → one hypothesis line (demote, not delete).
- **#2** reframe: magnitude + fingerprints = spine; EB-on-count = fast CSV-triage appendix.

### Tier 1 — one bounded download (~16 GB): real donor axis
- **#4** stream `by_donors.h5mu` → per-regulator cross-donor concordance (4 donors → 6 pairwise
  correlations); hard flag "concordant in ≥3/4 donors". Feeds **#6b** holdout-donor recovery.

### Tier 2 — real cloud, probe-gated: true E-distance (#1b)
- Probe first: does the cell h5ad carry `obsm['X_pca']`? If yes → stream the embedding only
  (~100s MB/file) and compute `edist` vs non-targeting within condition; pilot on Stim8hr hubs.
  If no → needs the 144 GB counts, only if Tier 0 left the metric question open.

## Judgment calls recorded
- **E-distance naming:** the local metric is effect *magnitude*, never called E-distance.
- **Harden vs cut (bullet #5):** borderline objects (disease bridge, effect network) are
  *demoted to appendix*, not deleted, unless the user prefers deletion.
- **Tier 2 is optional:** if Tier 0's magnitude ranking already reshuffles the top-30 and tracks
  the paper's hits, true cell-level E-distance is not required to answer critique #1.

## Success criteria
- Every headline claim survives: a magnitude metric (not a count), a within-condition null, and
  one external comparison — with one caveat each, max.
