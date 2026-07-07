# 80/20 EDA — CD4+ T cell Perturb-seq

## Scope

Analysis using only the **supplementary tables** (~15 MB). **No** `.h5ad`/`.h5mu` was loaded (1.8 TB).

| File | Rows | Use |
|---|---|---|
| `DE_stats.suppl_table.csv` | 33,983 | main table (the signal) |
| `sgrna_library_metadata.suppl_table.csv` | 26,504 | guide library |
| `sample_metadata.suppl_table.csv` | 12 | experimental design |

Reproducible: `python scripts/eda.py` → figures in `docs/figures/`, table in `docs/tables/`.

## Unit of analysis

**1 row = perturbed gene × culture condition** (`target_contrast` × `culture_condition`).
Rest / Stim8hr / Stim48hr.

## Key quality filters

| Filter | Available here | Source |
|---|---|---|
| `ontarget_significant` (effective KD, 10% FDR) | ✅ CSV | `DE_stats.suppl_table.csv` |
| `offtarget_flag` (possible off-target) | ✅ CSV | same |
| **cross-condition** reproducibility (proxy) | ✅ derivable | min/max of `n_downstream` across conditions |
| `single_guide_estimate` (2 concordant guides) | ❌ | `DE_stats.h5ad` `.obs` |
| `guide_correlation_all` (cross-guide) | ❌ | `DE_stats.h5ad` / `by_guide.h5mu` |
| `donor_correlation_hits_mean` (cross-donor) | ❌ | `DE_stats.h5ad` / `by_donors.h5mu` |

## Main findings

1. **DE effects are heavy-tailed.** Median **2 DEGs**, mean 60.5 (misleading), 15.4% with no effect,
   1.5% are hubs (>1000 DEGs). → summarize with **percentiles and rankings**, not the mean.
2. **Knockdown gates most of the signal.** 62% of contrasts have a significant on-target KD;
   these concentrate **85%** of all trans-effects. Filtering by `ontarget_significant` sharply raises
   signal density (not causal proof: a real KD may go undetected due to low baseline expression).
3. **Stimulated cells show broader effects.** Mean DEGs: Rest 53.1 · Stim8hr 68.9 · Stim48hr 59.4.
4. **Top hubs are plausible**, enriched in T cell signaling (CD3E/D/G, LAT, ZAP70, PLCG1,
   LCP2, VAV1) — suggesting the screen captures interpretable biological signal.
5. **The library covers ~2 guides/gene** (12,440 of 12,654 genes have exactly 2) → internal replication.
6. **Robust regulators ≠ raw hubs.** Requiring cross-condition stability + significant KD + no off-target
   surfaces **chromatin/transcription** regulators consistent across all 3 conditions (TADA2B, TADA1, SGF29,
   SUPT20H — SAGA complex; ELAVL1, NFRKB), distinct from the raw top dominated by Stim8hr-specific TCR signaling.

## Figures

`docs/figures/`
- `01_distribution_n_total_de_genes.png` — long tail
- `02_degs_by_condition.png`
- `03_top_hubs_by_condition.png`
- `04_ontarget_vs_downstream.png` — axis `kd_strength = −ontarget_effect_size`
- `05_guides_per_gene.png`
- `06_reproducibility_vs_effects.png` — cross-condition reproducibility vs. magnitude

## Actionable table

`docs/tables/top_robust_regulators.csv` (top 30). Score using **CSV columns only**:

```
robust_score = log1p(n_downstream)
             · ontarget_significant
             · (0.6 if offtarget_flag else 1.0)
             · (0.5 + 0.5 · cross_condition_reproducibility)
             · (n_signif_conditions / 3)
```

## Practical next steps

1. **Definitive robust ranking**: strengthen `robust_score` with `single_guide_estimate` +
   `donor_correlation_hits_mean` from `GWCD4i.DE_stats.h5ad` (17 GB).
2. **Gene-level downstream matrix**: load the `.layers` of `DE_stats.h5ad` (log_fc/zscore/padj)
   for the top regulators, already filtered.
3. **Per-condition regulatory graph**: build the regulator → downstream network per condition and
   compare Rest vs. Stim to identify context-specific regulators.
