# Data model — Genome-scale CD4+ T cell Perturb-seq (Marson 2025)

The dataset is a **star schema** whose axis is the triple
**(sgRNA guide → perturbed gene) × culture condition × donor**.
Expression is aggregated in a cascade: **cell → pseudobulk → DE statistics**.

## ER diagram

```mermaid
erDiagram
    DONOR ||--o{ SAMPLE : "has"
    SAMPLE ||--o{ CELL : "contains (lane/library)"
    GENE ||--o{ SGRNA : "is target of"
    SGRNA ||--o{ CELL : "detected in (guide_id)"
    GENE ||--o{ CELL : "perturbed_gene_id"

    SGRNA ||--o{ PSEUDOBULK : "aggregates"
    DONOR ||--o{ PSEUDOBULK : "aggregates"
    CELL }o--|| PSEUDOBULK : "aggregated by (guide,donor,cond)"

    GENE ||--o{ DE_RESULT : "target_contrast"
    PSEUDOBULK ||--o{ DE_RESULT : "input to DESeq2"
    DE_RESULT ||--o{ DE_VALUE : "per measured gene"
    GENE ||--o{ DE_VALUE : "measured gene (var)"

    SGRNA ||--o{ DE_BY_GUIDE : "DE per guide"
    DE_RESULT ||--o{ DE_BY_GUIDE : "disaggregates"
    DONOR ||--o{ DE_BY_DONOR : "donor pair"
    DE_RESULT ||--o{ DE_BY_DONOR : "disaggregates"

    SGRNA ||--o{ GUIDE_KD : "KD efficiency"
    GENE ||--o{ GUIDE_KD : "target gene"

    DONOR {
        string donor_id PK "CE0008162..."
        int    age
        string sex
        string ethnicity
        float  weight_kg
        float  height_cm
        bool   smoker
        string blood_type
    }

    SAMPLE {
        string cell_sample_id PK "CD4i_R1_D1_Rest"
        string donor_id FK
        string culture_condition "Rest|Stim8hr|Stim48hr"
        string x10xrun_id "R1|R2"
        string library_id "= cellranger output"
        string library_prep_kit
        string sequencing_platform
        date   harvest_date
    }

    GENE {
        string gene_id PK "ENSG..."
        string gene_name
        bool   mt "mitochondrial"
    }

    SGRNA {
        string sgRNA PK
        string chromosome
        int    pos
        string strand
        string seq
        string target_gene_id FK "curated target"
        string target_gene_name
        string designed_target_gene_id "designed target"
        string guide_type "targeting|non-targeting"
        int    distance_to_closest_target_tss
        bool   putative_bidirectional_promoter
    }

    CELL {
        string barcode PK "obs_names of the h5ad"
        string lane_id FK "-> library/sample"
        string guide_id FK "-> SGRNA | 'multi-guide'"
        string perturbed_gene_id FK "-> GENE"
        string guide_type
        int    total_counts
        int    n_genes_by_counts
        float  pct_counts_mt
        bool   low_quality
    }

    PSEUDOBULK {
        string guide_id FK "composite PK"
        string donor_id FK "composite PK"
        string culture_condition "composite PK"
        int    n_cells
        int    total_counts
        bool   keep_for_DE
        bool   keep_effective_guides
    }

    DE_RESULT {
        string target_contrast FK "gene_id, composite PK"
        string culture_condition "composite PK"
        string target_contrast_gene_name
        int    n_cells_target
        int    n_up_genes
        int    n_down_genes
        int    n_downstream "trans-effects"
        float  ontarget_effect_size
        bool   ontarget_significant
        int    n_guides "aggregated guides"
    }

    DE_VALUE {
        string target_contrast FK "composite PK"
        string culture_condition FK "composite PK"
        string gene_id FK "measured gene"
        float  log_fc
        float  p_value
        float  adj_p_value
        float  zscore
        float  lfcSE
        float  baseMean
    }

    DE_BY_GUIDE {
        string guide_id FK "modality guide_1|guide_2"
        string target_contrast FK
        string culture_condition FK
    }

    DE_BY_DONOR {
        string donor_pair FK "CE..._CE..."
        string target_contrast FK
        string culture_condition FK
    }

    GUIDE_KD {
        string sgRNA FK "composite PK"
        string culture_condition "composite PK"
        float  t_statistic
        float  adj_p_value
        bool   signif_knockdown
        string perturbed_gene_id FK
    }
```

## Entities and their physical origin

| Entity | File(s) | Grain (1 row =) |
|---|---|---|
| **DONOR** | `sample_metadata.suppl_table.csv` (denormalized) | one donor (4) |
| **SAMPLE** | `sample_metadata.suppl_table.csv` | donor × condition × run (11) |
| **GENE** | `.var` of any h5ad (reference) | one measured gene (~18k–36k) |
| **SGRNA** | `sgrna_library_metadata.suppl_table.csv` | one guide (31,109) |
| **CELL** | `D*_*.assigned_guide.h5ad` `.obs` | one cell |
| **PSEUDOBULK** | `GWCD4i.pseudobulk_merged.h5ad` `.obs` | guide × donor × condition |
| **DE_RESULT** | `GWCD4i.DE_stats.h5ad` `.obs` / `DE_stats.suppl_table.csv` | perturbed gene × condition (33,983) |
| **DE_VALUE** | `GWCD4i.DE_stats.h5ad` `.layers` | (perturbation×condition) × measured gene |
| **DE_BY_GUIDE** | `GWCD4i.DE_stats.by_guide.h5mu` | guide × condition |
| **DE_BY_DONOR** | `GWCD4i.DE_stats.by_donors.h5mu` | donor-pair × perturbation × condition |
| **GUIDE_KD** | `guide_kd_efficiency.suppl_table.csv` | guide × condition |

## Keys and main joins

- **Gene** is the central reference entity (`gene_id` = Ensembl `ENSG…`, `gene_name` = symbol).
  It appears in two roles: *perturbed gene* (the guide's target) and *measured gene* (a column of the expression matrix / `.var`).
- **SGRNA.target_gene_id → GENE.gene_id**: each guide points to a gene (note: `designed_target_gene_id`
  may differ from `target_gene_id` due to post-hoc curation; there are ~1–2 guides per gene).
- **CELL.guide_id → SGRNA.sgRNA** (special value `multi-guide` if more than one guide was detected).
  **CELL.lane_id → SAMPLE** (one 10x lane = one cellranger output = one library).
- **PSEUDOBULK** = aggregation of CELL by the composite key `(guide_id, donor_id, culture_condition)`.
- **DE_RESULT** = aggregation by `(target_contrast = gene_id, culture_condition)`; it joins the gene's `n_guides` guides.
  `DE_stats.suppl_table.csv` is exactly the `.obs` of this object in tabular form.
- **DE_VALUE** (in `.layers`: `log_fc`, `zscore`, `adj_p_value`, …) is the N:N relation between
  **DE_RESULT** (obs) and **GENE** (var): for each perturbation×condition, a vector over the measured genes.
- **DE_BY_GUIDE** and **DE_BY_DONOR** are the same structure as DE_RESULT but disaggregated
  (by individual guide, or by donor pair) — they feed the reproducibility metrics
  (`guide_correlation_*`, `donor_correlation_*`) that live in `DE_RESULT.obs`.

### Note on donor IDs
The short labels `D1..D4` (cell-level file names) resolve to the canonical `donor_id`
`CE…` via `sample_metadata` (`cell_sample_id` encodes `run_D#_condition`).
The modalities of `DE_stats.by_donors.h5mu` use the `CE…` IDs joined by `_`.
