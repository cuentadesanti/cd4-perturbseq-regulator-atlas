# Modelo de datos — Genome-scale CD4+ T cell Perturb-seq (Marson 2025)

El dataset es un **esquema en estrella** cuyo eje es la tripleta
**(guía sgRNA → gen perturbado) × condición de cultivo × donante**.
La expresión se agrega en cascada: **célula → pseudobulk → estadísticos de DE**.

## Diagrama ER

```mermaid
erDiagram
    DONOR ||--o{ SAMPLE : "tiene"
    SAMPLE ||--o{ CELL : "contiene (lane/library)"
    GENE ||--o{ SGRNA : "es diana de"
    SGRNA ||--o{ CELL : "detectada en (guide_id)"
    GENE ||--o{ CELL : "perturbed_gene_id"

    SGRNA ||--o{ PSEUDOBULK : "agrega"
    DONOR ||--o{ PSEUDOBULK : "agrega"
    CELL }o--|| PSEUDOBULK : "agregada por (guide,donor,cond)"

    GENE ||--o{ DE_RESULT : "target_contrast"
    PSEUDOBULK ||--o{ DE_RESULT : "input a DESeq2"
    DE_RESULT ||--o{ DE_VALUE : "por gen medido"
    GENE ||--o{ DE_VALUE : "gen medido (var)"

    SGRNA ||--o{ DE_BY_GUIDE : "DE por guía"
    DE_RESULT ||--o{ DE_BY_GUIDE : "desagrega"
    DONOR ||--o{ DE_BY_DONOR : "par de donantes"
    DE_RESULT ||--o{ DE_BY_DONOR : "desagrega"

    SGRNA ||--o{ GUIDE_KD : "eficiencia KD"
    GENE ||--o{ GUIDE_KD : "gen diana"

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
        string library_id "= output cellranger"
        string library_prep_kit
        string sequencing_platform
        date   harvest_date
    }

    GENE {
        string gene_id PK "ENSG..."
        string gene_name
        bool   mt "mitocondrial"
    }

    SGRNA {
        string sgRNA PK
        string chromosome
        int    pos
        string strand
        string seq
        string target_gene_id FK "diana curada"
        string target_gene_name
        string designed_target_gene_id "diana de diseño"
        string guide_type "targeting|non-targeting"
        int    distance_to_closest_target_tss
        bool   putative_bidirectional_promoter
    }

    CELL {
        string barcode PK "obs_names del h5ad"
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
        string guide_id FK "PK compuesta"
        string donor_id FK "PK compuesta"
        string culture_condition "PK compuesta"
        int    n_cells
        int    total_counts
        bool   keep_for_DE
        bool   keep_effective_guides
    }

    DE_RESULT {
        string target_contrast FK "gene_id, PK compuesta"
        string culture_condition "PK compuesta"
        string target_contrast_gene_name
        int    n_cells_target
        int    n_up_genes
        int    n_down_genes
        int    n_downstream "trans-efectos"
        float  ontarget_effect_size
        bool   ontarget_significant
        int    n_guides "guías agregadas"
    }

    DE_VALUE {
        string target_contrast FK "PK compuesta"
        string culture_condition FK "PK compuesta"
        string gene_id FK "gen medido"
        float  log_fc
        float  p_value
        float  adj_p_value
        float  zscore
        float  lfcSE
        float  baseMean
    }

    DE_BY_GUIDE {
        string guide_id FK "modalidad guide_1|guide_2"
        string target_contrast FK
        string culture_condition FK
    }

    DE_BY_DONOR {
        string donor_pair FK "CE..._CE..."
        string target_contrast FK
        string culture_condition FK
    }

    GUIDE_KD {
        string sgRNA FK "PK compuesta"
        string culture_condition "PK compuesta"
        float  t_statistic
        float  adj_p_value
        bool   signif_knockdown
        string perturbed_gene_id FK
    }
```

## Entidades y su origen físico

| Entidad | Archivo(s) | Grano (1 fila =) |
|---|---|---|
| **DONOR** | `sample_metadata.suppl_table.csv` (desnormalizado) | un donante (4) |
| **SAMPLE** | `sample_metadata.suppl_table.csv` | donante × condición × run (11) |
| **GENE** | `.var` de cualquier h5ad (referencia) | un gen medido (~18k–36k) |
| **SGRNA** | `sgrna_library_metadata.suppl_table.csv` | una guía (31.109) |
| **CELL** | `D*_*.assigned_guide.h5ad` `.obs` | una célula |
| **PSEUDOBULK** | `GWCD4i.pseudobulk_merged.h5ad` `.obs` | guía × donante × condición |
| **DE_RESULT** | `GWCD4i.DE_stats.h5ad` `.obs` / `DE_stats.suppl_table.csv` | gen perturbado × condición (33.983) |
| **DE_VALUE** | `GWCD4i.DE_stats.h5ad` `.layers` | (perturbación×condición) × gen medido |
| **DE_BY_GUIDE** | `GWCD4i.DE_stats.by_guide.h5mu` | guía × condición |
| **DE_BY_DONOR** | `GWCD4i.DE_stats.by_donors.h5mu` | par-de-donantes × perturbación × condición |
| **GUIDE_KD** | `guide_kd_efficiency.suppl_table.csv` | guía × condición |

## Claves y joins principales

- **Gen** es la entidad de referencia central (`gene_id` = Ensembl `ENSG…`, `gene_name` = símbolo).
  Aparece en dos roles: *gen perturbado* (diana de la guía) y *gen medido* (columna de la matriz de expresión / `.var`).
- **SGRNA.target_gene_id → GENE.gene_id**: cada guía apunta a un gen (ojo: `designed_target_gene_id`
  puede diferir de `target_gene_id` por curación post-hoc; hay ~1–2 guías por gen).
- **CELL.guide_id → SGRNA.sgRNA** (valor especial `multi-guide` si se detectó más de una guía).
  **CELL.lane_id → SAMPLE** (una lane 10x = un output de cellranger = una library).
- **PSEUDOBULK** = agregación de CELL por la clave compuesta `(guide_id, donor_id, culture_condition)`.
- **DE_RESULT** = agregación por `(target_contrast = gene_id, culture_condition)`; junta las `n_guides` guías del gen.
  `DE_stats.suppl_table.csv` es exactamente el `.obs` de este objeto en forma tabular.
- **DE_VALUE** (en `.layers`: `log_fc`, `zscore`, `adj_p_value`, …) es la relación N:N entre
  **DE_RESULT** (obs) y **GENE** (var): para cada perturbación×condición, un vector sobre los genes medidos.
- **DE_BY_GUIDE** y **DE_BY_DONOR** son la misma estructura que DE_RESULT pero desagregada
  (por guía individual, o por par de donantes) — sirven para métricas de reproducibilidad
  (`guide_correlation_*`, `donor_correlation_*`) que viven en `DE_RESULT.obs`.

### Nota sobre IDs de donante
Las etiquetas cortas `D1..D4` (nombres de archivo cell-level) se resuelven al `donor_id`
canónico `CE…` vía `sample_metadata` (`cell_sample_id` codifica `run_D#_condición`).
Las modalidades de `DE_stats.by_donors.h5mu` usan los IDs `CE…` unidos por `_`.
