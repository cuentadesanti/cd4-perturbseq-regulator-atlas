# EDA 80/20 — CD4+ T cell Perturb-seq

## Scope

Análisis con solo las **tablas suplementarias** (~15 MB). **No** se cargó ningún `.h5ad`/`.h5mu` (1.8 TB).

| Archivo | Filas | Uso |
|---|---|---|
| `DE_stats.suppl_table.csv` | 33,983 | tabla principal (la señal) |
| `sgrna_library_metadata.suppl_table.csv` | 26,504 | librería de guías |
| `sample_metadata.suppl_table.csv` | 12 | diseño experimental |

Regenerable: `python scripts/eda.py` → figuras en `docs/figures/`, tabla en `docs/tables/`.

## Unit of analysis

**1 fila = gen perturbado × condición de cultivo** (`target_contrast` × `culture_condition`).
Rest / Stim8hr / Stim48hr.

## Key quality filters

| Filtro | Disponible aquí | Fuente |
|---|---|---|
| `ontarget_significant` (KD efectivo, 10% FDR) | ✅ CSV | `DE_stats.suppl_table.csv` |
| `offtarget_flag` (posible off-target) | ✅ CSV | idem |
| reproducibilidad **cross-condición** (proxy) | ✅ derivable | min/max de `n_downstream` entre condiciones |
| `single_guide_estimate` (2 guías concordantes) | ❌ | `DE_stats.h5ad` `.obs` |
| `guide_correlation_all` (cross-guide) | ❌ | `DE_stats.h5ad` / `by_guide.h5mu` |
| `donor_correlation_hits_mean` (cross-donante) | ❌ | `DE_stats.h5ad` / `by_donors.h5mu` |

## Main findings

1. **Los efectos de DE son heavy-tailed.** Mediana **2 DEGs**, media 60.5 (engañosa), 15.4% sin efecto,
   1.5% son hubs (>1000 DEGs). → resume con **percentiles y rankings**, no media.
2. **El knockdown gatea la mayor parte de la señal.** 62% de contrastes con KD on-target significativo;
   estos concentran el **85%** de todos los trans-efectos. Filtrar por `ontarget_significant` sube mucho
   la densidad de señal (no es prueba causal: puede haber KD real no detectado por baja expresión basal).
3. **Las células estimuladas muestran efectos más amplios.** Media de DEGs: Rest 53.1 · Stim8hr 68.9 · Stim48hr 59.4.
4. **Los top hubs son plausibles**, enriquecidos en señalización de células T (CD3E/D/G, LAT, ZAP70, PLCG1,
   LCP2, VAV1) — sugiere que la pantalla captura señal biológica interpretable.
5. **La librería tiene cobertura ~2 guías/gen** (12,440 de 12,654 genes) → replicación interna.
6. **Reguladores robustos ≠ hubs crudos.** Al exigir estabilidad cross-condición + KD signif. + sin off-target,
   suben reguladores de **cromatina/transcripción** consistentes en las 3 condiciones (TADA2B, TADA1, SGF29,
   SUPT20H — complejo SAGA; ELAVL1, NFRKB), distintos del top crudo dominado por señalización TCR específica de Stim8hr.

## Figuras

`docs/figures/`
- `01_distribution_n_total_de_genes.png` — cola larga
- `02_degs_by_condition.png`
- `03_top_hubs_by_condition.png`
- `04_ontarget_vs_downstream.png` — eje `kd_strength = −ontarget_effect_size`
- `05_guides_per_gene.png`
- `06_reproducibility_vs_effects.png` — reproducibilidad cross-condición vs magnitud

## Tabla accionable

`docs/tables/top_robust_regulators.csv` (top 30). Score **solo con columnas del CSV**:

```
robust_score = log1p(n_downstream)
             · ontarget_significant
             · (0.6 si offtarget_flag else 1.0)
             · (0.5 + 0.5 · reproducibilidad_cross_condición)
             · (n_signif_conditions / 3)
```

## Practical next steps

1. **Ranking robusto definitivo**: reforzar `robust_score` con `single_guide_estimate` +
   `donor_correlation_hits_mean` desde `GWCD4i.DE_stats.h5ad` (17 GB).
2. **Matriz downstream a nivel de gen**: cargar los `.layers` de `DE_stats.h5ad` (log_fc/zscore/padj)
   para los reguladores top, ya filtrados.
3. **Grafo regulatorio por condición**: construir la red regulador → downstream por condición y
   comparar Rest vs Stim para identificar reguladores context-specific.
