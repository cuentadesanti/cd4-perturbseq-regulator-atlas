# Robust regulators of CD4+ T cell programs — Perturb-seq

> Hackathon submission · genome-scale CRISPRi Perturb-seq en células T CD4+ primarias (Marson Lab, 2025)

![pipeline](docs/figures/00_pipeline_overview.png)

## Pregunta

¿Qué genes son **reguladores robustos** de los programas de células T CD4+? El reto: la señal
de un screen genome-scale está dominada por ruido y por unos pocos hubs. Queremos **separar señal
de ruido con incertidumbre** y priorizar reguladores por efecto **grande y reproducible**, no por
conteos crudos ni `adj_p_value < 0.1`.

## Qué hace este repo

Un producto de investigación reproducible, **consciente de memoria y cómputo** (el dataset completo
son 1.8 TB; el disco de trabajo tiene ~10 GB). Todo el core corre **solo con las tablas CSV
suplementarias (~15 MB)**:

1. **Pipeline de descarga selectiva** desde el bucket público S3 (`scripts/download.sh`).
2. **Modelo de datos** del dataset (`docs/DATA_MODEL.md` + [artifact interactivo](docs/data-model.html)).
3. **EDA 80/20** (`scripts/eda.py`) — distribución de efectos, calidad de KD, hubs, reproducibilidad.
4. **Modelo 2 · empirical-Bayes** (`scripts/model_hubs.py`) — ranking de reguladores con incertidumbre.
5. **Modelo 1 · red de efectos con incertidumbre** (opcional, `scripts/model_edges.py`) — lee por *slice*
   el `.h5ad` de 17 GB desde S3 **sin descargarlo**.

## Dataset

**Primary Human CD4+ T Cell Perturb-seq** · [CZI Virtual Cells Platform](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq)
· bucket público `s3://genome-scale-tcell-perturb-seq/marson2025_data/` · 4 donantes × 3 condiciones
(Rest / Stim8hr / Stim48hr). El core usa 3 CSV (33,983 contrastes de DE, 26,504 guías, 12 muestras).

## Hallazgos principales

- **Efectos heavy-tailed**: perturbación mediana = 2 DEGs, 15% sin efecto, pero 1.5% son hubs (>1000 DEGs).
  → resumir con percentiles y rankings, no con la media.
- **El knockdown gatea la señal**: los contrastes con KD on-target significativo (62%) concentran el **85%**
  de todos los trans-efectos. Filtrar por `ontarget_significant` es el primer paso obligado.
- **Ranking robusto ≠ hubs crudos**: el modelo EB surface maquinaria de **cromatina/transcripción**
  (complejo SAGA: TADA1/TADA2B/SGF29/SUPT20H · Mediador: MED12/CCNC · KDM1A, SETD2) — reguladores con
  efecto grande **y** estable entre condiciones, por encima de los hubs de señalización TCR específicos de Stim8hr.
- **Red de efectos con incertidumbre (bonus)**: ~2,470 edges robustos (`P(|efecto|>1.5×)>0.8`, es decir
  probabilidad de que la *magnitud* del efecto supere 1.5×, no de que exista una arista causal) para los top
  reguladores, extraídos del h5ad remoto sin descargarlo.

Detalle: [`docs/report.md`](docs/report.md) · [`docs/EDA.md`](docs/EDA.md) · [`docs/MODELING.md`](docs/MODELING.md).

## Validación del ranking

Antes de creerle al ranking, lo auditamos (`scripts/audit_ranking.py`, sin deps nuevas).

**Naive hubs vs quality-aware regulators.** Rankear por `n_downstream` crudo premia hubs que no
sobreviven a los controles: de los 30 hubs crudos top, **2 caen por el gate de KD** (sin knockdown
validado) y **~9 se demotan** por ser condition-specific. El ranking EB surface reguladores con
efecto grande **y** estable. La estabilidad se midió por bootstrap (B=200): es *moderada* —léelo
como un conjunto de reguladores robustos, no como un orden exacto (`stability_frequency` por gen).

**Global versus context-specific regulators.** Separando por
`condition_specificity = max/sum de n_downstream`: los **globales** (SGF29, TADA2B, SUPT20H…) son
maquinaria de cromatina/transcripción activa en todas las condiciones; los **context-specific**
(ZAP70, LCK…) son señalización TCR, activa solo bajo estímulo. Ambas clases son biología real; la
distinción evita confundir un regulador universal con uno de contexto. Ver
`docs/tables/top_global_regulators.csv` y `top_condition_specific_regulators.csv`.

## Cómo reproducir

```bash
# 1. entorno
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2. datos (solo las tablas, ~15 MB)
scripts/download.sh tables

# 3. pipeline completo (EDA + modelo + reporte), verificado
make all
```

Targets: `make eda` · `make model` · `make report` · `make all` · `make clean`.
Opcional (remoto, requiere `pip install h5py s3fs fsspec`): `make spike` · `make edges`.

## Outputs

| Archivo | Qué es |
|---|---|
| `docs/report.md` | reporte judge-facing consolidado |
| `docs/tables/top_regulators_for_review.csv` | top 30 reguladores (judge-facing) |
| `docs/tables/hub_ranking_bayes.csv` | ranking EB completo (todos los genes) |
| `docs/tables/robust_edges.csv` | red de efectos con incertidumbre (bonus, Modelo 1) |
| `docs/figures/*.png` | EDA + ranking + overview |
| `docs/data-model.html` | explorador interactivo del modelo de datos |

## Limitaciones

- **Nomenclatura honesta**: los modelos son **empirical-Bayes / pseudo-bayesianos**, no NB jerárquico
  completo ni MCMC (sin PPL, sin random effects formales).
- `xcond_reproducibility` es una **feature exploratoria** (estabilidad cross-condición); **no** sustituye
  la reproducibilidad cross-donor/cross-guide, que vive en `DE_stats.h5ad` (`single_guide_estimate`,
  `donor_correlation_hits_mean`) — marcadas como `NA` en la tabla de review.
- **Modelo 1 es opcional**: el acceso remoto por slice es viable (~4.5 s/fila, medido) pero latency-bound;
  el entregable oficial se sostiene solo con el core local.

## Submission summary

Producto reproducible que convierte una matriz de DE de 1.8 TB en un **ranking de reguladores robustos
con incertidumbre**, ejecutable en una laptop con ~10 GB de disco usando solo 15 MB de datos, más una
**red de efectos con incertidumbre** (uncertainty-aware effect network) opcional leída en streaming del
h5ad de 17 GB sin descargarlo.
`make all` reproduce todo el core; ver `docs/report.md`.

---
Datos: CZI Virtual Cells Platform · Marson Lab 2025 · preprint biorxiv `10.64898/2025.12.23.696273`.
Código de análisis original: [emdann/GWT_perturbseq_analysis_2025](https://github.com/emdann/GWT_perturbseq_analysis_2025).
