# Side analysis — transcriptional fingerprints / perturbation program similarity

## Qué pregunta responde

El ranking EB y las auditorías dicen **quién** es un regulador fuerte. Este análisis responde otra
pregunta: **qué programa transcriptómico induce cada perturbación, y qué reguladores se parecen
entre sí**. Pasa de una *lista de genes* a un *mapa de programas*.

La **huella** de una perturbación = su vector de efectos sobre los ~10k genes medidos. Dos
reguladores con huellas parecidas actúan sobre el mismo programa. Sobre esa matriz se corre PCA,
similitud coseno y clustering. No reemplaza el ranking core, no agrega modelos pesados, no descarga
el h5ad completo.

## Datos y método

Matriz `layers/zscore` del h5ad remoto (densa, contigua): filas = perturbación×condición, columnas
= genes medidos. Se lee **sólo el panel** (~200 filas) por slice remoto, cacheado en
`data/cache/panel_zscore_*.npy`. (`--matrix log_fc` usa el cache local completo de 1.4 GB.)

**Panel BALANCEADO** — clave para no sesgar el mapa. Si tomas sólo el top-EB, sale puro
cromatina/SAGA/Mediador. En cambio:

| fuente | n | criterio |
|--------|---|----------|
| global            | 75 | top por regpower EB, clase `global` |
| context-specific  | 75 | top por regpower EB, clase `condition-specific` |
| promoted          | 25 | promovidos por la auditoría de reproducibilidad (rank_shift↑) |
| demoted           | 25 | demotados por la auditoría (rank_shift↓) |

Cada regulador entra en su condición pico. Columnas: **top-2000 genes por varianza** en el panel.
Se **estandarizan las columnas** (z-score por gen) antes de PCA/similitud, para que genes
comparables pesen igual.

Ejecutar: `make fingerprints`  (o `python scripts/analyze_fingerprints.py --n 200 --matrix zscore --top-genes 2000`).

## Resultados

### 1. El espacio latente captura *programas*, no *magnitud*

PCA sobre las huellas. Si PC1 sólo reordenara por "tamaño del efecto", el análisis no añadiría nada
sobre el ranking. **No es el caso:**

- `spearman(|PC1|, n_downstream) = 0.25` — correlación débil. PC1 **no** es magnitud.
- PC1 (28.9% var) es el eje de **señalización TCR / activación proximal**; PC2 (7.4%) separa los
  **coactivadores de cromatina** (SAGA/Mediador, arriba) del resto.
- Los top-loadings son inmunológicamente coherentes (citoquinas efectoras, naive/memoria,
  polarización Th, interferón).

→ **El ranking captura fuerza; el mapa latente captura la identidad del programa.** Complementarios.

### 2. El espacio recupera estructura biológica conocida (validación por permutación)

¿Los miembros de un complejo tienen huellas más parecidas entre sí que gene-sets aleatorios del
panel? (coseno intra-complejo vs null por permutación, N=5000):

| Complejo  | n | coseno intra | null | z | p |
|-----------|---|--------------|------|-----|-------|
| **TCR**       | 5 | **0.916** | 0.004 | **11.2** | 0.0002 |
| **SAGA**      | 8 | **0.478** | 0.004 | **9.2**  | 0.0002 |
| **Mediador**  | 7 | 0.196 | 0.007 | 3.2 | 0.013 |

Los tres son significativos. TCR es casi colineal (0.92): CD3E/ZAP70/LAT/LCK/PLCG1 inducen
prácticamente la misma huella. SAGA muy cohesivo. Mediador más difuso pero significativo.

### 3. Vecinos transcriptómicos — la vista demoable

Nearest neighbors por similitud coseno de huella (`fingerprint_neighbors.csv`, endpoint
`/programs/neighbors/{gene}`):

- **ZAP70** → LAT (0.95), CD3E (0.95), PLCG1 (0.95), LCK (0.84) — módulo TCR puro.
- **TADA2B** → TAF6L (0.94), SUPT20H (0.92), TADA1 (0.88) — módulo SAGA puro.
- **MED12** → MED19 (0.78), SUPT20H (0.75), TADA2B (0.75) — Mediador + SAGA (coactivadores).
- **SGF29** → CCNC (0.66) — SAGA más cercano a un miembro de Mediador.

El atlas ya no sólo rankea reguladores: **organiza perturbaciones por la huella que inducen** y
permite preguntar "¿a quién se parece este gen?".

## Salidas

**Figuras** (`docs/figures/`): `20_fingerprint_pca_by_condition.png`,
`21_fingerprint_pca_by_regulator_class.png` (principal), `22_fingerprint_similarity_heatmap.png`,
`23_fingerprint_neighbor_network.png`.

**Tablas** (`docs/tables/`): `fingerprint_panel.csv`, `fingerprint_pca_scores.csv`,
`fingerprint_neighbors.csv`, `fingerprint_similarity_edges.csv`, `fingerprint_clusters.csv`,
`fingerprint_pc_loadings.csv`, `fingerprint_complex_validation.csv`, `fingerprint_summary.json`.

**API** (`/programs/*`): `GET /programs/pca`, `GET /programs/neighbors/{gene}`,
`GET /programs/clusters`. **UI**: pestaña **Programas** (scatter PCA + buscador de vecinos).

## Limitaciones (honestas)

- **Una huella por regulador (condición pico)** → el mapa mezcla condiciones. La figura por
  condición ayuda a leer el efecto de contexto, pero un PCA separado por condición sería más limpio.
- **zscore = log_fc/lfcSE** ya pondera por error, pero se estandariza por columna encima; es una
  elección razonable, no la única.
- **PCA/SVD lineal**, no UMAP: más transparente y defendible, pero no captura no-linealidades.
- **Exploratorio**: es *perturbation program similarity*, **no** descubrimiento de vías causales ni
  de programas celulares nuevos. Los complejos usados son anotación conocida, no descubrimiento.
- **Panel de 200**, no los 33 983 contrastes: suficiente para el mapa, no exhaustivo.

## Veredicto

**Resultado fuerte.** Tres afirmaciones cuantificadas y defendibles: (1) la estructura latente es
programa, no magnitud (|PC1|~n_downstream = 0.25); (2) el espacio recupera complejos conocidos con
significancia por permutación (TCR z=11, SAGA z=9, Mediador z=3); (3) los vecinos transcriptómicos
son biológicamente coherentes y navegables desde la UI. El atlas pasa de "ranking de reguladores" a
"mapa de programas regulatorios".
