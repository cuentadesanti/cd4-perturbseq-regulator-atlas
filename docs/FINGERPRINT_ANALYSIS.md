# Side analysis — mapa de programas regulatorios (huellas transcriptómicas)

## Qué pregunta responde

El ranking EB y las auditorías dicen **quién** es un regulador fuerte. Este análisis responde una
pregunta distinta: **qué tipo de efecto** produce cada perturbación. Pasa de una *lista de genes*
a un *mapa de programas*.

La idea: cada perturbación deja una **huella transcriptómica** = el vector de log fold-changes de
los ~10k genes medidos cuando apagas ese gen. Dos reguladores con huellas parecidas actúan sobre
el mismo programa. Con eso se puede clusterizar, hacer PCA y validar contra biología conocida.

## De dónde salen los datos (sin descargar 1.8 TB)

La matriz de huellas **ya existe** en el h5ad remoto como `layers/log_fc`: densa, contigua y sin
comprimir, shape **(33 983 perturbación×condición) × (10 282 genes medidos)**. Al ser contigua se
puede leer por range-reads baratos. Se cacheó completa una sola vez a `float32`:

- `scripts/fetch_fingerprint_matrix.py` → `data/cache/log_fc.f32.npy` (1.40 GB) + `obs` + `var`.
- `scripts/analyze_fingerprints.py` → tablas `docs/tables/fingerprint_*` y figuras `docs/figures/fingerprint_*`.

> Nota de ingeniería: bajar la **precisión** (float64→float32) reduce RAM/disco pero **no** la
> transferencia — los bytes en S3 son float64 y sólo se pueden leer tal cual. La palanca de red real
> es leer **menos filas**, no menos bits. Como es contigua, la matriz completa cabe en una descarga
> secuencial (~2.8 GB), así que se bajó todo y el análisis corre offline.

**Alcance del análisis:** 220 reguladores (top-200 por `n_downstream` ∪ miembros de complejos
conocidos), cada uno en su condición pico, filtrados a KD on-target significativo.

## Resultados

### 1. El espacio latente captura *programas*, no *magnitud* — el resultado clave

PCA sobre las huellas. Si PC1 sólo reordenara por "tamaño del efecto", el análisis no añadiría nada
sobre el ranking. **No es el caso:**

- `spearman(|PC1|, n_downstream) = 0.21` — correlación débil. PC1 **no** es magnitud.
- Los ejes son **inmunológicamente interpretables** (top loadings):
  - **PC1** (16.4%): señalización TCR / activación proximal — eje dominado por el módulo TCR.
  - **PC2** (5.6%): naive/memoria vs efector (TCF7, PLAC8, PI16 vs FN1, CCL1, CXCL8).
  - **PC3**: polarización Th / interferón (IL22, IL5, IL10, CXCL10, IFIT3).
  - **PC4**: citotóxico vs interferón (GNLY, EOMES vs IFIT3, IFI27).

→ **El ranking captura fuerza; el mapa latente captura la identidad del programa.** Son
información complementaria, no redundante.

### 2. El espacio recupera estructura biológica conocida — validación por complejos

Test: ¿los miembros de un mismo complejo tienen huellas más parecidas entre sí que gene-sets
aleatorios del mismo tamaño? (coseno intra-complejo vs null por permutación, N=5000).

| Complejo  | n | coseno intra | null | z | p |
|-----------|---|--------------|------|-----|-------|
| **TCR**       | 16 | **0.507** | 0.064 | **15.6** | 0.0002 |
| **Mediador**  | 10 | 0.211 | 0.064 | 3.9 | 0.0020 |
| **SAGA**      | 14 | 0.195 | 0.064 | 4.2 | 0.0014 |
| SWI/SNF   | 8  | 0.144 | 0.063 | 1.9 | 0.047 |

TCR es espectacular (z=15.6): sus miembros forman un gradiente limpio a lo largo de PC1
(CD3E, ZAP70, CD3D/G, LAT, LCK, VAV1, PLCG1, CARD11…). SAGA y Mediador se validan sólidamente y se
**entremezclan** en el mapa — coherente, ambos son coactivadores transcripcionales generales.
SWI/SNF es marginal pero ocupa una región propia.

### 3. Global vs context-specific — formalizado con el patrón downstream completo

Para 1 388 reguladores significativos en ≥2 condiciones: coseno de su huella **entre condiciones**.
Coseno alto = mismo programa siempre (global). Coseno bajo = el programa cambia con el contexto.

- **Globales** (coseno alto): TADA2B (0.65), SGF29 (0.62), SUPT20H (0.55), MED12 (0.54),
  TADA1 (0.50) + maquinaria de cromatina (SETDB1, KDM1A, NSD1, PCGF1, L3MBTL2). Su huella es
  estable en Rest/Stim8hr/Stim48hr.
- **Context-specific** (coseno bajo): **ZAP70 (0.29), LCK (0.28)** — la señalización TCR cambia
  radicalmente su efecto downstream entre reposo y activación, exactamente la biología esperada.

Esto formaliza la distinción SGF29-global vs ZAP70-específico que antes se aproximaba con conteos:
ahora usa el patrón downstream completo, no sólo *cuántos* genes cambian.

## Figuras

- `fingerprint_pca.png` — mapa PCA, complejos resaltados. **La figura principal.**
- `fingerprint_pc1_vs_magnitude.png` — el test "¿PC1 es sólo magnitud?" (no lo es).
- `fingerprint_dendrogram.png` — clustering jerárquico de las huellas.
- `fingerprint_context_specificity.png` — distribución global vs context-specific.

## Tablas

`fingerprint_regulators.csv` (cluster + scores PC), `fingerprint_similarity.csv` (matriz coseno),
`fingerprint_pc_loadings.csv` (genes que definen cada eje), `fingerprint_complex_validation.csv`,
`fingerprint_context_specificity.csv`, `fingerprint_summary.json`.

## Limitaciones (honestas)

- **Una huella por regulador (condición pico)** para el mapa principal → el PCA mezcla condiciones.
  El análisis global-vs-context sí es por-condición y lo compensa, pero un mapa por-condición
  separado sería más limpio.
- **Coseno sobre log_fc crudo**, sin ponderar por `lfcSE` (no se cacheó esa capa) → genes ruidosos
  pesan igual que los confiables. Un coseno ponderado por precisión sería más robusto.
- **El clustering jerárquico es grueso** al umbral usado (un cluster "activación" grande). El valor
  cuantitativo está en la validación por permutación y el PCA, no en el corte de clusters.
- **Enrichment interpretado a ojo** por los top-loadings, no con un test de over-representation formal.

## Veredicto

**Resultado fuerte, no bonus.** A diferencia del análisis de edges (6 reguladores, PoC), este cubre
220 reguladores y produce tres afirmaciones defendibles y cuantificadas: (1) la estructura latente
es programa, no magnitud; (2) el espacio recupera complejos conocidos con significancia por
permutación; (3) global-vs-context queda formalizado con el patrón downstream completo. El proyecto
deja de ser "ranking de reguladores" y se vuelve "mapa de programas regulatorios".
