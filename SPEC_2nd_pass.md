# Spec para Claude Code — 2nd pass de community detection (validación + identidad)

Contexto: PR #23 shippeó el core (RMT-clean → Leiden×50 → consenso → gate 0.8). Resultado: 8 comunidades, 3 estables (n=177/132/87, comm 5/6/7), y las 6 subunidades SAGA co-clusterizan en comm 2 (n=526, s_c=0.56, sub-gate). Este 2nd pass hace la parte interpretativa/validación, en orden de dependencia. Todo laptop desde la DE heredada. Aditivo, sufijo _3106, nunca destructivo.

## Task A — Null de permutación (guardrail; hacelo PRIMERO)
El λ+=1.49 es cota inferior porque las columnas-gen no son i.i.d., así que "336 señales" probablemente sobre-cuenta. Fijá el borde de ruido con un null EMPÍRICO:
- Sobre la matriz fingerprint real X (3106×6000), generá ≥100 versiones null por **phase-randomization** (o permutación independiente de cada columna), computá corrcoef, quedate con el autovalor máximo de cada null → distribución del borde empírico λ+_emp (percentil 95).
- Reportá: n_señales sobre λ+_emp (vs 336 sobre el MP teórico). Ese es el número honesto para citar.
- ADEMÁS un **community null**: ¿la modularidad de la partición consenso supera la de particiones sobre matrices null? Reportá z-score de modularidad.
- Output: `docs/tables/operator_community_null_3106.csv` (edge teórico vs empírico, n_señales cada uno, modularity real vs null z).
- Criterio de aceptación: si n_señales_empírico << 336, se DOCUMENTA (el core ya lo flagueó como cota); no invalida las comunidades, recalibra el conteo de rank.

## Task B — K562 sobre el MÓDULO SAGA (el resultado; el de mayor payoff)
OJO: la unidad NO es "comm 2" entera (526 regs, mayoría no-SAGA). La unidad es el **set de las 6 subunidades SAGA que co-clusterizan en comm 2** {SUPT20H, TAF6L, TADA2B, SUPT7L, USP22, SGF29} — el sub-cluster coherente. (KAT2B NO: el clustering la puso en comm 0, así que no es miembro del módulo.) El test cross-cell-type reemplaza el z de un solo gen (SUPT20H) por el módulo n>1.
- Fuente K562: `docs/tables/k562_concordance.csv` (cols: symbol, pearson_z, well_powered, donor_robust, program, min_kd_rank...). Es concordancia por-regulador CD4↔K562, ya calculada.
- Test del módulo: de las 6 subunidades SAGA, ¿cuántas están en el shared set de K562 y con qué pearson_z? El resultado n>1: si ≥N subunidades muestran concordancia universal en K562 (no solo SUPT20H), el programa SAGA es cross-cell-type sobre el módulo, no sobre una ancla.
- Validación no-circular del módulo (la versión fuerte): computá la correlación regulador-regulador EN K562 (desde el bulk Replogle ya en disco, mismo pipeline que CD4) restringida a las subunidades SAGA, y testeá su modularidad interna vs un null de reguladores random del mismo tamaño. Un módulo que reaparece en K562 sobre >1 miembro es el resultado.
- Output: `docs/tables/k562_saga_module_3106.csv` (subunidad, en_k562_shared, pearson_z, well_powered, class universal/T-spec) + un statement de modularidad del módulo vs null.
- Escribilo HONESTO: si solo SUPT20H+1 reaparecen, es "el módulo se reduce a ~2 miembros en K562" — reportá el número real, no fuerces.

## Task C — Enriquecimiento de complejos por comunidad (identidad)
Da identidad biológica a las comunidades estables (incl. la n=87 sin nombre) de forma no-circular.
- Para cada comunidad (las 3 estables 5/6/7 + comm 2 SAGA), test de enriquecimiento hipergeométrico contra CORUM (complejos) y opcionalmente STRING/GO, con corrección FDR (BH).
- Fuente de anotaciones: bajar CORUM (allComplexes.txt) — si no hay red, marcá el task como bloqueado y pedí el archivo. NO inventes membresías.
- Output: `docs/tables/operator_community_enrichment_3106.csv` (community, complejo/término, overlap, p, q_FDR). Solo reportá q<0.05.
- Esto es lo que convierte "comm 7, n=87, 45% donor-robust, sin nombre" en "comm 7 = complejo X".

## Task D — Convergencia con CP / SVD (cross-check, el más barato)
¿Las comunidades contradicen o confirman §4 (factores CP) / los programas SVD?
- Cross-tab: para cada factor CP (loadings en `operator_cp_factors_3106.csv`), ¿sus top reguladores caen en una comunidad, o se reparten?
- Output: `docs/tables/operator_community_cp_overlap_3106.csv`. Es un cross-check de consistencia, no un resultado nuevo.

## Orden y escritura
1. A (null) → 2. B (K562-SAGA) → 3. C (enrichment) → 4. D (convergencia).
- A es guardrail de todo; B es el resultado; C/D son identidad + consistencia.
- Escribir en OPERATOR_ANALYSIS.md §Step 6 2nd pass, con los números reales. Honestidad sobre resultados negativos (módulo inmune fragmentado ya establecido; si SAGA se reduce en K562, decirlo).
- Aditivo, sufijo _3106, PR draft. No tocar el core de PR #23.
