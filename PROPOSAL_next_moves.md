# Qué hacer con las palancas abiertas — síntesis de tres frentes de research

Deep research en paralelo sobre las tres preguntas que dejaste abiertas (el R² absoluto, datasets con metadata como regresores, y el rank/clustering). Cada frente hizo research con literatura verificada. Abajo la síntesis con una recomendación única y priorizada.

---

## TL;DR — la recomendación en una línea

**Un solo paquete de trabajo, dos piezas, ambas en laptop desde la DE heredada:**
1. **Community detection RMT-limpio del operador** (clustering de reguladores) → convierte las anclas n=1 en módulos n>1 estables, y de paso barre el confound de potencia. **Mata el punto más fuerte del crítico (#5).**
2. **Bi-cross-validation del rank** → convierte el hedge "~7 vs ~86" en la descomposición defendible **"7 de 86 direcciones son transferibles"**, con barras de error. **Mata #2 y #4.**

Y una decisión de encuadre que **no cuesta código**: el R²≈0.07 **ya le gana al mean-predictor**, que es el baseline real del campo — no hay que perseguir un R² más alto ni meter un modelo más grande.

---

## Frente 1 — El R² absoluto: KEEP + reframe, NO demote, NO modelo más grande

**El hallazgo que da vuelta la crítica.** El crítico dijo "le ganás a un baseline sub-media (persistence −0.305), y a R²=0.074 apenas superás el mean-predictor trivial". **Confunde dos baselines.** En la literatura de predicción de perturbaciones, el baseline duro NO es persistence — es **el mean across training perturbations**, y es la barra que los foundation models NO logran superar:

- **Ahlmann-Eltze, Huber & Anders 2025** (Nat Methods, `10.1038/s41592-025-02772-6`): cinco foundation models + GEARS, comparados contra baselines deliberadamente simples. Para perturbaciones no vistas, los modelos NO le ganan a predecir la media.
- **Wei 2025 scPerturBench** (`10.1038/s41592-025-02980-0`) y **Peidli 2024** (E-distance, `10.1038/s41592-023-02144-y`): mismo mensaje — el campo se sostiene sobre baselines simples, y superarlos es el logro real.
- **Wong 2025** (`10.1093/bioinformatics/btaf317`): refuerza el mismo punto sobre modelos lineales vs complejos.

**En esa luz, R²=+0.074 > 0 (mean-predictor) es un resultado, no un fracaso** — es exactamente lo que los modelos grandes no consiguen. El error de la crítica es leer persistence (−0.305) como "nuestro baseline"; nuestro resultado bate DOS baselines (media Y persistence), y solo mostramos el segundo.

**Qué hacer (paquete de evaluación, todo laptop, DE heredada):**
- **P1 — Escalera de baselines + distribución de R² per-entidad.** Declarar explícito el piso del mean-predictor (R²=0) al lado de persistence (−0.305), y reportar el R² por-gen / por-regulador como distribución, no un número global. Muestra que *algunas* entidades se predicen bien y otras no — más honesto y más informativo que el R² agregado.
- **P2 — Margen cross-validado con barras de error + null de permutación.** Reemplazar el seed-0 único por K folds / múltiples seeds → margen ± SD, y un null de permutación para la significancia. Convierte "el signo es positivo" en "+0.38 ± algo, p<...".
- **P5 — Techo de reproducibilidad donor-split → R² normalizado.** Usar la reproducibilidad entre donantes como techo alcanzable y normalizar el R² contra él (no podés predecir mejor que el ruido donante-a-donante). Da un "% del techo alcanzable capturado", mucho más interpretable que 0.074 crudo.

**Recomendación del frente 1:** hacer P1+P2+P5 como UN paquete de upgrade de evaluación. Mata tres puntos del crítico a la vez (#1 R² casi-nulo, #2 single-seed, #1 baseline bajo). **NO perseguir el R² con un predictor más grande** — la literatura predice que sería overfitting.

---

## Frente 2 — Datasets con metadata como regresores: el operador inductivo (predice reguladores NUNCA vistos)

Esta es tu pregunta "investiga de datasets con más metadata como regresores", y es la palanca conceptualmente más novedosa — ataca el #3 del crítico (la tarea held-out es más fácil de lo que suena).

**El problema exacto.** Hoy el completion es low-rank puro sobre el tensor de efectos, así que solo puede extrapolar una **fibra de condición** de un regulador que ya viste en otras condiciones. NO puede predecir un regulador nunca medido. Para eso hace falta hacer el modelo **inductivo**: darle features del regulador y del gen para que pueda scorear un regulador solo desde sus anotaciones.

**El método (verificado):**
- **Inductive Matrix Completion — Natarajan & Dhillon 2014** (Bioinformatics, `10.1093/bioinformatics/btu269`, "Inductive matrix completion for predicting gene–disease associations"). VERIFICADO. Exactamente nuestro caso: completa una matriz incompleta usando features en filas y columnas, y predice filas/columnas nuevas desde sus features. El paper original es literalmente gene-disease; nosotros seríamos regulator-gene.
- Familia relacionada (teoría de IMC, collective/coupled MF, graph-regularized completion): **DOI unverified** — el subagente se atascó verificando estos; hay que confirmarlos antes de citar. El de Natarajan-Dhillon alcanza como ancla.

**Las features usables como regresores (lado regulador):**
- **TF census — Lambert et al. 2018** (Cell, "The Human Transcription Factors") — clase de TF, dominio de unión.
- **CORUM** — pertenencia a complejos proteicos (SAGA, Mediator son complejos → features directas).
- **STRING / BioGRID** — PPI, para regularización por grafo.
- Expresión basal, anotación GO, listas de reguladores de cromatina.
- Lado gen: mismas fuentes proyectadas a genes.

**El payoff novedoso.** Un operador que predice un regulador genuinamente out-of-sample **desde sus anotaciones** convierte el resultado hedgeado "fibra de condición out-of-panel" en predicción real de regulador-nunca-visto — mata #3 de raíz. Test honesto: **leave-regulator-out** (sacás un regulador entero, lo predecís solo desde sus features).

**Caveat honesto.** A escala ~3106 reguladores con DE heredada, IMC puede funcionar pero el techo lo pone cuánta señal llevan las features (si SAGA-membership predice el efecto trans, gana; si no, el R² inductivo va a ser bajo — posiblemente más bajo que el transductivo actual). Es la apuesta más novedosa **y** la más incierta. Yo la haría como experimento separado, no como reemplazo del resultado actual.

*(Nota: este frente lo tuve que sintetizar yo — el subagente verificó los DOIs de IMC pero se colgó en el loop de verificación antes de escribir el memo. El DOI de Natarajan-Dhillon está verificado; el resto de la familia IMC hay que verificarlo antes de citar.)*

---

## Frente 3 — Rank + clustering: el FLAGSHIP que mata el n=1

Tu "¿no clustereamos nada?!" es correcto y es la observación más filosa. Hicimos SVD/CP y una asignación supervisada por coseno (el "fingerprint"), pero **nunca clustering no-supervisado / community detection del operador mismo.** Y los nombres de programa descansan en anclas únicas (SAGA=SUPT20H, T-specific=DOCK2).

**P3 (FLAGSHIP) — Community detection RMT-limpio sobre el grafo de reguladores:**
- Construir la matriz de correlación regulador-regulador 3106×3106 desde el fingerprint aplanado (2000 genes × hasta 3 condiciones).
- **NO correr modularity directo** — MacMahon & Garlaschelli 2015 (`10.1103/PhysRevX.5.021006`) muestran que ese null es inconsistente para matrices de correlación. En su lugar: **restar la banda de ruido Marchenko-Pastur + el modo global** de la correlación. **El modo global ≈ nuestra dirección de fuerza-de-KD/potencia** → limpiar por RMT es a la vez el null correcto Y un barrido automático del confound de potencia (el problema del spearman=0.20). Dos pájaros de un tiro.
- Correr **Leiden** (Traag 2019, `10.1038/s41598-019-41695-z`) sobre la matriz limpia, barrido de resolución + muchas semillas.
- Envolver en **consenso** (Lancichinetti-Fortunato 2012, `10.1038/srep00336` / Monti 2003, `10.1023/A:1023949509487`), quedarte solo con comunidades de estabilidad de co-asignación ≥ 0.8 (lógica de stability-selection, Meinshausen-Bühlmann 2010, `10.1111/j.1467-9868.2010.00740.x`).
- **Precedente directo que hay que citar:** Replogle 2022 (`10.1016/j.cell.2022.05.013`) hizo exactamente esto en K562 — clusterizó perturbaciones en módulos y recuperó complejos (propuso C7orf26 = INTS15). Es el movimiento n>1 que no hicimos, en el dataset hermano.

**El payoff.** SAGA pasa a ser una comunidad estable {SUPT20H, TADA2B, TAF6L, …} cuya co-membresía la maneja el fingerprint trans **aunque** TADA2B/TAF6L sean cis-inflados individualmente; el módulo inmune-actina pasa a ser {NCKAP1L, WASF2, DOCK2, ARHGAP30, AHR/ARNT, …}. **La comunidad, no el ancla, es el objeto** — mata #5. Y el test de K562 se corre **sobre módulos, no sobre anclas** (validación no-circular: modularidad del módulo en K562 vs null de miembros permutados).

**P2 (al lado) — Bi-cross-validation del rank (Owen-Perry 2009, `10.1214/08-AOAS227`):**
- Reemplazar el holdout seed-0 único por BCV (holdout de bloque fila×columna, grilla de folds) → **distribución del rank óptimo con CI**, no un punto.
- Reportar los dos números juntos: rank predictivo BCV (~7) vs rank de señal GD/DDPA (~86). **El gap ES el hallazgo:** solo ~7 de ~86 direcciones de señal generalizan; el resto es idiosincrático del regulador — lo que explica honestamente el techo R²≈0.07 (es mayormente irreducible, no una falla del modelo).
- Sub-test de robustez: re-correr held-out sobre las fibras Rest y Stim8hr (no solo Stim48hr) + barrer holdout 10/20/40%. Si el rank óptimo se mueve según qué fibra sacás, el ~7 es config-específico y lo decimos.

**P1 (prerequisito casi gratis) — Piso de ruido RMT sobre el operador z-score:**
- Como las entradas son z = efecto/SE, bajo el null cada entrada es ≈N(0,1): **el ruido ya está estandarizado a σ=1**, input casi ideal para RMT. Overlay del threshold Gavish-Donoho (2014, `10.1109/TIT.2014.2323359`) con σ=1 + Dobriban-Owen DDPA (2019, `10.1111/rssb.12301`). Reemplaza "≈86 a ojo" por "k componentes sobre un piso σ=1 falsable".

**Recomendación del frente 3:** P3 primero (mata el n=1), P2 al lado (arregla la narrativa de rank), P1 como prerequisito barato. **NO** perseguir el R²≈0.07 con un predictor más elaborado.

---

## Cómo encaja todo — el plan priorizado

**Tier 1 (hacelo, mata 4 puntos del crítico, todo laptop):**
1. **P3 community detection RMT-limpio** → módulos n>1, barre confound de potencia (mata #5, mitiga #7).
2. **P2 BCV rank** → "7 de 86 transferibles" con barras de error (mata #2, #4).
3. **Reframe de evaluación (P1+P5 del frente 1)** → declarar el piso mean-predictor + normalizar al techo donor-split (mata #1).

Estos tres se refuerzan: el mismo grafo RMT-limpio del P3 da el modo-global que documenta el confound; el mismo piso σ=1 sirve al P1 de rank y al P1 de baselines.

**Tier 2 (apuesta novedosa, experimento separado):**
4. **Operador inductivo (IMC)** con features de regulador (Lambert TF census + CORUM + STRING) → predicción de regulador-nunca-visto (mata #3). Más incierto: el R² inductivo puede salir bajo. Vale como sección nueva "hacia predicción inductiva", no como reemplazo.

**Lo que NO hay que hacer** (consenso de los tres frentes): perseguir el R² con un modelo más grande. La literatura (frente 1) y la descomposición de rank (frente 3) coinciden en que el techo es mayormente irreducible bajo este holdout — subir R² sería overfitting.

**Novedad honesta.** El paquete Tier 1 responde tu "sé novedoso" sin inventar: es la primera vez que (a) se hace community detection no-supervisado y estabilidad-gated del operador (vs el fingerprint supervisado), (b) la limpieza RMT sirve doble como null correcto Y de-confounder de potencia, y (c) el rank se reporta como descomposición transferible-vs-idiosincrático con CI. Todo sobre datos que ya están en el laptop.

**Límite honesto que hereda todo:** las cuatro propuestas usan los efectos DE heredados y NO modelan incertidumbre de guide-assignment (el gap conceptual reconocido). Un módulo construido sobre z heredado hereda ese gap. Nada de esto necesita los 1.8 TB cell-level; propagar guide-assignment sí lo necesitaría.
