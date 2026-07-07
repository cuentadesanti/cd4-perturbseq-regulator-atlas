# Análisis de edges (Modelo 1) — ¿resultado fuerte o bonus?

## Qué se analizó

La red de efectos con incertidumbre `docs/tables/robust_edges.csv`: **2,470 edges** robustos
(`P(|efecto|>1.5×)>0.8`) de **6 reguladores** × **2 condiciones**, leídos por slice del
`DE_stats.h5ad` remoto. Se resumió por regulador y por gen downstream, con dirección y convergencia.

## ¿Se ven útiles?

**Sí como prueba de concepto, no como resultado fuerte.** Señales a favor:

- **Convergencia coherente**: 579 genes downstream son target de ≥2 reguladores. Los 6
  reguladores son co-miembros del **complejo SAGA** (TADA1/TADA2B/SGF29/SUPT20H) → que compartan
  targets es exactamente lo esperado, y **valida que el método recupera estructura biológica real**.
- **Dirección interpretable**: 83% de los edges son de activación, consistente con
  SAGA como coactivador de la transcripción.
- La magnitud de los efectos es modesta y bien acotada (|θ| mediana ≈ 0.92).

## ¿Por qué NO es un resultado fuerte (todavía)

- **Cobertura mínima**: 6 de 7,913 reguladores (0.08%), **seleccionados por el ranking EB**
  → muestra sesgada a un solo complejo, no representa el paisaje regulatorio.
- **Condiciones incompletas**: solo Stim48hr, Stim8hr (la demo tomó la
  condición pico por regulador; **falta Rest**).
- **Latencia remota**: escalar a ~150 reguladores son ~11 min de lectura del h5ad (4.5 s/fila medidos).
- **Semántica de la probabilidad**: `p_abs_effect_gt_1p5x` es **P(magnitud del efecto > 1.5×)**,
  NO P(existe una arista causal). No hay control de FDR a nivel de red ni de sparsity (spike-and-slab).

## Veredicto

Dejar como **bonus / proof-of-concept**: demuestra que el pipeline produce una red de efectos con
incertidumbre biológicamente sensata sin descargar 1.8 TB, pero la cobertura y la semántica no
sostienen todavía afirmaciones de red a escala. Para promoverlo a resultado fuerte: correr los ~150 top reguladores en
las 3 condiciones y añadir P(edge) con prior sparse.
