# Modelado — de la matriz de DE a reguladores con incertidumbre

Dos modelos pequeños, **dependency-light** (solo `scipy` + `statsmodels`), que separan señal
de ruido con incertidumbre en vez de rankear por conteos crudos y `adj_p_value < 0.1`.

> **Nomenclatura honesta:** ambos son **empirical-Bayes / pseudo-bayesianos**. No hay PPL,
> ni random effects formales, ni posterior conjunto muestreado por MCMC. Donde decimos
> "posterior" es la aproximación normal del EB con parámetros de prior estimados de los datos.

---

## Modelo 2 — ranking de reguladores (core, corre local)

**Script:** `scripts/model_hubs.py` · **Corre con:** `DE_stats.suppl_table.csv` (local, sin descargas).
**Grano:** 1 fila = gen perturbado × condición.

### Especificación

1. **Efectos fijos (media condicional).** GLM sobre `n_downstream`:

   ```
   n_downstream ~ C(culture_condition) + ontarget_significant + offtarget_flag
   ```

   Poisson y NB comparten el mismo modelo de media; como solo usamos la media ajustada `μᵢ`
   (no inferencia sobre coeficientes) ajustamos por IRLS estable: Poisson GLM → `α` de NB por
   método de momentos (`Var = μ + α·μ²`) → NB GLM con `α` fijo. Evita los problemas de
   convergencia del NB por MLE completo.

2. **Shrinkage empirical-Bayes del efecto por gen.** Desviación log-rate respecto al baseline:

   ```
   workᵢ = log(yᵢ + 0.5) − log(μᵢ + 0.5)
   ```

   Por gen g:  `d_g = mean(work)`,  `s²_g = σ²_e / n_g`.
   Prior `u_g ~ Normal(0, τ²)` con `τ²` por método de momentos (`Var(d_g) − mean(s²_g)`).
   Posterior aproximado:

   ```
   u_g | datos ~ Normal( shrink·d_g ,  shrink·s²_g ),   shrink = τ²/(τ²+s²_g)
   ```

   Genes con pocas condiciones / poca señal se encogen hacia 0.

### Salidas

`docs/tables/hub_ranking_bayes.csv` (todos los genes) y `docs/tables/top_regulators_for_review.csv`
(top 30, judge-facing). Columnas clave: `regpower_eb_mean/sd` (poder regulatorio log-rate),
`p_top_1pct` (probabilidad EB de exceder el umbral empírico del top-1%, no "P de estar en el top 1%"),
`expected_downstream`. Figura `07_hub_posterior_ranking.png`.

### Lectura del resultado

El ranking robusto surface maquinaria de **cromatina/transcripción** consistente entre condiciones
—complejo SAGA (TADA1/TADA2B/SGF29/SUPT20H/TAF6L), Mediador (MED12/CCNC), KDM1A, SETD2, CTBP1—
por encima de los hubs de señalización TCR crudos que eran específicos de Stim8hr. Es decir: el
shrinkage premia a los reguladores con efecto grande **y** estable.

### Caveats

- `xcond_reproducibility` es una **feature exploratoria** (estabilidad cross-condición). **No**
  sustituye la reproducibilidad cross-donor / cross-guide, que requiere `DE_stats.h5ad`.
- El baseline de efectos fijos se trata como conocido (plug-in) → pseudo-bayesiano, no full-Bayes.
- `single_guide_estimate` y `n_guides` NO están en el CSV; en la tabla de review del core aparecen como
  `NA (requiere DE_stats.h5ad)` — sí están en la auditoría de sensibilidad de abajo.

### Auditoría de sensibilidad guide/donor-aware (opcional)

Cuando existe `de_obs_reproducibility_metadata.csv` (extraído del `.obs` de `DE_stats.h5ad`,
sin `.layers`), `model_hubs.py` corre una auditoría: **repondera** el score EB con reproducibilidad
real (`reweighted_score = regpower_eb_mean · repro_weight`) y reporta qué reguladores sobreviven
(`reproducibility_audit.csv`, fig 19).

- **Es un análisis de sensibilidad, no un posterior nuevo**: NO se reestima el modelo EB.
- **Cobertura parcial**: `guide_correlation_all` ~78% de contrastes, `donor_correlation_hits_mean`
  solo ~19% → en la práctica es más *guide-aware* que *donor-aware*. Donde falta la métrica se usa un
  **peso neutral** (0.75), de modo que **un gen no se penaliza solo por no tener metadata de donante**.
- El ranking **core no depende** de este archivo (`make all` corre sin él).

---

## Modelo 1 — red de efectos con incertidumbre (ESTRICTAMENTE OPCIONAL)

**Scripts:** `scripts/model_edges_spike.py` (validación) y `scripts/model_edges.py` (escalado).
**Regla:** si el spike remoto falla o es lento, el entregable oficial es el Modelo 2 + docs.

### Idea

EB normal-normal exacto sobre `log_fc` / `lfcSE` de los `.layers` del h5ad:

```
yᵢ | θᵢ ~ Normal(θᵢ, seᵢ²)          # observado
θᵢ     ~ Normal(0, τ²)              # prior con shrinkage
θᵢ|yᵢ  ~ Normal(mᵢ, vᵢ),  vᵢ = 1/(1/τ² + 1/seᵢ²),  mᵢ = vᵢ·yᵢ/seᵢ²
```

Salidas por edge: `theta_post_mean/sd`, `p_effect_positive`, `p_abs_effect_gt_1p5x`.
**Regla de decisión** (más interpretable que FDR): `p_abs_effect_gt_1p5x > 0.8 AND ontarget_significant`.

### Estrategia consciente de memoria/cómputo

Disco: **9.8 GB libres < 17 GB** del h5ad → no se descarga. En vez de eso:
- Solo se necesitan las edges de los **reguladores candidatos** (top del Modelo 2), no las ~350M.
- `model_edges_spike.py` **mide** (no asume) el layout/chunking y el coste real de leer una fila
  por slice desde S3 (`fsspec` anónimo + `h5py`). Si es viable, `model_edges.py` baja solo esas
  filas y corre el EB vectorizado (segundos, ~15 MB de RAM).
- `τ²` se estima de una muestra de filas, no de toda la matriz (aproximación documentada).

Ver el veredicto real del spike en `docs/report.md` (sección Modelo 1).

---

## Cómo correr

```bash
make model          # Modelo 2 (core)
make spike          # Modelo 1 spike (opcional, requiere: pip install h5py s3fs fsspec)
```

## Next steps (no incluidos)

- Reforzar `regpower` con reproducibilidad cross-donor/cross-guide real desde `DE_stats.h5ad`.
- Término condición-específico `γ_{p,c,g}` y prior spike-and-slab (`z ~ Bernoulli(π)`) para la red.
- Full-Bayes (NumPyro/PyMC) si el EB deja de ser suficiente.
