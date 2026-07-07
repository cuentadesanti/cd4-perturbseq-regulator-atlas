# Edge analysis (Model 1) — strong result or bonus?

## What was analyzed

The uncertainty-aware effect network `docs/tables/robust_edges.csv`: **2,470** robust edges
(`P(|effect|>1.5×)>0.8`) from **6 regulators** × **2 conditions**, read by slice from the
remote `DE_stats.h5ad`. Summarized per regulator and per downstream gene, with direction and convergence.

## Do they look useful?

**Yes as a proof-of-concept, not as a strong result.** Signals in favor:

- **Coherent convergence**: 579 downstream genes are targeted by ≥2 regulators. All 6
  regulators are co-members of the **SAGA complex** (TADA1/TADA2B/SGF29/SUPT20H) → sharing
  targets is exactly what's expected, and **validates that the method recovers real biological structure**.
- **Interpretable direction**: 83% of the edges are activating, consistent with
  SAGA as a transcriptional coactivator.
- Effect magnitudes are modest and well bounded (median |θ| ≈ 0.92).

## Why it is NOT a strong result (yet)

- **Minimal coverage**: 6 of 7,913 regulators (0.08%), **selected by the EB ranking**
  → a sample biased toward a single complex, not representative of the regulatory landscape.
- **Incomplete conditions**: only Stim48hr, Stim8hr (the demo took the
  peak condition per regulator; **Rest is missing**).
- **Remote latency**: scaling to ~150 regulators is ~11 min of h5ad reads (4.5 s/row measured).
- **Probability semantics**: `p_abs_effect_gt_1p5x` is **P(effect magnitude > 1.5×)**,
  NOT P(a causal edge exists). There is no network-level FDR control or sparsity (spike-and-slab).

## Verdict

Keep it as a **bonus / proof-of-concept**: it shows the pipeline produces a biologically sensible
uncertainty-aware effect network without downloading 1.8 TB, but the coverage and semantics do not
yet support network-scale claims. To promote it to a strong result: run the ~150 top regulators over
all 3 conditions and add P(edge) with a sparse prior.
