# Side analysis — transcriptional fingerprints / perturbation program similarity

## What question it answers

The EB ranking and the audits say **who** is a strong regulator. This analysis answers a different
question: **what transcriptomic program each perturbation induces, and which regulators resemble
each other**. It moves from a *list of genes* to a *map of programs*.

A perturbation's **fingerprint** = its vector of effects over the ~10k measured genes. Two
regulators with similar fingerprints act on the same program. On that matrix we run PCA,
cosine similarity, and clustering. It does not replace the core ranking, adds no heavy models, and
does not download the full h5ad.

## Data and method

The `layers/zscore` matrix from the remote h5ad (dense, contiguous): rows = perturbation×condition,
columns = measured genes. We read **only the panel** (~200 rows) by remote slice, cached in
`data/cache/panel_zscore_*.npy`. (`--matrix log_fc` uses the full 1.4 GB local cache.)

**BALANCED panel** — key to not biasing the map. If you take only the top-EB, you get pure
chromatin/SAGA/Mediator. Instead:

| source | n | criterion |
|--------|---|-----------|
| global            | 75 | top by EB regpower, class `global` |
| context-specific  | 75 | top by EB regpower, class `condition-specific` |
| promoted          | 25 | promoted by the reproducibility audit (rank_shift↑) |
| demoted           | 25 | demoted by the audit (rank_shift↓) |

Each regulator enters in its peak condition. Columns: **top-2000 genes by variance** in the panel.
Columns are **standardized** (z-score per gene) before PCA/similarity, so comparable genes carry
equal weight.

Run: `make fingerprints`  (or `python scripts/analyze_fingerprints.py --n 200 --matrix zscore --top-genes 2000`).

## Results

### 1. The latent space captures *programs*, not *magnitude*

PCA over the fingerprints. If PC1 merely reordered by "effect size", the analysis would add nothing
over the ranking. **That is not the case:**

- `spearman(|PC1|, n_downstream) = 0.25` — weak correlation. PC1 is **not** magnitude.
- PC1 (28.9% var) is the axis of **TCR signaling / proximal activation**; PC2 (7.4%) separates the
  **chromatin coactivators** (SAGA/Mediator, top) from the rest.
- The top loadings are immunologically coherent (effector cytokines, naive/memory,
  Th polarization, interferon).

→ **The ranking captures strength; the latent map captures program identity.** Complementary.

### 2. The space recovers known biological structure (permutation validation)

Do members of a complex have more similar fingerprints to each other than random gene-sets from the
panel? (intra-complex cosine vs. permutation null, N=5000):

| Complex   | n | intra cosine | null | z | p |
|-----------|---|--------------|------|-----|-------|
| **TCR**       | 5 | **0.916** | 0.004 | **11.2** | 0.0002 |
| **SAGA**      | 8 | **0.478** | 0.004 | **9.2**  | 0.0002 |
| **Mediator**  | 7 | 0.196 | 0.007 | 3.2 | 0.013 |

All three are significant. TCR is nearly colinear (0.92): CD3E/ZAP70/LAT/LCK/PLCG1 induce
practically the same fingerprint. SAGA is very cohesive. Mediator is more diffuse but significant.

### 3. Transcriptomic neighbors — the demoable view

Nearest neighbors by fingerprint cosine similarity (`fingerprint_neighbors.csv`, endpoint
`/programs/neighbors/{gene}`):

- **ZAP70** → LAT (0.95), CD3E (0.95), PLCG1 (0.95), LCK (0.84) — pure TCR module.
- **TADA2B** → TAF6L (0.94), SUPT20H (0.92), TADA1 (0.88) — pure SAGA module.
- **MED12** → MED19 (0.78), SUPT20H (0.75), TADA2B (0.75) — Mediator + SAGA (coactivators).
- **SGF29** → CCNC (0.66) — SAGA closest to a Mediator member.

The atlas no longer just ranks regulators: it **organizes perturbations by the fingerprint they
induce** and lets you ask "who does this gene resemble?".

## Outputs

**Figures** (`docs/figures/`): `20_fingerprint_pca_by_condition.png`,
`21_fingerprint_pca_by_regulator_class.png` (main), `22_fingerprint_similarity_heatmap.png`,
`23_fingerprint_neighbor_network.png`.

**Tables** (`docs/tables/`): `fingerprint_panel.csv`, `fingerprint_pca_scores.csv`,
`fingerprint_neighbors.csv`, `fingerprint_similarity_edges.csv`, `fingerprint_clusters.csv`,
`fingerprint_pc_loadings.csv`, `fingerprint_complex_validation.csv`, `fingerprint_summary.json`.

**API** (`/programs/*`): `GET /programs/pca`, `GET /programs/neighbors/{gene}`,
`GET /programs/clusters`. **UI**: the **Programs** tab (PCA scatter + neighbor search).

## Limitations (honest)

- **One fingerprint per regulator (peak condition)** → the map mixes conditions. The per-condition
  figure helps read the context effect, but a separate PCA per condition would be cleaner.
- **zscore = log_fc/lfcSE** already weights by error, but we standardize per column on top; it is a
  reasonable choice, not the only one.
- **Linear PCA/SVD**, not UMAP: more transparent and defensible, but does not capture nonlinearities.
- **Exploratory**: it is *perturbation program similarity*, **not** causal pathway discovery nor
  discovery of new cell programs. The complexes used are known annotation, not discovery.
- **Panel of 200**, not all 33,983 contrasts: enough for the map, not exhaustive.

## Verdict

**Strong result.** Three quantified, defensible claims: (1) the latent structure is program, not
magnitude (|PC1|~n_downstream = 0.25); (2) the space recovers known complexes with permutation
significance (TCR z=11, SAGA z=9, Mediator z=3); (3) the transcriptomic neighbors are biologically
coherent and navigable from the UI. The atlas moves from "regulator ranking" to "map of regulatory
programs".
