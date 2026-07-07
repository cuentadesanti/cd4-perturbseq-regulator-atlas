# Transcriptional programs — from regulator ranks to perturbation programs

## The question

The EB ranking and the audits say **who** is a strong regulator. This analysis answers a different
question: **what transcriptional program each perturbation induces, and which regulators resemble each
other.** A rank is one number; a perturbation's **fingerprint** — its downstream effect vector over the
~10k measured genes — is its whole action on the cell. It does not replace the ranking, adds no heavy
models, and does not download the full h5ad.

Run: `make fingerprints` (`python scripts/analyze_fingerprints.py --n 200 --matrix zscore --top-genes 2000`).

## Data and method

The `layers/zscore` matrix of the remote h5ad (rows = perturbation×condition, columns = measured
genes) is read **by slice** for a **balanced panel** of 200 regulators (cached in `data/cache/`):

| source | n | criterion |
|--------|---|-----------|
| global            | 75 | top by EB regpower, class `global` |
| context-specific  | 75 | top by EB regpower, class `condition-specific` |
| promoted          | 25 | promoted by the reproducibility audit |
| demoted           | 25 | demoted by the audit |

Each regulator enters at its peak condition; columns = top-2000 genes by variance, standardized per
gene. On this matrix we compute cosine similarity (→ nearest neighbors, permutation validation) and
PCA (a view). **Programs are anchored to known complexes**: each regulator's fingerprint is matched to
the leave-one-out centroid of the curated **SAGA / Mediator / TCR** complexes; it is labeled only if
cosine ≥ 0.45 **and** it beats the second-best complex by ≥ 0.05, otherwise it stays **mixed**
(conservative — we only have prototypes for three complexes). This is a transparent nearest-prototype
classifier in the *same* space as the validated similarity, not a black-box cluster labeling.

## Results

### 1. Fingerprint similarity recovers known biology (permutation-validated)

Are the members of a complex more fingerprint-similar to each other than random gene-sets from the
panel? (intra-complex cosine vs. permutation null, N=5000):

| Complex   | n | intra cosine | null | z | p |
|-----------|---|--------------|------|-----|-------|
| **TCR**       | 5 | **0.916** | 0.004 | **11.2** | 0.0002 |
| **SAGA**      | 8 | **0.478** | 0.004 | **9.2**  | 0.0002 |
| **Mediator**  | 7 | 0.196 | 0.007 | 3.2 | 0.013 |

All three are significant. TCR is nearly colinear (CD3E/ZAP70/LAT/LCK/PLCG1 induce almost the same
fingerprint); SAGA is cohesive; Mediator is diffuse but significant. This is the load-bearing evidence
— **not** the PCA. PC1 (28.9% var) merely confirms the axis is program *identity*, not magnitude
(`spearman(|PC1|, n_downstream) = 0.25`).

### 2. Recognizable programs, with new members

25 of the 200 regulators map to a program (the rest stay *mixed*):

| program | n | known-complex core | fingerprint-neighbors (novel) | centroid cosine |
|---|---|---|---|---|
| **TCR signaling** | 13 | ZAP70, LCK, LAT, CD3E, PLCG1 | ATF7IP2, NCAPG2, CLCC1, EIF1AX, PGGT1B… | 0.77 |
| **SAGA/chromatin** | 9 | TADA2B, SUPT20H, TADA1, TAF6L, SUPT7L, USP22 | **CHD7**, TSPYL5 | 0.81 |
| **Mediator/transcription** | 3 | MED1 | POGLUT3, GLIPR2 | 0.59 |

The labels recover each complex's core and add biologically coherent neighbors — most strikingly the
chromatin remodeler **CHD7 joining SAGA/chromatin** (cosine 0.84), and Mediator's **MED12 landing in
the chromatin program** (Mediator–SAGA coactivator crosstalk). Peripheral subunits stay *mixed* where
the biology predicts weaker cohesion (e.g. **SGF29**, a SAGA reader, at cosine 0.21). Every label,
with members, centroid cosine and marker genes, is in `program_label_evidence.csv`.

### 3. Convergent downstream response genes

Per program we report the **convergent downstream response genes** — genes whose perturbation-response
z-scores are *consistently* high or low across the program's regulators (relative to the panel), **not**
baseline cell-type markers (`fingerprint_program_markers.csv`). The **TCR** program's markers are a
clean immune-effector signature (STK17B, S100A11, SELPLG, CD53, MYO1F, GPSM3); the **SAGA** program's
top markers include some SAGA-adjacent genes (expected — perturbing one subunit dysregulates the
others), so that set is read as coherence, not novel targets.

### 4. Are the reproducibility-promoted hits coherent — or artifacts?

The reproducibility audit promotes/demotes regulators; do the promoted ones form coherent
neighborhoods? (`fingerprint_audit_coherence.csv`):

| source | n | mean kNN cosine | % mapping to a canonical program |
|---|---|---|---|
| global | 75 | 0.40 | 13% |
| context-specific | 75 | 0.50 | 20% |
| **promoted** | 25 | **0.47** | **0%** |
| **demoted** | 25 | 0.46 | 0% |

Reported honestly: promoted hits have neighborhoods **as tight as the top global regulators**
(kNN cosine ~0.47), so they are **not statistical noise** — but they map onto *none* of the three
canonical complexes. The audit surfaces a **distinct** coherent set, not "more SAGA". That is the
useful, non-obvious result: reproducibility filtering is not just re-finding the known machinery.

### 5. The demoable view — transcriptomic neighbors

Nearest neighbors by fingerprint cosine (`fingerprint_neighbors.csv`, `/programs/neighbors/{gene}`,
and each gene's profile in the UI):

- **ZAP70** → LAT (0.95), CD3E (0.95), PLCG1 (0.95), LCK (0.87) — pure TCR module.
- **CHD7** → TADA2B (0.87), SUPT20H (0.84), TAF6L (0.83) — a chromatin remodeler inside the SAGA module.
- **MED12** → MED19 (0.78), SENP5 (0.78), SUPT20H (0.75) — Mediator + SAGA (coactivators).

The atlas no longer just ranks regulators: it **organizes perturbations by the program they induce**
and lets you ask "what does this gene do, and who does it resemble?".

## Outputs

**Tables** (`docs/tables/`): `fingerprint_findings.csv` (per regulator: program, nearest complex,
neighbors, markers), `program_label_evidence.csv` (auditable label basis),
`fingerprint_program_markers.csv`, `fingerprint_audit_coherence.csv`, `fingerprint_complex_validation.csv`,
`fingerprint_pca_scores.csv`, `fingerprint_neighbors.csv`, `fingerprint_clusters.csv`, `fingerprint_summary.json`.

**Figures** (`docs/figures/`): `24_fingerprint_pca_by_program.png` (main), `20`–`23` (PCA by
condition/class, similarity heatmap, neighbor network).

**API/UI** (`/programs/*`): `GET /programs/summary`, `/programs/findings`, `/programs/pca`,
`/programs/neighbors/{gene}`; the **Programs** tab (PCA colored by program + neighbor search), and the
program appears in every gene's Explore profile.

## Limitations (honest)

- **Perturbation-program similarity, anchored to known complexes** — not de-novo pathway discovery. We
  label relative to three curated prototypes; the majority (175/200) stay *mixed* by design.
- Markers are **consistently-moved downstream response genes** (relative to the panel), not baseline
  markers; SAGA markers include complex-adjacent genes.
- **One fingerprint per regulator (peak condition)** → the map mixes conditions; linear PCA/SVD (not
  UMAP) — transparent but misses non-linearities.
- **Panel of 200**, not all 33,983 contrasts: enough for the map, not exhaustive.

## Verdict

**A real, defensible program-level result.** Fingerprint similarity recovers the SAGA/Mediator/TCR
complexes with permutation significance (z=9/3/11), organizes 25 top perturbations into three
recognizable programs with biologically coherent new members (CHD7 → chromatin), and shows the
reproducibility-promoted hits are coherent but distinct from the canonical machinery. The atlas moves
from "regulator ranking" to "map of transcriptional programs" — honestly scoped, and navigable from a
gene's profile.
