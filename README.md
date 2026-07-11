# Robust regulators of CD4+ T-cell programs from genome-scale Perturb-seq

*A quality-aware, donor-audited framework for identifying robust regulators and transcriptional response programs in primary human CD4+ T cells*

## Overview

Genome-scale Perturb-seq screens contain thousands of statistically significant effects, but significance alone does not identify reliable regulators.

Raw differential-expression counts are strongly influenced by:

* perturbation quality
* statistical power
* a small number of extreme hubs
* condition-specific responses
* instability across guides and donors

This project develops a reproducible framework for distinguishing large, stable and replicable regulatory effects from these sources of variation.

The analysis produces three primary outputs:

1. A quality-aware ranking of robust regulators
2. Audits of guide-, condition- and donor-level reproducibility
3. Transcriptional response programs induced by selected perturbations

The core analysis runs from approximately 15 MB of supplementary tables, despite the complete dataset occupying roughly 1.8 TB.

## Main contribution

The central contribution is not a new differential-expression test.

It is a framework for deciding which perturbation effects remain credible after accounting for knockdown quality, effect uncertainty, context specificity and biological replication.

The resulting regulator atlas separates:

* robust global regulators from condition-specific hubs
* reproducible effects from guide- or donor-fragile signals
* perturbation magnitude from response-program identity
* exploratory biological hypotheses from evidence supported by the screen

## Dataset

**Primary Human CD4+ T Cell Perturb-seq**

* Genome-scale CRISPR interference screen
* Primary human CD4+ T cells
* 4 donors
* 3 conditions: Rest, Stim8hr and Stim48hr
* 33,983 differential-expression contrasts
* 26,504 guides
* 12 donor-condition samples

Data source: CZI Virtual Cells Platform, Marson Lab, 2025

The local core uses three supplementary CSV files. Larger `.h5ad` and `.h5mu` objects are accessed selectively for analyses that require gene-level perturbation profiles or donor-specific effects.

## Analysis workflow

### 1. Perturbation quality control

A significant on-target knockdown is required before interpreting downstream effects.

Only 62% of contrasts pass this gate, but they account for approximately 85% of all detected trans-effects. This makes perturbation quality control a necessary first step rather than an optional sensitivity analysis.

### 2. Quality-aware regulator ranking

An empirical-Bayes model combines effect magnitude with uncertainty and cross-condition behavior.

This prevents the ranking from being dominated by:

* noisy large estimates
* raw DEG counts
* weakly validated knockdowns
* effects observed in only one favorable condition

The output is best interpreted as a set of high-confidence regulators, not as an exact ordering of genes.

### 3. Reproducibility audits

The ranking is pressure-tested using:

* bootstrap rank stability
* alternative baseline rankings
* guide-level agreement
* cross-condition consistency
* donor-level replication
* within-condition null comparisons

Donor robustness is reported as an audit column rather than being used to silently reorder the ranking. Large but donor-fragile effects therefore remain visible.

### 4. Transcriptional response programs

For selected regulators, the full downstream expression fingerprint is used to characterize what each perturbation does to the cell.

Fingerprint similarity recovers coherent response neighborhoods associated with:

* TCR signaling
* SAGA and chromatin regulation
* Mediator and transcriptional control

Program labels indicate similarity of perturbation response. They do not imply physical protein-complex membership.

### 5. Empirical regulatory operator

The complete regulator-by-response-gene log-fold-change matrix is also analyzed as an empirical operator.

This extension investigates:

* low-dimensional gene-response programs
* condition gating
* regulator communities
* out-of-panel reconstruction
* inductive prediction for unseen regulators

The leave-regulator-out experiments provide an important negative result: response structure is strongly compressible within the observed matrix, but current regulator-side features do not reliably predict the effects of entirely unseen perturbations.

This distinguishes transductive structure recovery from genuine inductive generalization.

## Principal findings

### Robust regulators differ from raw hubs

Ranking genes by the number of downstream DEGs overemphasizes perturbations that are condition-specific, noisy or insufficiently validated.

The quality-aware ranking instead prioritizes stable chromatin and transcriptional regulators, including members of:

* the SAGA complex
* the Mediator complex
* chromatin-modifying machinery

TCR-signaling genes remain biologically meaningful, but many are correctly classified as stimulation-dependent rather than universal regulators.

### Donor replication changes interpretation

The strongest regulators are donor-robust as a class, with 29 of the top 30 passing the donor-concordance criterion.

The audit also identifies prominent exceptions. For example, a regulator may have thousands of downstream effects yet fail cross-donor replication. Such genes are retained in the tables but explicitly flagged.

### Perturbation fingerprints recover recognizable programs

Response similarity recovers coherent SAGA/chromatin, Mediator and TCR neighborhoods.

The classifier is deliberately conservative. Most regulators remain unassigned rather than being forced into a known program.

Candidate neighbors are treated as hypotheses based on response similarity. For example, CHD7 is placed near the SAGA/chromatin response program, but this is not presented as evidence of physical complex membership.

### Strong perturbations share a stimulation-linked interferon response

SAGA-family knockdown produces a consistent de-repressive interferon response.

However, effect-size-matched controls show that much of the enrichment magnitude is shared by strong perturbations under stimulation. The more specific SAGA-associated observation is the consistency of the response direction, not uniquely high enrichment.

### Inductive prediction remains unresolved

Low-rank and nonlinear models can reconstruct held-out entries or regulators when information from the same observed response system remains available.

They do not yet generalize convincingly to entirely unseen regulators using the available regulator annotations.

This negative result is retained because it defines the current boundary of what can be inferred from the dataset.

## Repository structure

```
.
├── api/                    # Read-only FastAPI service and atlas UI
├── docs/
│   ├── figures/            # Generated figures
│   ├── tables/             # Generated analysis tables
│   ├── report.md           # Consolidated analysis report
│   ├── MODELING.md         # Ranking methodology
│   ├── DONOR_ROBUSTNESS.md # Donor-level audit
│   ├── FINGERPRINT_ANALYSIS.md
│   └── OPERATOR_ANALYSIS.md
├── scripts/                # Reproducible analysis pipeline
├── Makefile                # Main workflow entry points
└── requirements.txt
```

## Reproducing the core analysis

```bash
# Create an environment
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Download the supplementary tables
scripts/download.sh tables

# Run the verified core pipeline
make all
```

The core pipeline includes:

```bash
make eda
make model
make audit
make report
```

Optional analyses that access larger remote objects require additional dependencies:

```bash
pip install h5py s3fs fsspec scikit-learn
```

Selected optional targets include:

```bash
make fingerprints
make edges
make operator
```

See the `Makefile` and the corresponding documents under `docs/` for the complete workflow.

## Regulator Atlas

The repository includes a read-only API and browser interface for exploring the generated artifacts.

The application does not train models or download the full dataset. It serves versioned outputs created by the analysis pipeline.

```bash
make all
make api
```

Then open:

```
http://localhost:8000/
```

Swagger documentation is available at:

```
http://localhost:8000/docs
```

The atlas supports:

* regulator search and ranking
* per-condition effect profiles
* knockdown and reproducibility audits
* donor-robustness flags
* transcriptional-program assignments
* fingerprint neighbors and response markers
* downstream effect-network exploration
* regulator-class comparisons

See `api/README.md` for details.

## Key outputs

| Output | Description |
|---|---|
| `docs/report.md` | Consolidated analysis report |
| `docs/tables/top_regulators_for_review.csv` | Review-facing table of top regulators |
| `docs/tables/hub_ranking_bayes.csv` | Full empirical-Bayes regulator ranking |
| `docs/tables/fingerprint_findings.csv` | Program assignments, neighbors and markers |
| `docs/tables/program_label_evidence.csv` | Evidence supporting each program label |
| `docs/tables/class_convergent_targets.csv` | Convergent response genes by regulator class |
| `docs/tables/class_isg_enrichment.csv` | Interferon enrichment and specificity controls |
| `docs/tables/robust_edges.csv` | Uncertainty-aware effect-network edges |
| `docs/figures/` | Generated figures and diagnostic plots |
| `docs/data-model.html` | Interactive data-model explorer |

## Interpretation guide

Several distinctions are important when reading the results.

**Robust does not mean universally causal**

The ranking quantifies evidence for a large and stable perturbation response within this experiment. It is not a posterior probability that a gene is a causal regulator in every biological context.

**Program similarity is not complex membership**

A gene assigned to a response program produces a transcriptomic fingerprint similar to that program's reference regulators. This does not establish a physical interaction or direct pathway membership.

**Donor robustness is an audit, not a hidden filter**

Genes that fail donor replication are retained and flagged. This avoids presenting a cleaner ranking by silently removing inconvenient observations.

**Negative results are part of the deliverable**

Failed specificity tests and unsuccessful inductive-prediction experiments are reported alongside positive findings. They constrain the biological and predictive claims that the dataset can support.

## Limitations

* The ranking model is empirical Bayes or pseudo-Bayesian, not a full hierarchical count model
* Differential-expression summaries discard some cell-level information
* Guide-level reproducibility is unavailable for part of the original table-only core
* Program assignments depend on the selected perturbation panel and reference complexes
* Response similarity cannot establish direct molecular interactions
* Remote slicing of large expression objects is latency-bound
* Inductive prediction for unseen regulators remains weak with the currently available regulator-side features
* Disease and GWAS overlaps are hypothesis-generating and are not tested mechanistically here

## Documentation

Start with:

* `docs/report.md` for the complete analysis
* `docs/MODELING.md` for regulator ranking
* `docs/DONOR_ROBUSTNESS.md` for donor replication
* `docs/FINGERPRINT_ANALYSIS.md` for response programs
* `docs/OPERATOR_ANALYSIS.md` for the empirical operator
* `docs/REVIEW_RESPONSE.md` for robustness checks and reviewer-style critiques

## Scope

This repository provides a reproducible, audit-first reanalysis of a genome-scale primary-cell Perturb-seq screen.

Its strongest claims concern:

1. the identification of quality-aware and donor-audited regulators
2. the separation of global from context-specific effects
3. the organization of perturbations into reproducible response programs

The operator, network, disease-overlap and inductive-prediction analyses extend this core result, but are not required to support it.

---

Data: CZI Virtual Cells Platform · Marson Lab 2025
Preprint: 10.64898/2025.12.23.696273
Original analysis code: emdann/GWT_perturbseq_analysis_2025
