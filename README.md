# Robust regulators of CD4+ T cell programs — Perturb-seq

> Hackathon submission · genome-scale CRISPRi Perturb-seq in primary human CD4+ T cells (Marson Lab, 2025)

![pipeline](docs/figures/00_pipeline_overview.png)

## The question

Which genes are **robust regulators** of CD4+ T cell programs? The challenge: the signal from a
genome-scale screen is dominated by noise and by a handful of hubs. We want to **separate signal
from noise with uncertainty** and prioritize regulators by a **large and reproducible** effect —
not by raw counts or `adj_p_value < 0.1`.

## What this repo does

This project turns a massive CD4+ T cell Perturb-seq screen into **three explorable objects: robust
regulators, reproducibility audits, and transcriptional programs.** The CSV-only core identifies
quality-aware regulators; optional h5ad slicing maps selected perturbations into transcriptomic
fingerprint neighborhoods. It is **memory- and compute-aware** (the full dataset is 1.8 TB; the
working disk has ~10 GB free) — the entire core runs **from the supplementary CSV tables alone
(~15 MB)**:

1. **Selective download pipeline** from the public S3 bucket (`scripts/download.sh`).
2. **Data model** of the dataset (`docs/DATA_MODEL.md` + [interactive artifact](docs/data-model.html)).
3. **80/20 EDA** (`scripts/eda.py`) — effect-size distribution, knockdown quality, hubs, reproducibility.
4. **Model 2 · empirical Bayes** (`scripts/model_hubs.py`) — regulator ranking with uncertainty.
5. **Ranking audits** (`scripts/audit_ranking.py`) — baselines, bootstrap stability, global vs.
   context-specific, and a guide/donor-aware reproducibility audit.
6. **Transcriptional programs** (`scripts/analyze_fingerprints.py`) — organizes top perturbations by
   the fingerprint they induce, reading `.h5ad` layers *by slice* from S3 **without downloading it**.
7. **Model 1 · uncertainty-aware effect network** (optional, `scripts/model_edges.py`, bonus).

## Dataset

**Primary Human CD4+ T Cell Perturb-seq** · [CZI Virtual Cells Platform](https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq)
· public bucket `s3://genome-scale-tcell-perturb-seq/marson2025_data/` · 4 donors × 3 conditions
(Rest / Stim8hr / Stim48hr). The core uses 3 CSVs (33,983 DE contrasts, 26,504 guides, 12 samples).

## Key findings

- **Heavy-tailed effects**: the median perturbation yields 2 DEGs, 15% have no effect, but 1.5% are
  hubs (>1000 DEGs). → summarize with percentiles and rankings, not the mean.
- **Knockdown gates the signal**: contrasts with a significant on-target knockdown (62%) concentrate
  **85%** of all trans-effects. Filtering by `ontarget_significant` is the mandatory first step.
- **Robust ranking ≠ raw hubs**: the EB model surfaces **chromatin/transcription machinery**
  (SAGA complex: TADA1/TADA2B/SGF29/SUPT20H · Mediator: MED12/CCNC · KDM1A, SETD2) — regulators with a
  large **and** stable effect across conditions, ranked above the Stim8hr-specific TCR-signaling hubs.
- **Fingerprint similarity organizes perturbations into recognizable programs**: matching each
  regulator's downstream fingerprint to the curated complexes recovers **TCR signaling, SAGA/chromatin,
  and Mediator/transcription** (permutation-validated at z=11, z=9, z=3) and adds fingerprint-neighbors
  — e.g. the chromatin remodeler **CHD7 joins the SAGA/chromatin program** (cosine 0.84). The latent
  axis is program *identity*, not effect magnitude (|PC1| vs. n_downstream = 0.25).
- **Uncertainty-aware effect network (bonus)**: ~2,470 robust edges (`P(|effect|>1.5×)>0.8`, i.e. the
  probability that the effect *magnitude* exceeds 1.5×, not that a causal edge exists) for the top
  regulators, extracted from the remote h5ad without downloading it.

Detail: [`docs/report.md`](docs/report.md) · [`docs/EDA.md`](docs/EDA.md) · [`docs/MODELING.md`](docs/MODELING.md).

## Validating the ranking

Before trusting the ranking, we audit it (`scripts/audit_ranking.py`, no new dependencies).

**Naive hubs vs. quality-aware regulators.** Ranking by raw `n_downstream` rewards hubs that don't
survive the controls: of the top 30 raw hubs, **2 fall out at the knockdown gate** (no validated
knockdown) and **~15 are demoted** by EB shrinkage for being condition-specific. The EB ranking
surfaces regulators with a large **and** stable effect. Stability was measured by bootstrap (B=200):
it is *moderate* — read it as a set of robust regulators, not an exact ordering (`stability_frequency`
per gene).

**Global vs. context-specific regulators.** Splitting by
`condition_specificity = max/sum of n_downstream`: the **global** ones (SGF29, TADA2B, SUPT20H…) are
chromatin/transcription machinery active in all conditions; the **context-specific** ones (ZAP70,
LCK…) are TCR signaling, active only under stimulation. Both classes are real biology; the
distinction avoids confusing a universal regulator with a context-dependent one. See
`docs/tables/top_global_regulators.csv` and `top_condition_specific_regulators.csv`.

## Transcriptional programs

A rank is one number; a **fingerprint** is a regulator's whole effect on the cell. On a balanced
panel of 200 top perturbations (global · context-specific · reproducibility-promoted · demoted), each
regulator's downstream fingerprint (zscore vector) is compared to the curated **SAGA / Mediator / TCR**
complexes; matches above threshold get a program label, the rest stay *mixed* (conservative — we only
have prototypes for three complexes). `scripts/analyze_fingerprints.py` · `make fingerprints`.

- **The space recovers known biology.** By permutation test the three complexes are each significantly
  cohesive: **TCR z=11, SAGA z=9, Mediator z=3** (N=5000). The latent PC1 is program *identity*, not
  effect magnitude (|PC1| vs. n_downstream Spearman = 0.25).
- **Recognizable programs + new members.** 25 of 200 regulators map to a program: **TCR signaling (13),
  SAGA/chromatin (9), Mediator/transcription (3)**. The labels recover each complex's core and add
  fingerprint-neighbors — e.g. the chromatin remodeler **CHD7** joins SAGA/chromatin (cosine 0.84), and
  Mediator's **MED12** sits in the chromatin program (Mediator–SAGA crosstalk). Every label is auditable
  in `docs/tables/program_label_evidence.csv` with its members, centroid cosine, and marker genes.
- **The reproducibility-promoted hits are coherent but distinct.** Promoted/demoted regulators have
  transcriptomic neighborhoods as tight as the top global regulators (kNN cosine ~0.47 vs. 0.40), so
  they are not statistical noise — yet **none map onto the canonical complexes**: the audit surfaces a
  *separate* coherent set, not "more SAGA".

Tables: `fingerprint_findings.csv` (per-regulator program, neighbors, markers) ·
`fingerprint_program_markers.csv` · `program_label_evidence.csv` · `fingerprint_audit_coherence.csv`.
Figures 20–24. Detail in [`docs/FINGERPRINT_ANALYSIS.md`](docs/FINGERPRINT_ANALYSIS.md).

> *Honest scope:* this is **perturbation-program similarity** anchored to known complexes, not de-novo
> pathway discovery. The convergent "response genes" are genes consistently moved by a program's
> regulators (relative to the panel), not baseline cell-type markers; PCA is a view, not the proof.

## How to reproduce

```bash
# 1. environment
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 2. data (tables only, ~15 MB)
scripts/download.sh tables

# 3. full pipeline (EDA + model + report), verified
make all
```

Targets: `make eda` · `make model` · `make audit` · `make report` · `make all` · `make clean`.
Optional (remote, requires `pip install h5py s3fs fsspec` + `scikit-learn`): `make fingerprints` ·
`make spike` · `make edges`.

## Outputs

| File | What it is |
|---|---|
| `docs/report.md` | consolidated judge-facing report |
| `docs/tables/top_regulators_for_review.csv` | top 30 regulators (judge-facing) |
| `docs/tables/hub_ranking_bayes.csv` | full EB ranking (all genes) |
| `docs/tables/fingerprint_findings.csv` | per-regulator transcriptional program + neighbors + markers |
| `docs/tables/program_label_evidence.csv` | auditable basis for each program label |
| `docs/tables/robust_edges.csv` | uncertainty-aware effect network (bonus, Model 1) |
| `docs/figures/*.png` | EDA · ranking · programs (20–24) · overview |
| `docs/data-model.html` | interactive explorer (data model + EDA + study) |

## Product: Regulator Atlas (API + UI)

*The pipeline produces the artifacts; a read-only API makes them explorable.* The API does **not**
run models, does **not** download the h5ad, and does **not** touch S3 — it only serves the versioned
CSVs from `make all`.

```bash
pip install -r requirements.txt      # includes fastapi + uvicorn
make all                             # generates the tables the API serves
make api                             # uvicorn :8000 (serves API + UI on the same origin)
open http://localhost:8000/          # the UI · Swagger at http://localhost:8000/docs
```

The UI is served by the API itself (same origin), so there's **no file to open** and no CORS to deal
with. If port 8000 is taken, use another (`--port 8010`) and open `http://localhost:8010/`.

Key endpoints: `/summary`, `/regulators?q=&regulator_class=&sort_by=`, **`/regulators/{gene}`**
(full profile: class, per-condition profile, audits, **transcriptional program + neighbors + markers**,
top edges, interpretation), `/audit/reproducibility`, **`/programs/summary`**, `/programs/findings`,
`/edges/downstream`. The UI has **5 screens**: Overview, Explore, Audit, Effect network, and Programs.
Detail in [`api/README.md`](api/README.md).

## Limitations

- **Honest naming**: the models are **empirical-Bayes / pseudo-Bayesian**, not a full hierarchical NB
  or MCMC (no PPL, no formal random effects).
- `xcond_reproducibility` is an **exploratory feature** (cross-condition stability). It is **audited**
  with a **guide/donor-aware sensitivity analysis** (`make repro-meta` →
  `scripts/extract_de_obs_metadata.py`, which extracts only the `.obs`, no `.layers`) that
  **reweights** the EB score with **real** cross-guide (`guide_correlation_all`) and cross-donor
  (`donor_correlation_hits_mean`) reproducibility — it is a **sensitivity analysis, not a new
  posterior**. Optional, and it does not replace the core; see *Sensitivity audit* in
  `docs/report.md`. Honest coverage: cross-guide ~78% but cross-donor only ~19% of contrasts → a
  **neutral weight** where it's missing (a gene is not penalized for lacking donor metadata).
- **Model 1 is optional**: remote per-slice access is viable (~4.5 s/row, measured) but latency-bound;
  the official deliverable stands on its own with the local core.

## Submission summary

A reproducible product that turns a 1.8 TB CD4 Perturb-seq screen into **three explorable objects —
robust regulators, reproducibility audits, and transcriptional programs** — runnable on a laptop with
~10 GB of disk using only 15 MB of data for the core. The **CSV-only** core ranks regulators with
uncertainty (empirical Bayes) and audits them (bootstrap stability + a guide/donor-aware
reproducibility audit). On top of that, **fingerprint similarity organizes the top perturbations into
recognizable programs** — recovering the SAGA, Mediator and TCR complexes (permutation z=9/3/11) and
adding new members like CHD7 → chromatin — plus a bonus **uncertainty-aware effect network**, both
streamed from the 17 GB h5ad without downloading it. An explorable **Regulator Atlas** (read-only
FastAPI + UI) ties it together: search a gene and see its rank, audit survival, transcriptional
program, transcriptomic neighbors, and defining response genes in one view. `make all` reproduces the
core; `make fingerprints` builds the programs; `make api` launches the atlas. See `docs/report.md`.

---
Data: CZI Virtual Cells Platform · Marson Lab 2025 · biorxiv preprint `10.64898/2025.12.23.696273`.
Original analysis code: [emdann/GWT_perturbseq_analysis_2025](https://github.com/emdann/GWT_perturbseq_analysis_2025).
