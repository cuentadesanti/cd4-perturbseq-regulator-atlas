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
   context-specific, and a guide-level reproducibility audit (with partial, 19%, donor coverage).
6. **Transcriptional programs** (`scripts/analyze_fingerprints.py`) — organizes top perturbations by
   the fingerprint they induce, reading `.h5ad` layers *by slice* from S3 **without downloading it**.
7. **Model 1 · uncertainty-aware effect network** (optional, `scripts/model_edges.py`, bonus).
8. **Empirical regulatory operator** (`make operator`) — treats the log-FC matrix as one operator and
   recovers what the fingerprint PCA discarded: gene programs, condition-gating (CP with bootstrap CIs),
   and out-of-panel prediction, all built in precision-decoupled z-score space. See
   `docs/OPERATOR_ANALYSIS.md`.

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
  and Mediator/transcription** (permutation-validated at z=11, z=9, z=3) and surfaces candidate
  neighbors — e.g. the chromatin remodeler **CHD7 is assigned to the SAGA/chromatin program** by
  fingerprint similarity (cosine 0.84; a related perturbation response, not a complex-membership claim).
  The latent axis is program *identity*, not effect magnitude (|PC1| vs. n_downstream = 0.25).
- **Uncertainty-aware effect network (bonus)**: ~2,470 robust edges (`P(|effect|>1.5×)>0.8`, i.e. the
  probability that the effect *magnitude* exceeds 1.5×, not that a causal edge exists) for the top
  regulators, extracted from the remote h5ad without downloading it.

Detail: [`docs/report.md`](docs/report.md) · [`docs/EDA.md`](docs/EDA.md) · [`docs/MODELING.md`](docs/MODELING.md).
Rigor audit / response to reviewer: [`docs/REVIEW_RESPONSE.md`](docs/REVIEW_RESPONSE.md)
(effect-size metric audit, within-condition null, external concordance).

## Validating the ranking

Before trusting the ranking, we audit it (`scripts/audit_ranking.py`, no new dependencies). It is
further pressure-tested against a senior-researcher critique in [`docs/REVIEW_RESPONSE.md`](docs/REVIEW_RESPONSE.md)
(effect-size metric vs. DE-count, within-condition null, external concordance) and against the true
replication unit in [`docs/DONOR_ROBUSTNESS.md`](docs/DONOR_ROBUSTNESS.md).

**Donor robustness (the replication unit that matters).** From the dedicated per-donor DE object
(`by_donors.h5mu`, 100% cross-donor coverage on KD-gated contrasts) we add a `donor_robust` flag
(worst of 6 pairwise donor correlations ≥ 0.5) as a **column, not a re-sort** — a large-effect hub
that fails replication stays visible at its rank (e.g. **SMG1, rank 24, 2,683 DEGs, fails donor
concordance**). The top regulators are donor-robust as a class (29/30). The same check over the
fingerprint programs flags **4 assigned-neighbor false positives** (TCR: ATF7IP2/NCAPG2/EIF1AX;
Mediator: GLIPR2) while the SAGA/chromatin neighbor **CHD7 survives** — so the programs are now
internally consistent with the donor-audited ranking.

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

A rank is one number; a **fingerprint** is what the perturbation actually does to the cell. On a
balanced panel of 200 top perturbations (global · context-specific · reproducibility-promoted ·
demoted), each regulator's downstream fingerprint (zscore vector) is compared to the curated **SAGA /
Mediator / TCR** complexes. These are **candidate program assignments by fingerprint similarity — not
claims of physical complex membership.** `scripts/analyze_fingerprints.py` · `make fingerprints`.

- **Fingerprint similarity recovers the known complexes.** By permutation test the three complexes are
  each significantly cohesive: **TCR z=11, SAGA z=9, Mediator z=3** (N=5000). The latent PC1 is program
  *identity*, not effect magnitude (|PC1| vs. n_downstream Spearman = 0.25).
- **The classifier is conservative — only 25 of 200 are assigned; the rest stay *mixed*, by design.**
  The assigned set: **TCR signaling (13), SAGA/chromatin (9), Mediator/transcription (3)** — of which
  **4 are flagged donor-fragile** by the per-donor check (TCR 3: ATF7IP2/NCAPG2/EIF1AX; Mediator 1:
  GLIPR2; SAGA 0), so SAGA/chromatin is the most donor-robust program. Each program
  recovers its curated core and adds **newly assigned neighbors** (non-curated genes placed in the same
  fingerprint neighborhood) — e.g. the chromatin remodeler **CHD7 is assigned to the SAGA/chromatin
  program** (cosine 0.84; a related response, not complex membership), and Mediator's **MED12** lands in
  the chromatin neighborhood (Mediator–SAGA-like response). Every assignment is auditable in
  `docs/tables/program_label_evidence.csv`.
- **The reproducibility-promoted hits are coherent but distinct.** Promoted/demoted regulators have
  transcriptomic neighborhoods as tight as the top global regulators (kNN cosine ~0.47 vs. 0.40), so
  they are not statistical noise — yet **none map onto the canonical complexes**: the audit surfaces a
  *distinct high-confidence set* rather than simply rediscovering known complexes.

Tables: `fingerprint_findings.csv` (per-regulator program, neighbors, markers) ·
`fingerprint_program_markers.csv` · `program_label_evidence.csv` · `fingerprint_audit_coherence.csv`.
Figures 20–24. Detail in [`docs/FINGERPRINT_ANALYSIS.md`](docs/FINGERPRINT_ANALYSIS.md).

> *Honest scope:* fingerprint-based, program-level re-analysis anchored to known complexes — candidate
> assignments and hypotheses, **not** de-novo pathway discovery or novel complex membership. The
> convergent "response genes" are genes consistently moved by a program's regulators (relative to the
> panel), not baseline cell-type markers; PCA is a view, not the proof.

## Convergent programs by regulator class

A rank tells you *who* is strong; the fingerprint tells you *what a perturbation does*. But is
"chromatin machinery recovers as top hubs" just the expected result of perturbing coactivators? To
test that, a **balanced 30-regulator panel** — chosen by *class* (SAGA/chromatin, Mediator, TCR,
other-robust, reproducibility-promoted, demoted control), not by rank — asks whether different classes
converge on *different* downstream programs. Fully offline (`make class-programs`, uses the cached
fingerprint panel; `scripts/analyze_class_programs.py`).

- **Classes converge on distinct programs.** The pairwise Jaccard of per-class convergent-target sets
  has a **median off-diagonal of ~0.05** — the classes barely share targets, so the "programs" are
  real, not an artifact of every strong perturbation moving the same genes.
- **A convergent interferon module** (SAGA-family knockdown *de-represses* interferon). Genes hit by
  ≥4 of the 6 robust SAGA-family regulators form a **163-gene module** enriched for interferon-
  stimulated genes, **all de-repressive** (knockdown *raises* ISGs).
  > **Specificity control (important requalification).** A dedicated control
  > (`make specificity-control`, fig 29) shows the ISG-enrichment *magnitude* is **largely a general
  > strong-perturbation-under-stimulation effect**: effect-size-matched random regulators already reach
  > **~3×** and SAGA is only marginally above (**5.4×**, p=0.04 on the uniform-threshold method). SAGA's
  > real distinction is the **consistency of the de-repressive direction (95% ISG-up on KD vs 76%
  > random)**, not a unique fold. So the earlier "19–24×" headline is *not* strong SAGA-specific
  > evidence — read it as a prominent but largely general program that SAGA consistently restrains.
- **Programs are condition-dependent.** Across 3 conditions (phase-2 comparison), TCR programs are
  **stimulation-gated** (~4× Rest→Stim) while chromatin programs are **constitutive**; the interferon
  program is stimulation-gated in every class.
- **Disease link (one hypothesis, not a headline).** The module overlaps the clinically-tracked
  lupus / interferonopathy IFN signature and 31/163 of its genes are autoimmune GWAS risk genes — but
  an IFN module overlapping a clinical IFN signature is confirmatory by construction, so we state it as
  a single untested hypothesis: *coactivator knockdown de-represses ISGs; whether this is a control
  point in disease is not tested here.* Details in `docs/disease_and_specificity.md`.

These are **candidate convergent-target programs** (ISG-flagged), not causal pathways. UI: the
**Programs by class** tab (per-class cards, Jaccard heatmap, ISG-flagged target lists). Tables:
`class_isg_enrichment.csv`, `class_convergent_targets.csv`, `chromatin_stress_control.csv`,
`module_disease_overlap.csv`, `module_gwas_hits.csv`, `convergent_module_*`, `phase2_*`. Figures 26–31.
Detail: [`docs/disease_and_specificity.md`](docs/disease_and_specificity.md) ·
[`docs/literature_positioning.md`](docs/literature_positioning.md).

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
| `docs/tables/class_isg_enrichment.csv` · `class_convergent_targets.csv` | per-class convergent programs + interferon test |
| `docs/tables/convergent_module_*` · `phase2_*` | interferon module + condition-dependence analyses |
| `docs/tables/robust_edges.csv` | uncertainty-aware effect network (bonus, Model 1) |
| `docs/figures/*.png` | EDA · ranking · programs (20–25) · class programs (26–28) · overview |
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
**`/programs/classes`**, `/programs/class-targets?class=`, `/edges/downstream`. The UI has **6 screens**:
Overview, Explore, Audit, Effect network, Programs, and **Programs by class**.
Detail in [`api/README.md`](api/README.md).

## Limitations

- **Honest naming**: the models are **empirical-Bayes / pseudo-Bayesian**, not a full hierarchical NB
  or MCMC (no PPL, no formal random effects).
- `xcond_reproducibility` is an **exploratory feature** (cross-condition stability). It is **audited**
  with a **guide-level sensitivity analysis (partial donor coverage)** (`make repro-meta` →
  `scripts/extract_de_obs_metadata.py`, which extracts only the `.obs`, no `.layers`) that
  **reweights** the EB score with **real** cross-guide (`guide_correlation_all`) and cross-donor
  (`donor_correlation_hits_mean`) reproducibility — it is a **sensitivity analysis, not a new
  posterior**. Optional, and it does not replace the core; see *Sensitivity audit* in
  `docs/report.md`. Honest coverage: cross-guide ~78% but cross-donor only ~19% of contrasts (the
  KD-gated subset) → a **neutral weight** where it's missing (a gene is not penalized for lacking
  donor metadata). *This is the old sensitivity audit; the dedicated per-donor object
  (`by_donors.h5mu`, Tier 1) later gave **100%** cross-donor coverage — see
  [`docs/DONOR_ROBUSTNESS.md`](docs/DONOR_ROBUSTNESS.md).*
- **Model 1 is optional**: remote per-slice access is viable (~4.5 s/row, measured) but latency-bound;
  the official deliverable stands on its own with the local core.

## Submission summary

A reproducible product that turns a 1.8 TB CD4 Perturb-seq screen into **three explorable objects —
robust regulators, reproducibility audits, and transcriptional programs** — runnable on a laptop with
~10 GB of disk using only 15 MB of data for the core. The **CSV-only** core ranks regulators with
uncertainty (empirical Bayes) and audits them (bootstrap stability + a guide-level reproducibility
audit with partial donor coverage). On top of that, **fingerprint similarity organizes the top perturbations into
recognizable programs** — recovering the SAGA, Mediator and TCR complexes (permutation z=9/3/11,
robust to a within-condition null) and
surfacing candidate neighbors (e.g. CHD7 assigned to the chromatin program by fingerprint) — plus a
bonus **uncertainty-aware effect network**, both
streamed from the 17 GB h5ad without downloading it. An explorable **Regulator Atlas** (read-only
FastAPI + UI) ties it together: search a gene and see its rank, audit survival, transcriptional
program, transcriptomic neighbors, and defining response genes in one view. `make all` reproduces the
core; `make fingerprints` builds the programs; `make api` launches the atlas. See `docs/report.md`.

---
Data: CZI Virtual Cells Platform · Marson Lab 2025 · biorxiv preprint `10.64898/2025.12.23.696273`.
Original analysis code: [emdann/GWT_perturbseq_analysis_2025](https://github.com/emdann/GWT_perturbseq_analysis_2025).
