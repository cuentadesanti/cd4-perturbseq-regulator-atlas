# Hackathon submission — form answers

> Copy-paste text for the submission form. Fields marked **[decide]** need your choice.
> Placeholders in ⟨angle brackets⟩ are things to fill in (repo URL, video URL).

---

## Team Name
⟨your team name⟩

## Team Members
Santiago Silva (cuentadesanti)

## Project Name
**Regulator Atlas — robust regulators of CD4+ T cell programs**

*(Alternatives if a shorter name is wanted: "CD4 Perturb-seq Regulator Atlas" · "Robust Regulators".)*

## Track — **[decide]**
Pick the track closest to *computational biology / applied ML on omics data*. This project is a
reproducible analysis + product built on a genome-scale CRISPRi Perturb-seq dataset, so a
**biology / life-sciences** or **data/ML** track fits best. (Set to the exact option the form offers.)

## Link to your work
- GitHub repo: ⟨public repo URL⟩
- Judge-facing report: `docs/report.md` · interactive data model: `docs/data-model.html`
- One command to reproduce the core: `make all` (runs in ~8 s from 15 MB of CSVs)
- One command to launch the product: `make api` → http://localhost:8000/

## Demo Video (≤3 min)
⟨YouTube/Loom URL⟩ — script and shot list in [`docs/DEMO_SCRIPT.md`](docs/DEMO_SCRIPT.md).

---

## Project description
*(what you built or investigated, what you found, and why it matters)*

**What we built.** We turn a massive CD4+ T cell Perturb-seq screen (Marson Lab, 2025 — 4 donors × 3
conditions, 33,983 differential-expression contrasts) into an atlas of **three explorable objects:
robust regulators, reproducibility audits, and transcriptional programs.** It is memory- and
compute-aware: the full dataset is 1.8 TB and the working laptop had ~10 GB free, so the entire core
runs from the **15 MB supplementary CSV tables alone**, and the heavier objects stream `.h5ad` layers
*by slice* from S3 without downloading them.

- **Robust regulators.** An **empirical-Bayes ranking** (fixed-effects negative-binomial GLM +
  normal-normal shrinkage) ranks regulators with posterior uncertainty, gated on a validated knockdown.
- **Reproducibility audits.** Baseline comparison, bootstrap stability, a global-vs-context-specific
  split, and a **guide/donor-aware audit** that reweights the ranking with real cross-guide/cross-donor
  reproducibility extracted from just the `.obs` of the 17 GB `.h5ad`.
- **Transcriptional programs.** A regulator's rank is one number; its **fingerprint** — its downstream
  effect vector — is its whole action on the cell. On a balanced panel of 200 top perturbations we
  match each fingerprint to the curated SAGA/Mediator/TCR complexes and label the recognizable programs.
- **Regulator Atlas.** A read-only FastAPI service + single-page UI: search a gene and see its rank,
  audit survival, transcriptional program, transcriptomic neighbors, and defining response genes at once.

**What we found.**
- Effects are heavy-tailed: the median perturbation moves 2 genes and 15% move none, but 1.5% are
  hubs (>1000 DEGs) — so percentiles and rankings beat means.
- The knockdown gate is decisive: the 62% of contrasts with a significant on-target knockdown carry
  **85%** of all downstream effects.
- Robust regulators ≠ raw hubs. After the gates and shrinkage, the top is **chromatin/transcription
  machinery** — the SAGA complex (TADA1/TADA2B/SGF29/SUPT20H), Mediator (MED12/CCNC), KDM1A, SETD2 —
  a large *and* stable effect across conditions, above the Stim8hr-specific TCR-signaling hubs.
- **Fingerprint similarity organizes the top perturbations into recognizable programs.** The three
  curated complexes are each significantly cohesive by permutation test (**TCR z=11, SAGA z=9, Mediator
  z=3**), and the latent axis is program *identity*, not magnitude (|PC1| vs. n_downstream = 0.25). 25
  of 200 regulators map to a program — **TCR signaling (13), SAGA/chromatin (9), Mediator/transcription
  (3)** — recovering each complex's core and adding new members, most strikingly the chromatin remodeler
  **CHD7 joining the SAGA/chromatin program** (cosine 0.84). The reproducibility-promoted hits have
  neighborhoods as tight as the top regulators (so they aren't noise) but map onto *none* of the
  canonical complexes — the audit surfaces a distinct coherent set, not "more SAGA".

**Why it matters.** It answers not just *who is a strong regulator* but *what program each perturbation
induces and who resembles whom* — turning a 1.8 TB screen into a defensible, uncertainty-aware atlas a
bench scientist can act on. Reproducible on a laptop, explorable gene-by-gene, and honest about its
scope (program similarity anchored to known complexes, empirical-Bayes rather than full MCMC).

---

## How did you use Claude? Which products, and where did they matter most?
*(honest draft — adjust to match exactly what you used)*

I built this end-to-end with **Claude Code** as the development environment and reasoning partner.
Where it mattered most:

- **Data modeling & scoping.** Claude helped reverse-engineer the dataset's star schema (cell →
  pseudobulk → DE statistics, and the two roles of the GENE table) and decide the 80/20 cut: run the
  whole core from the 15 MB CSVs and treat the 1.8 TB of `.h5ad`/`.h5mu` as optional, streamed-by-slice
  extras. That memory-aware framing is the backbone of the project.
- **Statistical modeling.** Claude helped design the empirical-Bayes ranking (a stable IRLS route to
  the negative-binomial mean, then normal-normal shrinkage of a per-gene log-rate effect) and, just
  as importantly, insisted on **honest naming** — calling it pseudo-Bayesian rather than a full
  hierarchical model, and separating the sensitivity audit from the core posterior.
- **The programs analysis.** Claude built the transcriptional-program layer and, at a review gate,
  caught that the first attempt (labeling agnostic clusters) merged SAGA/TCR/Mediator into one blob;
  it pivoted to a transparent nearest-known-complex-centroid classifier, kept the labels conservative
  (175/200 stay *mixed*), and reported the promoted-hit coherence result honestly even though it came
  back *null* against the canonical complexes.
- **Building the product.** Claude wrote the FastAPI service, the single-page Atlas UI, the audit and
  program pipelines, and the figure-generation code, and kept the reproducibility contract intact
  (`make all` is the source of truth; the API only serves versioned outputs).
- **Verification & polish.** Claude re-checked every quantitative claim against the actual tables,
  and did a full-repo pass to make the submission judge-ready.

*(If you used **Claude Science** / a specific Claude model or feature for any of the domain reasoning,
name it here — e.g. which model, and whether Claude Science surfaced the SAGA/Mediator biology.)*

---

## Thoughts / feedback on building with Claude Science
*(draft prompts — keep what's true for you)*

- The biggest win was **compute-aware research design**: being pushed to extract value from 15 MB
  instead of reaching for 1.8 TB made the whole thing reproducible on a laptop, and reviewers can
  actually re-run it.
- The **honesty guardrails** were valuable — being nudged to label the models "empirical-Bayes, not
  full-Bayes" and to keep the effect-network as a clearly-scoped bonus kept the claims defensible.
- ⟨add anything that surprised you — where it saved the most time, where you had to correct it, what
  you'd want next time⟩
