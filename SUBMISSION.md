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

**The question.** In a genome-scale CRISPRi Perturb-seq screen of primary human CD4+ T cells
(Marson Lab, 2025 — 4 donors × 3 conditions, 33,983 differential-expression contrasts), which genes
are *robust regulators* of T cell programs? The raw signal is dominated by noise and by a handful of
hubs, so we wanted to separate signal from noise **with uncertainty** and rank regulators by a
**large and reproducible** effect — not by raw differentially-expressed-gene counts or `p < 0.1`.

**What we built.** A reproducible, memory- and compute-aware research product. The full dataset is
1.8 TB and the working laptop had ~10 GB free, so the entire core runs from the **15 MB supplementary
CSV tables alone**:
1. An 80/20 EDA that characterizes the effect-size distribution and quality gates.
2. An **empirical-Bayes regulator ranking** (fixed-effects negative-binomial GLM + normal-normal
   shrinkage of a per-gene log-rate effect) that ranks regulators with posterior uncertainty.
3. A **ranking audit** — baseline comparison, bootstrap stability, and a global-vs-context-specific
   split — that shows *why* the quality gates matter.
4. A **guide/donor-aware sensitivity audit** that reweights the ranking with real cross-guide and
   cross-donor reproducibility extracted from just the `.obs` of the 17 GB `.h5ad` (no `.layers`).
5. A bonus **uncertainty-aware effect network** streamed by slice from the remote 17 GB `.h5ad`
   *without downloading it*, and a **transcriptional-fingerprint** side analysis (PCA + cosine
   similarity) that maps perturbations to programs.
6. A **Regulator Atlas** on top: a read-only FastAPI service + single-page UI to search a gene, see
   its full profile, filter global vs. context-specific, and browse every audit.

**What we found.**
- Effects are heavy-tailed: the median perturbation moves 2 genes and 15% move none, but 1.5% are
  hubs (>1000 DEGs) — so percentiles and rankings beat means.
- The knockdown gate is decisive: the 62% of contrasts with a significant on-target knockdown carry
  **85%** of all downstream effects.
- Robust regulators ≠ raw hubs. After the gates and shrinkage, the top is **chromatin/transcription
  machinery** — the SAGA complex (TADA1/TADA2B/SGF29/SUPT20H), Mediator (MED12/CCNC), KDM1A, SETD2 —
  a large *and* stable effect across conditions, above the Stim8hr-specific TCR-signaling hubs.
- The fingerprint space independently recovers known biology: SAGA, Mediator, and TCR complexes are
  each significantly cohesive by permutation test (TCR z=11, SAGA z=9, Mediator z=3), and PC1 is
  program identity, not effect magnitude (|PC1| vs n_downstream Spearman = 0.25).

**Why it matters.** It turns a 1.8 TB screen into a defensible, uncertainty-aware shortlist of
regulators that a bench scientist can actually act on — reproducible on a laptop, honest about its
limits (empirical-Bayes, not full MCMC; partial cross-donor coverage), and explorable through a UI
rather than a static table.

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
- **Building the product.** Claude wrote the FastAPI service, the single-page Atlas UI, the audit
  pipeline, and the figure-generation code, and kept the reproducibility contract intact (`make all`
  is the source of truth; the API only serves versioned outputs).
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
