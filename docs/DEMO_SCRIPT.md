# Demo video script (≤3 minutes)

A tight walkthrough for the ≤3-min submission video. Two columns: **what to show** and **what to
say**. Total narration ~430 words ≈ 2:50 at a calm pace. Record the screen; a talking head is
optional.

**Before you record**
```bash
make all          # regenerate outputs (≈8 s)
make api          # serve API + UI on http://localhost:8000/
```
Open two things: the browser at `http://localhost:8000/` and `docs/report.md` (or the figures).

---

### 0:00–0:25 — The problem  ·  *show: title slide, then figure `01_distribution_n_total_de_genes.png`*
> "This is a genome-scale CRISPRi Perturb-seq screen of primary human CD4+ T cells — 34,000
> differential-expression contrasts. The question: which genes are *robust regulators* of T cell
> programs? The catch is that the signal is heavy-tailed — the median perturbation moves just two
> genes, 15% move none, but a handful of hubs move thousands. So ranking by raw counts is a trap."

### 0:25–0:50 — Compute-aware design  ·  *show: `README.md` top + `make all` running in a terminal*
> "The full dataset is 1.8 terabytes; my laptop had ten gigabytes free. So the entire core runs from
> the 15-megabyte supplementary tables alone. One command — `make all` — reproduces everything in
> about eight seconds: the EDA, the model, the audits, and the report."

### 0:50–1:20 — The model & the finding  ·  *show: figure `08_kd_gate_changes_ranking.png`, then `07_hub_posterior_ranking.png`*
> "Two ideas do the work. First, a knockdown gate: only 62% of contrasts have a real on-target
> knockdown, and those carry 85% of all downstream signal. Second, an empirical-Bayes model that
> shrinks each gene's effect toward zero unless it's large *and* stable across conditions. Watch what
> happens to the ranking: the raw TCR-signaling hubs get demoted, and what rises is chromatin and
> transcription machinery — the SAGA complex, Mediator, KDM1A, SETD2. Robust regulators, with
> uncertainty, not just big hubs."

### 1:20–1:45 — It survives real reproducibility  ·  *show: figure `19_reproducibility_aware_ranking_shift.png`*
> "I stress-tested that ranking against *real* cross-guide and cross-donor reproducibility, pulled
> from just the metadata of the 17-gigabyte file — no full download. It's a sensitivity analysis, not
> a new model. The SAGA and Mediator regulators survive; a few single-guide hits get demoted. Honest,
> and defensible."

### 1:45–2:40 — The product: Regulator Atlas  ·  *show: the live UI at localhost:8000*
> "All of that becomes an explorable atlas."
- *Overview tab:* "7,913 ranked regulators, global versus context-specific."
- *Explore → search `SGF29`:* "Search a gene and you get its full profile — regulatory power,
  bootstrap stability, per-condition effect, and whether it survived every audit."
- *Explore → search `ZAP70`:* "Or a context-specific TCR regulator — active only under stimulation."
- *Audit tab:* "The reproducibility audit, transparent, row by row."
- *Programs tab → neighbors of `SGF29`:* "And a fingerprint map: it independently recovers the SAGA,
  Mediator, and TCR complexes as statistically cohesive groups — biology the model was never told
  about."

### 2:40–2:55 — Close  ·  *show: README "Submission summary" or the pipeline overview figure*
> "So: a 1.8-terabyte screen turned into a laptop-reproducible, uncertainty-aware shortlist of robust
> regulators — and a UI to explore it. `make all` reproduces the science; `make api` launches the
> atlas. Thanks for watching."

---

**Recording tips**
- Keep the terminal font large; pre-run `make all` once so the on-camera run is fast and clean.
- For the UI, use a maximized browser window and move deliberately — hover, click, pause a beat on
  each screen so a judge can read it.
- If you run long, the two cuttable sections are 1:20–1:45 (reproducibility) and the `ZAP70` beat.
