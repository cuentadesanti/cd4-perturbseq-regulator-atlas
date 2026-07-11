# Narration script

Beat-level: each line is one `with self.beat(...)` block wrapping the animation it describes. Edit in `demo.py`; this mirrors it. Numbers from `docs/OPERATOR_ANALYSIS.md`.

- We begin with a genome-scale screen of primary human T-cells: thousands of gene regulators, each switched off in turn, while single-cell sequencing reads the response.
- Take one regulator. Inhibit it, and follow a single downstream gene.
- Across many cells, the screen compares control to perturbed — the whole distribution shifts.
- That shift collapses into one standardized effect — direction and evidence in a single number.
- Do that for every gene — thousands of standardized effects stack into one transcriptional fingerprint.
- One fingerprint per regulator — stack thousands into the regulator-by-gene map.
- Correlate every pair, and the map folds into one square object — the operator.
- A word on method. Every claim ran the same loop —
- question, hypothesis, executable test, artifact — each one audited adversarially.
- Concretely: is SAGA a stable community? The gate weakened it. Predict an unseen regulator? The control failed it.
- Now the questions. Inside the operator, which directions are real — and which are noise?
- Random-matrix theory draws a line. The closed-form Marchenko-Pastur edge admits three hundred and thirty-six directions.
- But that edge is optimistic. A permutation null — the same data with the signal shuffled out — pushes it outward, and the honest count drops to ninety-two.
- Real, reproducible structure — of which only about seven carry across held-out conditions.
- Cleaned and clustered, the operator falls into communities — eight; three stable enough to trust.
- Watch this one.
- Its members name themselves — the N-D-U-F genes, subunits of Mitochondrial Respiratory Chain Complex One.
- CORUM, an external database, confirms it — a false- discovery rate near one in ten million. The clustering used no annotations; we asked what it contained only afterward.
- Good science also reports what fails — so Claude argued against its own results.
- Take the SAGA module. Its subunits co-cluster — but as a stable community it misses the gate. It stands only as a convergent module, supported independently.
- The hardest test: can a regulator's own features predict a regulator we have never perturbed?
- No. Real features and shuffled features both land on zero. A clean null — reported, not buried.
- None of this asks you to take our word for it.
- The same code that runs the analysis writes the tables, the figures, and the manuscript itself.
- Seventeen tests, open source, and the whole core rebuilding in about eight seconds on a laptop.
- Which leaves one honest lesson. Recoverable structure is not the same as inductive predictability.
- Claude helped us find the structure — and, just as importantly, where the claim stops.
- The empirical regulatory operator of CD4 T-cells — built with Claude.

_28 beats · 431 words._


## Voice
- Default preview: macOS `say` (Samantha), offline.
- Target: ElevenLabs **Elise** (`EST9Ui6982FZPSi7gCHi`), wired as default `eleven` voice — blocked on the free plan (HTTP 402 on library voices); `./render.sh final eleven` after upgrading.
