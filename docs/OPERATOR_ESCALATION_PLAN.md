# Operator escalation plan — expanding the regulator axis

**Audience:** an engineer picking this up cold. You do not need to have followed the earlier
discussion. Everything you need — context, commands, the one trap that already bit us, and how to
know when you're done — is below.

---

## 1. What this is, in one paragraph

We built an "empirical regulatory operator": a 3-way tensor `T[regulator, gene, condition]` holding
the transcriptional effect of each CRISPRi knockdown, in **precision-decoupled z-score space**. A
**pilot** ran it on **800 regulators** and produced two results worth escalating: (i) three
condition-gated CP factors with bootstrap CIs (proves the representation is genuinely pooled across
conditions), and (ii) a low-rank completion model that **predicts held-out regulators' late-stim
responses better than a persistence baseline** (R² 0.154 vs −0.176, on 160 regulators outside the
original panel). The full method and results are in [`docs/OPERATOR_ANALYSIS.md`](OPERATOR_ANALYSIS.md).
**Your job: run the same pipeline on a larger regulator axis** to see how far the predictive result
extends. This document is the how.

---

## 2. The one trap that will waste your time if you skip this

There is a **fail-closed confound guard** in `scripts/build_operator_tensor.py`. It refuses to build
the tensor (`[REFUSE-TO-CACHE]`, non-zero exit) if the regulator selection re-introduces a known
statistical confound: knockdowns done in **fewer cells** produce **larger** z-score norms (noise, not
signal), so `spearman(‖slab‖, n_cells)` must stay `|ρ| < 0.15`. If you pick regulators badly, the
leading factors of the operator become noise-magnitude axes and every downstream result is
meaningless. The guard exists to stop exactly that.

**Two facts about the guard you must internalize before choosing an axis size:**

1. **The power floor alone does NOT make the tensor guard-clean.** The selection rule has two parts —
   a median cell-count floor *and* ranking by `n_downstream` (breadth). The floor removes the
   worst low-power tail, but the confound still grows as you add more regulators. At "all regulators
   above the floor" (~3106) the *raw* confound is high and rising (≈ −0.28 on panel genes, −0.36 on
   all genes) — so 3106 is **very likely to fail the guard**, but note we established this from the
   raw numbers only; we never built the 3106 tensor and ran the z-guard on it. (This is consistent
   with fact 2 below: raw predicts *roughly*, the guard decides. Don't take "3106 fails" as measured —
   build it and let the guard rule if you want the real answer.) It is the breadth-ranking to a
   *top-N* that keeps the confound down. So the question is never "floor or no floor," it is
   **"what is the largest top-N that still passes the guard."**

2. **You cannot predict the guard result from raw log-FC numbers. Do not try.** We have exactly one
   calibration point: at N=800, the raw-space confound was ≈ −0.03 to −0.09, but the *actual z-space*
   confound the guard measured was **+0.08** — a sign flip, not a scaling. Any table that extrapolates
   "raw ρ × some factor → predicted z-ρ" is unreliable (it mis-predicts even that one point we can
   check). **The guard on the built tensor is the only authority.** The good news: running the guard
   is cheap — one rebuild, it reads ρ and either caches or refuses. So we *search* with the guard
   instead of guessing.

---

## 3. The decision you make before touching a keyboard

There are two genuinely different escalations. Pick one consciously; they have different payoffs and
different costs.

**Option A — "largest guard-clean panel, no new code" (recommended first).**
Sweep the regulator count upward, let the guard find the ceiling, run the pipeline there. Zero new
method — just a larger `--n-total`. Fast (local, hours at most). The honest limitation: the
guard-clean ceiling is somewhere **between 800 and ~2500** regulators (we don't know exactly — that's
what the sweep finds). This is a **modest** step up from the 800 pilot, not the whole genome.

**Option B — "the real 6209, with a new method" (only if A plateaus and you need more).**
To include the low-power tail the floor removes (getting toward the full 6209 regulators), you must
replace the crude power floor with **per-cell inverse-variance weighting using real `lfcSE`**. That
data is **not local** — it needs a per-cell-layer fetch (same fetch class as
`scripts/fetch_fingerprint_matrix.py`, and as the donor step needs). This is a real implementation
task, not a flag change. Only take it if Option A's completion curve is still improving at its ceiling
*and* you specifically want the sub-floor regulators.

**Recommendation: do Option A now.** It answers "does the predictive result extend to more
regulators" cheaply and honestly. Treat Option B as a separate, later project with its own plan.
Do **not** start Option B speculatively — the point of A is to find out whether B is even worth it.

---

## 4. Option A, step by step

All commands run from the repo root. The Python interpreter is the repo venv
(`.venv/bin/python`, which the Makefile calls `$(PY)`). Each `build_operator_tensor.py` run
**re-fetches the z-score slice once** (~0.15 GB per 1000 regulators) to
`data/cache/panel_zscore_<hash>.npy` and caches it; repeat runs at the same size are offline.

### Step 4.1 — Find the guard-clean ceiling (the actual new work)

Rebuild the tensor at increasing `--n-total`, reading the guard's verdict each time. **Start at 2000**
(higher than the earlier over-cautious guess of 1200 — the guard, not a prediction, decides):

```bash
.venv/bin/python scripts/build_operator_tensor.py --n-total 2000 --top-genes 2000
```

- **If it prints a tensor summary and exits 0** → 2000 passed. Try larger: rebuild at `--n-total 2500`.
- **If it exits with `[REFUSE-TO-CACHE] confound guard: |rho|=… >= 0.15`** → that size failed. Step
  back: rebuild at `--n-total 1600`, then `1200` if needed, until one passes.

You are binary-searching for the largest `--n-total` that passes. Expect 2–4 rebuilds. The **last
size that passed** is your escalation panel; note its regulator count and the reported
`n_cells_confound_rho` (it's written to `docs/tables/operator_tensor_summary.json`). Leave the tensor
built at that size — the next steps read `data/cache/operator_tensor.npz`.

> Sanity check after it passes: open `docs/tables/operator_tensor_summary.json` and confirm
> `passed_confound_guard: true`, `representation: "pooled_zscore"`, and that `n_new_regulators`
> (regulators outside the original fingerprint 200) is a large majority of `n_regulators` — the
> out-of-panel prediction test in Step 4.3 depends on those.

### Step 4.2 — Recompute the SVD programs (cheap, no tuning)

```bash
.venv/bin/python scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate
```

Writes `docs/tables/operator_svd_{programs,enrichment,power}.csv` and
`docs/figures/32_operator_svd_scree.png`. **Check `operator_svd_power.csv`:** every PC's
`|spearman(u_k, n_cells)|` must stay small (pilot max was 0.121). If a leading PC is now
power-confounded, the larger panel re-introduced the problem the guard is supposed to prevent — stop
and report it, don't paper over it.

### Step 4.3 — The completion sweep (the flagship — run rank well past 12)

```bash
.venv/bin/python scripts/operator_completion.py --max-rank 30 --holdout 0.2
```

Writes `docs/tables/operator_completion_condition.csv` and
`docs/figures/35_operator_completion_curve.png`. The columns that matter are **`r2_model_novel`**
(low-rank model, out-of-panel regulators only) vs **`r2_persistence_novel`** (the baseline: predict
Stim48hr := Stim8hr). **Why `--max-rank 30`:** the pilot's R² was still *rising* at rank 11 with no
plateau, so 12 was too short to see where it saturates. On a larger matrix it may need even more —
if R² is still climbing at rank 30, raise `--max-rank` further and re-run.

### Step 4.4 — The CP decomposition at full config

```bash
.venv/bin/python scripts/decompose_operator_cp.py --rank auto --scale-control rms --max-rank 12 --stab-subsample 800 --boot-n 100
```

Writes `operator_cp_{factors,stability,cosine}.csv` and `docs/figures/33,34_operator_cp_*.png`. The
pilot ran a reduced config (`--stab-subsample 250 --boot-n 50`) to fit local runtime; this is the
full one. **This is the slow step** — the split-half stability sweep dominates. If it exceeds ~10 min
locally at your panel size, that is the signal to dispatch this one step to cluster compute (see §6).

### Shortcut

If you'd rather not run them individually, edit the `--n-total` in the Makefile's `operator-tensor`
target to your chosen size, then:

```bash
make operator-tensor operator-svd operator-completion operator-cp
```

(Order matters — the last three all read the tensor `operator-tensor` writes.)

---

## 5. How to read the result — the plateau gate

This decides what you're allowed to claim, and it is the one place the run can silently overclaim.

Open `docs/tables/operator_completion_condition.csv` and look at `r2_model_novel` across the rank
sweep:

- **It beats `r2_persistence_novel` at the ranks tested** → the operator is **predictive** out-of-panel.
  This is the real, defensible claim: "low-rank structure predicts unseen late-stim responses for
  regulators never characterized." Report it.
- **It saturates (flattens) at some rank `r*`** → you may *additionally* claim the operator is
  **effectively rank-`r*`** (low-dimensional). Report `r*`.
- **It is still climbing at the maximum rank you swept** → you may **NOT** claim an effective rank.
  The honest statement is "predictive, but not low-dimensional at the ranks tested." Do not cite an
  effective-rank number in this case — that is the overclaim to avoid.

Predictivity and low-dimensionality are *different* claims. The pilot earned the first, not the
second. Your escalation must clear the plateau gate to earn the second.

---

## 6. When to move a step to the cluster

Everything here runs locally except possibly Step 4.4 (CP) at large panels. There is no need to
pre-emptively dispatch to a cluster. The rule: if a single step exceeds ~10 minutes locally — in
practice only the CP stability sweep at large `--n-total` — that step is the dispatch candidate. Steps
4.1–4.3 stay local regardless. Do not move the whole pipeline remote; only the one slow step, and only
if it's actually slow at your size.

---

## 7. Acceptance criteria — what "done" looks like

You are done with Option A when **all** of these hold:

1. `operator_tensor_summary.json` shows `passed_confound_guard: true` at the largest `--n-total` you
   could get past the guard, and you've recorded that regulator count and its ρ.
2. `operator_svd_power.csv` shows no leading PC is power-confounded (max `|ρ|` well under 0.15).
3. `operator_completion_condition.csv` has `r2_model_novel > r2_persistence_novel` at the tested
   ranks (predictivity holds), and you've applied the §5 plateau gate to decide whether an
   effective-rank claim is licensed.
4. `operator_cp_{factors,stability,cosine}.csv` written at full config, with ≥1 clean gated factor
   (bootstrap CI excludes flat) — this re-confirms pooling at the new size.
5. One paragraph added to `docs/OPERATOR_ANALYSIS.md` reporting: the escalation panel size, the guard
   ρ, the completion R²-vs-persistence at the plateau (or "still climbing"), and — explicitly —
   whether the predictive result **extended, held, or weakened** relative to the 800-pilot's 0.154.

That last point is the scientific payoff. The whole reason to escalate is to answer "does the operator
stay predictive as the axis grows." Whatever the answer — extends / holds / weakens — it is a result,
and reporting it honestly (including "weakened") is the job.

---

## 8. Failure modes, and what each means

| Symptom | Cause | What to do |
|---|---|---|
| `[REFUSE-TO-CACHE]` at your chosen N | breadth-ranked confound exceeds 0.15 at that size | step down `--n-total`; the guard is right, your N was too big |
| Guard passes but a leading SVD PC is power-confounded | subtle power leak the guard's slab-norm metric didn't catch | stop, report; do not interpret that PC as a program |
| Completion R² still rising at rank 30 | rank sweep too short, OR operator genuinely not low-dimensional | raise `--max-rank`; if still climbing, report "predictive, not low-dimensional" (no effective-rank claim) |
| CP stability sweep runs > 10 min | large panel × `--boot-n 100` | dispatch that one step to cluster (§6) |
| Want to go past the guard ceiling to 6209 | the floor discards the low-power tail | that's Option B — needs real `lfcSE` per-cell weighting, a non-local fetch, a separate plan |

---

## 9. What NOT to do

- **Do not weaken or remove the confound guard to fit a bigger panel.** If a large panel trips the
  guard (3106 very likely will, given its raw confound), that is the guard working. If you want 6209,
  that's Option B (a new method), not a looser threshold.
- **Do not pick `--n-total` from a raw→z extrapolation table.** We tried; it mis-predicts the one
  point we can check. Let the guard decide.
- **Do not claim an effective rank without a plateau.** Predictive ≠ low-dimensional.
- **Do not start Option B speculatively.** Run Option A first; its result tells you whether B is worth
  building.
