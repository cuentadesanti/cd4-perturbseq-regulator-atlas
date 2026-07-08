# The interferon module — specificity control and the medical bridge

Two analyses probe the convergent interferon module (the 163 downstream genes hit by ≥4 of the 6
robust SAGA-family regulators). **Analysis 2 was run first** because it decides how strongly the
medical framing in Analysis 1 can be stated. Reproduce with `make convergence-extras`.

---

## Analysis 2 — the specificity control (decisive)

**Question:** is interferon (ISG) de-repression a SAGA property, or does *any* strong
transcription/chromatin perturbation de-repress interferon (generic stress)?

Three effect-comparable regulator groups, identical pipeline (convergent targets at
`|log₂FC|>log₂1.5`, self-edges removed; ISG hypergeometric test, 48-gene core, 10,282-gene
background), from the fully-local `log_fc` cache (`scripts/analyze_chromatin_stress_control.py`,
fig 29):

| group | n | median n_downstream | ISG fold | ISG-up on KD (direction) |
|---|---|---|---|---|
| **SAGA/chromatin** | 8 | 4,138 | **5.4×** (p=0.042 vs random) | **95%** de-repressive |
| other chromatin/transcription | 14 | 150 *(weak — see caveat)* | 2.9× | 44% (mixed) |
| **random (effect-size matched)** | 8 | 2,292 | **3.1×** (null, 500 draws) | 76% |

**Honest verdict — specificity is partial:**
- **The magnitude is largely generic.** Effect-size-matched random strong regulators already reach
  **~3.1×** ISG enrichment (their convergent targets are the most stimulation-responsive genes, which
  include ISGs). SAGA's 5.4× is only **marginally above** that (z=1.8, p=0.042). So the large "19–24×"
  headline from the earlier class-programs analysis reflects, in good part, a **general property of
  strong perturbations under stimulation**, not a SAGA-unique magnitude. *(That earlier number also
  used a stricter EB-gated edge set, which shrinks the target set and inflates the fold; here the
  uniform-threshold method puts every group on the same footing.)*
- **SAGA's real distinction is the *direction*.** SAGA knockdown de-represses ISGs **95%** of the time
  — more consistently than random (76%) and far more than the heterogeneous chromatin set (44%, which
  includes Polycomb repressors). SAGA-family coactivators **restrain** the interferon program.
- **Caveat on the other-chromatin control:** in this screen most non-SAGA chromatin/transcription
  genes knock down *weakly* (median 150 downstream vs SAGA's 4,138), so a strong, effect-matched
  "specific machinery" control isn't buildable — the **matched-random null is the fair comparator**.

**Takeaway:** the interferon axis is a prominent, *largely general* convergent readout of strong
perturbations in stimulated CD4⁺ T cells; SAGA-family coactivators are distinguished by **consistently
de-repressing** it, not by a unique enrichment magnitude.

---

## Analysis 1 — the medical bridge (nomination, framed by the control)

**Do the module's interferon genes coincide with the clinically-tracked, drugged autoimmune IFN
axis, and does the map nominate upstream control points?** (`scripts/analyze_disease_overlap.py` +
`scripts/analyze_module_gwas.py`, figs 30–31; background = the 10,282 measured genes.)

**Disease-signature overlap (confirmatory by construction, but IFN-specific):**

| signature | overlap | fold | p |
|---|---|---|---|
| Lupus type-I IFN (21-gene core) | 8/20 | 25× | 3.6e-10 |
| Interferonopathy score (AGS/SAVI) | 2/8 | 16× | 6.6e-3 |
| Type-I IFN response (canonical core) | 15/38 | 25× | 5.9e-18 |
| **T-cell exhaustion (contrast)** | **0/14** | **0×** | ns |

The module overlaps the lupus/interferonopathy IFN axis strongly (**expected** — its ISG core *is*
that clinical signature), and the **exhaustion contrast is not enriched**, so the module is
interferon-specific, not generically "inflammatory." The **payload is the direction**: module edges
are de-repressive (91% positive; ISG edges 100% positive) → the SAGA/Mediator coactivators **restrain**
the clinically-tracked interferon axis, which **nominates** them as upstream control points.

**Genetics corollary (Open Targets, autoimmune SLE/RA/MS/Sjögren/T1D):** **31/163** module genes are
autoimmune genetic-association (GWAS) risk genes — headlined by **STAT4** (associated with all 5
diseases, score 0.82; a canonical SLE/RA risk gene), plus IKZF3, RASGRP3, ANXA1, NOD2, and the ISGs
OAS1/OAS3. The SAGA regulators themselves carry modest autoimmune genetic signal (TADA2B 0.42,
SUPT20H 0.30). A coactivator that **restrains** an autoimmune-risk interferon gene is a **mechanistic
hypothesis** worth testing — nothing more.

---

## Honest bottom line

- The interferon/ISG axis is a real, IFN-specific convergent program that SAGA-family coactivators
  **de-repress**, and it maps cleanly onto the clinically-tracked lupus/interferonopathy signature and
  autoimmune GWAS genes (STAT4 chief among them).
- **But the *magnitude* of ISG enrichment is not SAGA-specific** — it is largely a general
  strong-perturbation-under-stimulation effect (matched random ~3×). SAGA's distinguishing feature is
  the **consistency of de-repression**, not a unique fold.
- Therefore this is framed as **nomination, not discovery**: SAGA/Mediator coactivators are candidate
  de-repressive control points of a disease-relevant interferon program; the disease/GWAS overlaps are
  **expected/confirmatory**; and no novel disease association is claimed from overlap alone.

Outputs: `chromatin_stress_control.csv` · `module_disease_overlap.csv` · `module_gwas_hits.csv`
(+ `*_summary.json`); figures 29–31.
