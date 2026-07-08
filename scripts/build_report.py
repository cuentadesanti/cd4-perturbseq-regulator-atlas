#!/usr/bin/env python3
"""Generates the consolidated judge-facing report and the pipeline overview figure.

    python scripts/build_report.py

Produces:
    docs/figures/00_pipeline_overview.png
    docs/report.md   (consolidates DATA_MODEL.md + EDA.md + MODELING.md + figures + top table)
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "docs" / "figures"
TAB = ROOT / "docs" / "tables"
FIG.mkdir(parents=True, exist_ok=True)

ACCENT, AMBER, VIOLET, INK, MUT = "#0a8f9c", "#b47818", "#6b5fc0", "#1b2130", "#586074"


def pipeline_overview():
    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.set_xlim(0, 100); ax.set_ylim(0, 34); ax.axis("off")

    def box(x, y, w, h, text, color, sub=""):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.4",
                    linewidth=1.5, edgecolor=color, facecolor=color + "18"))
        ax.text(x + w / 2, y + h / 2 + (1.2 if sub else 0), text, ha="center", va="center",
                fontsize=9.5, fontweight="bold", color=INK)
        if sub:
            ax.text(x + w / 2, y + h / 2 - 2.4, sub, ha="center", va="center", fontsize=7.2, color=MUT)

    def arrow(x0, x1, y, label=""):
        ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=13,
                    linewidth=1.4, color=MUT))
        if label:
            ax.text((x0 + x1) / 2, y + 2.1, label, ha="center", fontsize=6.8, color=MUT, style="italic")

    box(1, 10, 20, 14, "Supplementary CSVs", AMBER, "DE_stats · sgRNA · samples\n~15 MB · local")
    arrow(21, 29, 17, "features")
    box(29, 10, 18, 14, "80/20 EDA", ACCENT, "distribution, hubs,\nQC · scripts/eda.py")
    arrow(47, 55, 17, "ranking")
    box(55, 10, 20, 14, "Model 2 · EB", ACCENT, "robust regulators\nscripts/model_hubs.py")
    arrow(75, 83, 17, "candidates")
    box(83, 4, 16, 26, "Model 1 (optional)", VIOLET, "uncertainty-aware\neffect network\nh5ad · not downloaded")

    ax.text(50, 30.5, "Pipeline — from local CSV to regulators with uncertainty",
            ha="center", fontsize=12, fontweight="bold", color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "00_pipeline_overview.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("  figure → docs/figures/00_pipeline_overview.png")


def df_to_md(df):
    """Markdown table without depending on `tabulate`."""
    cols = list(df.columns)
    head = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = "\n".join("| " + " | ".join(str(v) for v in row) + " |"
                     for row in df.itertuples(index=False))
    return "\n".join([head, sep, body])


def section(path):
    """Returns the body of a .md without its first H1 title (for nesting)."""
    lines = (ROOT / "docs" / path).read_text().splitlines()
    out, skipped = [], False
    for ln in lines:
        if not skipped and ln.startswith("# "):
            skipped = True
            continue
        out.append(ln)
    return "\n".join(out).strip()


def build_report():
    review = pd.read_csv(TAB / "top_regulators_for_review.csv")
    top = review.head(15)[["rank", "gene", "condition", "regpower_eb_mean",
                           "p_top_1pct", "observed_n_downstream", "interpretation_note"]]
    tbl = df_to_md(top)

    edges_line = ""
    ep = TAB / "robust_edges.csv"
    if ep.exists():
        ed = pd.read_csv(ep)
        n, nreg = len(ed), ed["perturbed_gene"].nunique()
        poc = " It was assessed as a **proof-of-concept** (see `docs/EDGE_ANALYSIS.md`): coherent" \
              " and biologically sensible, but of minimal coverage — kept as a bonus, not a" \
              " strong result." if (ROOT / "docs" / "EDGE_ANALYSIS.md").exists() else ""
        edges_line = (f"\nWe also built an **uncertainty-aware effect network** (bonus, Model 1) "
                      f"with **{n:,} robust edges** (`P(|effect|>1.5×)>0.8` — the probability that the "
                      f"*magnitude* exceeds 1.5×, not that a causal edge exists) from {nreg} regulators in "
                      f"`docs/tables/robust_edges.csv`.{poc}\n")

    # --- audit sections (conditional on the audit existing) ---
    audit_md = ""
    comp_p = TAB / "ranking_baseline_comparison.csv"
    gtab_p = TAB / "top_global_regulators.csv"
    if comp_p.exists() and gtab_p.exists():
        comp = pd.read_csv(comp_p)
        top30 = comp.sort_values("raw_rank").head(30)
        n_drop = int(top30["dropped_by_kd_gate"].sum())
        demoted = top30[(~top30["dropped_by_kd_gate"]) & (top30["eb_rank"] > 30)]
        g_tab = pd.read_csv(gtab_p)
        c_tab = pd.read_csv(TAB / "top_condition_specific_regulators.csv")
        glist = ", ".join(g_tab.head(6)["gene"])
        clist = ", ".join(c_tab.head(6)["gene"])
        audit_md = f"""
## Naive hubs vs. quality-aware regulators

Ranking by raw `n_downstream` rewards hubs that don't survive the quality controls.
Of the top 30 raw hubs: **{n_drop} drop at the KD gate** (no validated on-target knockdown)
and **{len(demoted)} are demoted by EB shrinkage** for being condition-specific (the signal lives in
a single condition). The EB ranking surfaces regulators with a large **and** stable effect.

![gate](figures/08_kd_gate_changes_ranking.png)

Stability was audited with a bootstrap (B=200) over the eligible rows: how often each gene falls
in the top-30 (`stability_frequency`) is in `top_regulators_for_review.csv`. The ranking is
moderately stable — best read as a *set* of robust regulators, not an exact ordering.

![stability](figures/11_ranking_stability.png)

## Global vs. context-specific regulators

Splitting by `condition_specificity = max/sum of n_downstream` across conditions with significant KD:

- **Global** (stable effect in ≥2 conditions): {glist}… — chromatin/transcription machinery.
- **Context-specific** (effect concentrated in one condition): {clist}… — includes TCR signaling
  (ZAP70, LCK), active only under stimulation.

Both classes are real biology; the distinction avoids confusing a universal regulator with a
context-dependent one. Tables: `top_global_regulators.csv`, `top_condition_specific_regulators.csv`.

![globalvs](figures/10_global_vs_context_specific.png)
"""

    # --- guide/donor-aware sensitivity audit (real metadata, optional) ---
    repro_p = TAB / "reproducibility_audit.csv"
    cov_p = TAB / "reproducibility_coverage.csv"
    if comp_p.exists() and repro_p.exists() and cov_p.exists():
        au = pd.read_csv(repro_p)
        cov = pd.read_csv(cov_p).iloc[0]
        dem = au[au["status"] == "demoted"].sort_values("old_rank")
        pro = au[au["status"] == "promoted"].sort_values("new_rank")
        surv = au[au["status"] == "survives"].sort_values("new_rank")
        survlist = ", ".join(surv.head(5)["gene"])
        # interpretable table: demoted + promoted
        mv = pd.concat([dem, pro]).sort_values(["status", "old_rank"])
        mv = mv[["gene", "old_rank", "new_rank", "status", "guide_corr", "donor_corr", "reason"]]
        mv_tbl = df_to_md(mv)
        audit_md += f"""
## Sensitivity audit: does the ranking survive real (guide/donor) reproducibility?

The core ranking's caveat was using `xcond_reproducibility` (a cross-condition proxy). We audit it with
**real** reproducibility from the `.obs` of `DE_stats.h5ad` (`scripts/extract_de_obs_metadata.py`, `.obs`
only, ~4 s, no `.layers`): `guide_correlation_all` (agreement between the 2 guides) and
`donor_correlation_hits_mean` (cross-donor agreement), plus a penalty for single-guide targets.

**We do not re-estimate the EB posterior and it is not a new model**: we reweight the EB score
(`reweighted_score = regpower_eb_mean · repro_weight`) as a **sensitivity analysis** — which
regulators survive. In practice it is more **guide-aware** than **donor-aware**: `guide_corr` covers
{cov.pct_guide_corr_available:.0f}% of contrasts but `donor_corr` only **{cov.pct_donor_corr_available:.0f}%**
(the cross-donor analysis was done on a subset); where it's missing, a **neutral weight** is used (no
penalty for absent data). Of the top-30 EB, **{len(dem)} are demoted** and **{len(pro)} are promoted**:

{mv_tbl}

- **Survivors** (large effect + reproducible): {survlist}… — SAGA + Mediator.

![reproshift](figures/19_reproducibility_aware_ranking_shift.png)

**The core ranking still works without this file** — the audit lives separately in
`hub_ranking_bayes_reproducibility_aware.csv` / `reproducibility_audit.csv`.
"""

    # --- transcriptional programs (optional, `make fingerprints`) ---
    programs_md = ""
    fnd_p, cval_p, evid_p = (TAB / "fingerprint_findings.csv", TAB / "fingerprint_complex_validation.csv",
                             TAB / "program_label_evidence.csv")
    if fnd_p.exists() and cval_p.exists() and evid_p.exists():
        import json
        fnd = pd.read_csv(fnd_p)
        cval = pd.read_csv(cval_p)
        evid = pd.read_csv(evid_p)
        sp = TAB / "fingerprint_summary.json"
        summ = json.loads(sp.read_text()) if sp.exists() else {}
        coh, pc1 = summ.get("audit_coherence", {}), summ.get("pc1_abs_vs_ndownstream_spearman")
        n_in = int((fnd["program_label"] != "mixed").sum())
        counts = fnd[fnd.program_label != "mixed"]["program_label"].value_counts().to_dict()
        counts_txt = ", ".join(f"{k} ({v})" for k, v in sorted(counts.items(), key=lambda x: -x[1]))
        # donor-robustness of the assigned set (Tier 1 fold-in): how many assignments fail the donor check
        assigned = fnd[fnd.program_label != "mixed"]
        if "donor_robust" in fnd.columns:
            frag_by = assigned[assigned.donor_robust == False].groupby("program_label").size().to_dict()
            n_frag = int(sum(frag_by.values()))
            frag_txt = ("; ".join(f"{k} {v}" for k, v in sorted(frag_by.items(), key=lambda x: -x[1]))
                        or "none")
        else:
            n_frag, frag_txt = 0, "none"
        zline = " · ".join(f"{r.complex} z={r.z_cross}→{r.z_within} (cross→within-condition)"
                            for r in cval.itertuples() if pd.notna(r.z_within))
        ev = evid[evid.program_label != "mixed"][["program_label", "n_regulators", "n_known_complex_members",
                                                  "assigned_neighbors", "mean_centroid_cosine", "top_marker_genes"]]
        ev_tbl = df_to_md(ev)

        def _c(s):
            d = coh.get(s, {})
            return f"{d.get('mean_knn_sim', '?')}" if d else "—"
        programs_md = f"""
## Transcriptional programs

A rank is one number; a **fingerprint** — a regulator's downstream effect vector — is what the
perturbation actually does to the cell. On a balanced panel of 200 top perturbations we match each
regulator's fingerprint to the curated **SAGA / Mediator / TCR** complexes (nearest-centroid in the
same space as the validated cosine similarity). These are **candidate program assignments by
fingerprint similarity — not claims of physical complex membership.** The classifier is conservative:
only **{n_in} of 200** perturbations are assigned a program (**{n_frag} flagged donor-fragile**: {frag_txt});
the rest remain *mixed*, by design.

- **Fingerprint similarity recovers the known complexes** (permutation test, N=5000): {zline}. The
  latent PC1 is program *identity*, not effect magnitude (|PC1| vs. n_downstream Spearman = {pc1}).
- **{n_in} assigned** ({n_frag} donor-fragile — {frag_txt}): {counts_txt}. Each program recovers its
  curated core and adds **newly assigned neighbors** (non-curated genes placed in the same fingerprint
  neighborhood) — e.g. the chromatin remodeler **CHD7** is assigned to the SAGA/chromatin program (a
  related perturbation response, not complex membership) and is **donor-robust**; the donor check flags
  ATF7IP2/NCAPG2/EIF1AX (TCR) and GLIPR2 (Mediator) as fingerprint artifacts that do not replicate
  across donors (see `donor_fragile_neighbors` in `program_label_evidence.csv`).

{ev_tbl}

![programs](figures/24_fingerprint_pca_by_program.png)
![neighbors](figures/23_fingerprint_neighbor_network.png)

**Do the reproducibility-promoted hits form coherent programs?** They have neighborhoods as tight as
the top global regulators (mean kNN cosine: promoted {_c('promoted')}, demoted {_c('demoted')} vs.
global {_c('global')}), so they are **not statistical noise** — yet they map onto *none* of the
canonical complexes. Read as: the audit surfaces a **distinct high-confidence set** rather than simply
rediscovering the known complexes.

*Scope: fingerprint-based, program-level re-analysis anchored to known complexes — candidate
assignments and hypotheses, not de-novo pathway discovery or novel complex membership. "Response
genes" are consistently-moved downstream genes (relative to the panel), not baseline markers; PCA is a
view, not the proof. `make fingerprints` · detail in `docs/FINGERPRINT_ANALYSIS.md`.*
"""

    # --- convergent programs by regulator class (optional, `make class-programs`) ---
    class_md = ""
    cls_p = TAB / "class_isg_enrichment.csv"
    if cls_p.exists():
        ci = pd.read_csv(cls_p).sort_values("fold", ascending=False)
        saga = ci[ci["class"] == "SAGA/chromatin"]
        saga_fold = float(saga["fold"].iloc[0]) if len(saga) else None
        jm = None
        csp = TAB / "class_program_summary.json"
        if csp.exists():
            import json as _j
            jm = _j.loads(csp.read_text()).get("_jaccard_offdiag_median")
        cm = TAB / "convergent_module_summary.json"
        mod = {}
        if cm.exists():
            import json as _j
            mod = _j.loads(cm.read_text())
        ct = ci[["class", "targets", "ISGs", "fold", "p", "frac_up_on_KD"]].copy()
        ct["fold"] = ct["fold"].map(lambda x: f"{x}×")
        ct["p"] = ct["p"].map(lambda x: f"{x:.1e}")
        class_md = f"""
## Convergent programs by regulator class

Is "chromatin machinery recovers as top hubs" just the expected result of perturbing coactivators? A
**balanced 30-regulator panel** — chosen by *class*, not by rank — tests whether classes converge on
*distinct* downstream programs (`make class-programs`, fully offline).

- **Classes converge on distinct programs**: median off-diagonal Jaccard of per-class convergent-target
  sets ≈ **{jm:.2f}** — classes barely share targets.
- **A convergent interferon module** answers the "SAGA is expected" critique: genes hit by ≥4 of the 6
  robust SAGA-family regulators form a **{mod.get('module_size_ge4of6','?')}-gene module**,
  **{mod.get('ISG_fold_enrichment','?')}× enriched for interferon-stimulated genes**
  (P≈{mod.get('ISG_hypergeom_p', float('nan')):.1e}), **all de-repressive** (knockdown raises ISGs).
- Interferon repression is **most concentrated under SAGA/chromatin ({saga_fold}× in the class panel)**.

{df_to_md(ct)}

![class programs](figures/27_regulator_class_programs.png)

*Candidate convergent-target programs (ISG-flagged), not causal pathways. Detail:
`docs/literature_positioning.md`; per-class target lists in the **Programs by class** UI tab.*
"""

    md = f"""# Report — Genome-scale CD4+ T cell Perturb-seq

*Consolidated report for review. Reproducible with `make all` (local CSVs only).*

![pipeline](figures/00_pipeline_overview.png)

## Question

Which genes are **robust regulators** of CD4+ T cell programs — separating real
signal from noise and prioritizing by a **large and reproducible** effect, not by raw counts?

## Executive summary

- Perturbation effects are **heavy-tailed**: median 2 DEGs, but 1.5% are hubs
  with >1000. Summarize with percentiles and rankings, not the mean.
- **Effective knockdown gates the signal**: contrasts with a significant on-target KD (62%)
  concentrate **85%** of all trans-effects.
- An **empirical-Bayes** (pseudo-Bayesian) model ranks regulators by their latent regulatory
  power with uncertainty. The robust top is **chromatin/transcription** machinery
  (SAGA complex, Mediator, KDM1A, SETD2) — a large **and** stable effect across conditions.
- **Fingerprint similarity organizes the top perturbations into recognizable programs** — recovering
  TCR signaling, SAGA/chromatin and Mediator/transcription (permutation z=11/9/3) and surfacing
  candidate neighbors (e.g. the chromatin remodeler CHD7 assigned to the chromatin program by
  fingerprint, not by complex membership). Same honesty, different object: not just *who* is strong,
  but *what program* each perturbation resembles and *who resembles whom*.
{edges_line}
## Top regulators (for review)

{tbl}

Full table (30, with all columns): `docs/tables/top_regulators_for_review.csv`.

![ranking](figures/07_hub_posterior_ranking.png)
{audit_md}
{programs_md}
{class_md}
## EDA findings

![degs](figures/01_distribution_n_total_de_genes.png)
![hubs](figures/03_top_hubs_by_condition.png)
![kd](figures/04_ontarget_vs_downstream.png)
![repro](figures/06_reproducibility_vs_effects.png)

---

## Appendix A — Data model

{section("DATA_MODEL.md")}

---

## Appendix B — EDA

{section("EDA.md")}

---

## Appendix C — Modeling

{section("MODELING.md")}
"""
    (ROOT / "docs" / "report.md").write_text(md)
    print("  report → docs/report.md")


if __name__ == "__main__":
    print("== build_report ==")
    pipeline_overview()
    build_report()
    print("✓ Report generated.")
