#!/usr/bin/env python3
"""Reproducible re-render of the two convergence figures that were previously committed as
precomputed PNGs with no generating script — fig 26 (convergent interferon module) and fig 28
(phase-2 condition dependence). Threaded through the shared palette (scripts/_figstyle.py) so the
whole deck is colour-consistent, and fig 26 carries the same specificity-control caveat as fig 27.

Fig 26 panels A/B come from the committed module tables; panel C recomputes per-regulator ISG
edges from the local log_fc cache (same |log2FC|>log2(1.5) rule as the rest of the convergence work).
Fig 28 comes entirely from the committed phase-2 tables.

    python scripts/render_convergence_figures.py    # fully offline

Outputs: docs/figures/26_convergent_interferon_module.png · docs/figures/28_phase2_condition_comparison.png
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import _figstyle as S
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
LOG2_1P5 = np.log2(1.5)
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def fig26():
    genes = pd.read_csv(TAB / "convergent_module_genes.csv")
    summ = json.loads((TAB / "convergent_module_summary.json").read_text())
    regs = summ["robust_regulators"]
    isg_members = set(summ["ISG_members_in_module"])
    n_mod, n_isg = summ["module_size_ge4of6"], summ["ISGs_in_module"]
    bg_frac = summ["ISG_reference_in_universe"] / summ["universe_size"]
    mod_frac = n_isg / n_mod

    S.apply_rc()
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))

    # A — module genes by number of robust regulators that hit them
    vc = genes["n_robust_regulators"].value_counts().sort_index()
    axes[0].bar(vc.index.astype(int).astype(str), vc.values, color=S.SAGA, width=0.68, zorder=2)
    for xi, v in enumerate(vc.values):
        axes[0].text(xi, v + max(vc.values) * 0.01, str(int(v)), ha="center", fontsize=9, fontweight="bold")
    axes[0].set_xlabel(f"# of the {len(regs)} robust SAGA-family regulators hitting the gene")
    axes[0].set_ylabel("convergent module genes")
    axes[0].set_title(f"A · convergent module: {n_mod} genes hit by ≥4 of {len(regs)}", loc="left", fontsize=9.5)

    # B — ISG fraction: module vs measured background (the fold)
    fold = mod_frac / bg_frac
    axes[1].bar([0, 1], [bg_frac * 100, mod_frac * 100], color=[S.GENERIC, S.ISG], width=0.6, zorder=2)
    axes[1].set_xticks([0, 1]); axes[1].set_xticklabels(["measured\nbackground", "convergent\nmodule"], fontsize=9)
    axes[1].set_ylabel("% interferon-stimulated genes (ISG)")
    for xi, fr, lab in [(0, bg_frac, f"{bg_frac*100:.1f}%"), (1, mod_frac, f"{mod_frac*100:.0f}%")]:
        axes[1].text(xi, fr * 100 + mod_frac * 100 * 0.02, lab, ha="center", fontsize=10, fontweight="bold")
    S.callout(axes[1], f"{fold:.0f}×", xy=(1, mod_frac * 100 * 0.5), xytext=(0.5, mod_frac * 100 * 0.62),
              color=S.INK, fs=13)
    axes[1].set_ylim(0, mod_frac * 100 * 1.25)
    axes[1].set_title("B · module is interferon-enriched", loc="left", fontsize=9.5)
    axes[1].text(0, -0.20, "⚠ this magnitude is largely a general strong-perturbation effect — see the\n"
                 "specificity control (fig 29); the SAGA signal is the consistent de-repressive direction (panel C)",
                 transform=axes[1].transAxes, fontsize=5.8, color="#b03a2e", va="top")

    # C — per-regulator ISG edges among the module's ISG members (recomputed), coloured by direction
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv").reset_index(drop=True)
    obs["row"] = np.arange(len(obs))
    var = pd.read_csv(CACHE / "fingerprint_var.csv")["gene_name"].values
    lf = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
    isg_idx = {g: i for i, g in enumerate(var) if g in isg_members}
    gcol = "target_contrast_gene_name"
    # count at Stim8hr — the condition where the interferon program is active (phase-2); the module
    # is a stimulation-gated program, so peak-breadth conditions (e.g. Stim48hr) would understate it
    rows = []
    for g in regs:
        sub = obs[(obs[gcol] == g) & (obs["culture_condition"] == "Stim8hr")]
        if sub.empty:
            rows.append((g, 0, 0)); continue
        r = int(sub.iloc[0]["row"])
        y = np.asarray(lf[r], dtype=float)
        up = dn = 0
        for gene, i in isg_idx.items():
            if abs(y[i]) > LOG2_1P5:
                up += y[i] > 0; dn += y[i] < 0
        rows.append((g, up, dn))
    rd = pd.DataFrame(rows, columns=["reg", "up", "dn"]).sort_values("up")
    yp = np.arange(len(rd))
    axes[2].barh(yp, rd["up"], color=S.ISG, height=0.6, zorder=2, label="de-repressed (KD ↑ ISG)")
    axes[2].barh(yp, -rd["dn"], color=S.GENERIC, height=0.6, zorder=2, label="repressed (KD ↓)")
    axes[2].axvline(0, color="0.6", lw=0.8)
    axes[2].set_yticks(yp); axes[2].set_yticklabels(rd["reg"], fontsize=8)
    axes[2].set_xlabel(f"ISG edges at Stim8hr (of {n_isg} module ISGs)")
    axes[2].set_title("C · knockdown de-represses ISGs — consistently", loc="left", fontsize=9.5)
    axes[2].set_xlim(left=-5.5)
    axes[2].legend(fontsize=6.0, loc="lower left", frameon=False, handlelength=1.1,
                   borderaxespad=0.4)

    fig.suptitle("A convergent interferon module — SAGA-family knockdown de-represses interferon genes",
                 fontsize=11.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=2.4)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "26_convergent_interferon_module.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  fig 26 → 26_convergent_interferon_module.png  (module {n_mod}g, {n_isg} ISG, {fold:.0f}×)")


def fig28():
    p = pd.read_csv(TAB / "phase2_condition_comparison.csv")
    summ = json.loads((TAB / "phase2_condition_summary.json").read_text())
    sj = summ["self_jaccard_rest_vs_stim8hr"]
    byc = summ["by_class_isg_fold"]
    short = {"SAGA/chromatin": "SAGA", "Mediator": "Mediator", "TCR (context-specific)": "TCR",
             "Other robust": "Other"}
    S.apply_rc()
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
    xpos = {c: i for i, c in enumerate(COND_ORDER)}

    # A — downstream breadth across conditions: thin gene lines + bold class median (TCR steep, rest flat)
    for gene, sub in p.groupby("gene"):
        sub = sub.set_index("condition").reindex(COND_ORDER)
        cls = sub["reg_class"].iloc[0]
        axes[0].plot(range(3), sub["n_downstream"], color=S.CLASS_COLORS.get(cls, S.GENERIC),
                     alpha=0.28, lw=1, zorder=1)
    for cls, sub in p.groupby("reg_class"):
        med = sub.groupby("condition")["n_downstream"].median().reindex(COND_ORDER)
        axes[0].plot(range(3), med.values, color=S.CLASS_COLORS.get(cls, S.GENERIC), lw=2.6,
                     zorder=3, label=short.get(cls, cls))
    # label the two steepest TCR genes at their endpoints
    for gene in ("ZAP70", "CD3E"):
        g = p[(p["gene"] == gene) & (p["condition"] == "Stim8hr")]
        if not g.empty:
            axes[0].text(1.03, float(g["n_downstream"].iloc[0]), gene, fontsize=6.5, color=S.TCR, va="center")
    axes[0].set_xticks(range(3)); axes[0].set_xticklabels(COND_ORDER)
    axes[0].set_ylabel("downstream genes (breadth)")
    axes[0].set_title("A · TCR breadth is stimulation-gated;\nchromatin/Mediator constitutive", loc="left", fontsize=9.5)
    axes[0].legend(fontsize=6.5, loc="upper left", frameon=False)

    # B — per-class interferon fold across conditions (the interferon program is stim-gated everywhere)
    classes = [c for c in ["SAGA/chromatin", "Mediator", "Other robust", "TCR (context-specific)"] if c in byc]
    w = 0.26
    for j, cond in enumerate(COND_ORDER):
        vals = [byc[c].get(cond, 0) for c in classes]
        axes[1].bar(np.arange(len(classes)) + (j - 1) * w, vals, width=w,
                    color=S.CONDITION_RAMP[cond], label=cond, zorder=2)
    axes[1].axhline(1, color="0.6", lw=0.7, ls=":")
    axes[1].set_xticks(range(len(classes))); axes[1].set_xticklabels([short.get(c, c) for c in classes], fontsize=8)
    axes[1].set_ylabel("interferon (ISG) fold-enrichment")
    axes[1].set_title("B · interferon program peaks under stimulation", loc="left", fontsize=9.5)
    axes[1].legend(fontsize=6.5, frameon=False, title="condition", title_fontsize=6.5)

    # C — how much each gene's target set turns over Rest→Stim8hr (low Jaccard = context-specific)
    cls_of = dict(zip(p["gene"], p["reg_class"]))
    sjd = pd.DataFrame({"gene": list(sj), "jac": list(sj.values())})
    sjd["cls"] = sjd["gene"].map(cls_of)
    sjd = sjd.sort_values("jac")
    yp = np.arange(len(sjd))
    axes[2].barh(yp, sjd["jac"], color=[S.CLASS_COLORS.get(c, S.GENERIC) for c in sjd["cls"]], height=0.66, zorder=2)
    axes[2].set_yticks(yp); axes[2].set_yticklabels(sjd["gene"], fontsize=7)
    axes[2].set_xlabel("Jaccard(Rest, Stim8hr) target overlap")
    axes[2].set_title("C · TCR target sets turn over most\n(lowest Rest↔Stim overlap)", loc="left", fontsize=9.5)

    fig.suptitle("Phase-2 condition dependence — TCR programs are stimulation-gated, chromatin constitutive; "
                 "interferon peaks under stimulation", fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95), w_pad=2.2)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "28_phase2_condition_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  fig 28 → 28_phase2_condition_comparison.png")


if __name__ == "__main__":
    fig26()
    fig28()
