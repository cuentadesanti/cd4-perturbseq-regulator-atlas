#!/usr/bin/env python3
"""Balanced 30-regulator class-program analysis (offline, local caches only).

Expands the Model-1 edge network from the 6 top-EB regulators to a *balanced*
30-regulator panel chosen by CLASS (SAGA/Mediator/TCR/other/promoted/demoted),
then asks whether different regulator classes converge on different downstream
programs. Runs with NO network: reads log_fc from the full local cache and
recovers lfcSE = log_fc / zscore from the cached fingerprint panel.

    python scripts/analyze_class_programs.py

Outputs (docs/tables/): balanced_panel_30.csv, robust_edges_balanced30.csv,
class_convergent_targets.csv, class_isg_enrichment.csv, class_program_summary.json
and figure docs/figures/27_regulator_class_programs.png
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm, hypergeom

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
LOG2_1P5 = np.log2(1.5)

# Requested balanced panel: (class, ordered candidate genes). First assignment wins.
CLASS_DEFS = [
    ("SAGA/chromatin", ["SGF29", "TADA2B", "SUPT20H", "TADA1", "TAF6L", "USP22", "SUPT7L", "ATXN7L3"]),
    ("Mediator", ["MED12", "MED1", "MED24", "CCNC"]),
    ("TCR (context-specific)", ["ZAP70", "LCK", "LAT", "CD3E", "PLCG1"]),
    ("Other robust", ["SENP5", "KDM1A", "NFRKB", "ARNT"]),
    ("Repro-promoted", ["SETDB1", "WDR82", "CPSF6", "WAC", "SMARCE1"]),
    ("Demoted control", ["ELOB", "ELOF1", "SMG1", "EIF4G2"]),
]

# Curated type-I interferon core (for the ISG-specificity test).
ISG_REF = {
    "IFI44", "IFI44L", "IFI6", "IFI35", "IFI27", "OAS1", "OAS2", "OAS3", "OASL", "MX1", "MX2",
    "BST2", "RIGI", "DDX58", "XAF1", "GBP1", "GBP2", "GBP4", "CMPK2", "HELZ2", "ISG15", "IRF7",
    "STAT1", "STAT2", "HERC5", "HERC6", "USP18", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3",
    "RSAD2", "SAMD9", "SAMD9L", "EPSTI1", "PARP9", "PARP12", "DTX3L", "LY6E", "PLSCR1", "APOL1",
    "APOL2", "APOL3", "APOL4", "APOL6", "TRIM22", "SP100", "IFI16",
}


def main():
    panel = pd.read_csv(TAB / "fingerprint_panel.csv")
    var = pd.read_csv(CACHE / "fingerprint_var.csv")
    var_names = var["gene_name"].values
    logfc = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
    zpath = next(CACHE.glob("panel_zscore_*.npy"))
    zscore = np.load(zpath)

    panel_idx = {row["gene"]: i for i, row in panel.iterrows()}
    assigned = {}
    for cls, genes in CLASS_DEFS:
        for g in genes:
            if g in assigned or g not in panel_idx:
                continue
            assigned[g] = cls

    sel = panel[panel["gene"].isin(assigned)].copy()
    sel["reg_class"] = sel["gene"].map(assigned)
    sel = sel.rename(columns={"row": "logfc_row"})
    sel["zrow"] = sel.index
    sel = sel[["gene", "condition", "logfc_row", "zrow", "reg_class",
               "regpower_eb_mean", "n_downstream"]].reset_index(drop=True)
    sel.to_csv(TAB / "balanced_panel_30.csv", index=False)
    print(f"balanced panel: {len(sel)} regulators across {sel.reg_class.nunique()} classes")

    # Recover (log_fc, lfcSE) per selected row, all local.
    recs = []
    for _, r in sel.iterrows():
        y = np.asarray(logfc[int(r["logfc_row"])], dtype=float)
        zs = zscore[int(r["zrow"])].astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            se = np.where(zs != 0, y / zs, np.nan)
        recs.append((r["gene"], r["condition"], r["reg_class"], y, se))

    allY = np.concatenate([r[3] for r in recs])
    allSE = np.concatenate([r[4] for r in recs])
    ok = np.isfinite(allY) & np.isfinite(allSE) & (allSE > 0)
    tau2 = max(float(np.var(allY[ok]) - np.median(allSE[ok] ** 2)), 1e-6)

    rows = []
    for g, c, cls, y, se in recs:
        valid = np.isfinite(y) & np.isfinite(se) & (se > 0)
        v = 1.0 / (1.0 / tau2 + 1.0 / se ** 2)
        m = v * (y / se ** 2)
        sd = np.sqrt(v)
        p_big = norm.sf(LOG2_1P5, loc=m, scale=sd) + norm.cdf(-LOG2_1P5, loc=m, scale=sd)
        p_pos = norm.sf(0, loc=m, scale=sd)
        keep = valid & (p_big > 0.8)
        for j in np.where(keep)[0]:
            rows.append({"perturbed_gene": g, "condition": c, "reg_class": cls,
                         "measured_gene": var_names[j], "log_fc": round(float(y[j]), 3),
                         "lfcSE": round(float(se[j]), 3), "theta_post_mean": round(float(m[j]), 3),
                         "theta_post_sd": round(float(sd[j]), 3),
                         "p_effect_positive": round(float(p_pos[j]), 3),
                         "p_abs_effect_gt_1p5x": round(float(p_big[j]), 3)})
    edges = pd.DataFrame(rows)
    edges.to_csv(TAB / "robust_edges_balanced30.csv", index=False)
    print(f"robust edges: {len(edges):,} across {edges.perturbed_gene.nunique()} regulators")

    # Per-class convergent targets (>= half the class members, self-edges removed).
    all_regs = set(sel["gene"])
    reg_by_class = sel.groupby("reg_class")["gene"].apply(set).to_dict()
    class_targets = {}
    ct_rows = []
    universe = set(var_names)
    isg_u = ISG_REF & universe
    for cls, members in reg_by_class.items():
        sub = edges[(edges.reg_class == cls) & (~edges.measured_gene.isin(all_regs))]
        thr = max(2, int(np.ceil(len(members) / 2)))
        conv = sub.groupby("measured_gene")["perturbed_gene"].nunique()
        tgts = set(conv[conv >= thr].index)
        class_targets[cls] = tgts
        for gtar in sorted(tgts):
            ct_rows.append({"reg_class": cls, "target_gene": gtar, "is_ISG": gtar in isg_u})
    pd.DataFrame(ct_rows).to_csv(TAB / "class_convergent_targets.csv", index=False)

    # ISG specificity (hypergeometric vs measured-gene background).
    N, K = len(universe), len(isg_u)
    iso_rows = []
    for cls, tgts in class_targets.items():
        t = tgts & universe
        n = len(t)
        k = len(t & isg_u)
        p = hypergeom.sf(k - 1, N, K, n) if n else 1.0
        fold = (k / n) / (K / N) if n else 0.0
        subisg = edges[(edges.reg_class == cls) & (edges.measured_gene.isin(isg_u))]
        up = float((subisg.theta_post_mean > 0).mean()) if len(subisg) else None
        iso_rows.append({"class": cls, "targets": n, "ISGs": k, "fold": round(fold, 1),
                         "p": p, "isg_edges": len(subisg),
                         "frac_up_on_KD": (round(up, 3) if up is not None else None)})
    iso = pd.DataFrame(iso_rows)
    iso.to_csv(TAB / "class_isg_enrichment.csv", index=False)

    # Pairwise Jaccard of class target sets.
    classes = list(class_targets)
    jac = {a: {b: (len(class_targets[a] & class_targets[b]) /
                   len(class_targets[a] | class_targets[b])
                   if (class_targets[a] | class_targets[b]) else 0.0)
               for b in classes} for a in classes}
    offdiag = [jac[a][b] for a in classes for b in classes if a != b]

    summary = {cls: {"n_regulators": int((sel.reg_class == cls).sum()),
                     "n_convergent_targets": int(len(class_targets[cls])),
                     "isg_fold": float(iso.set_index("class").loc[cls, "fold"]),
                     "isg_p": float(iso.set_index("class").loc[cls, "p"])}
               for cls in classes}
    summary["_jaccard_offdiag_median"] = float(np.median(offdiag))
    summary["_n_edges_total"] = int(len(edges))
    json.dump(summary, open(TAB / "class_program_summary.json", "w"), indent=2)
    print(f"jaccard off-diagonal median = {np.median(offdiag):.3f}")
    print(iso.to_string(index=False))

    _draw_figure(edges, class_targets, jac, iso, sel, classes)
    print(f"\ndone. figure -> {FIG / '27_regulator_class_programs.png'}")


def _draw_figure(edges, class_targets, jac, iso, sel, classes):
    """Three-panel summary. Self-contained styling (no external helpers)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = ["SAGA/chromatin", "Mediator", "TCR (context-specific)",
             "Other robust", "Repro-promoted", "Demoted control"]
    order = [c for c in order if c in classes]
    col = {"SAGA/chromatin": "#8e44ad", "Mediator": "#2980b9",
           "TCR (context-specific)": "#16a085", "Other robust": "#7f8c8d",
           "Repro-promoted": "#f39c12", "Demoted control": "#bdbdbd"}
    short = {"SAGA/chromatin": "SAGA", "Mediator": "Mediator",
             "TCR (context-specific)": "TCR", "Other robust": "Other",
             "Repro-promoted": "Promoted", "Demoted control": "Demoted"}
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 6, "ytick.labelsize": 6})

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))

    # A: edges per regulator, colored by class
    epr = edges.groupby(["reg_class", "perturbed_gene"]).size().rename("n").reset_index()
    epr["reg_class"] = pd.Categorical(epr["reg_class"], order, ordered=True)
    epr = epr.sort_values(["reg_class", "n"], ascending=[True, False]).reset_index(drop=True)
    yp = np.arange(len(epr))[::-1]
    axes[0].barh(yp, epr["n"], color=[col[c] for c in epr["reg_class"]], height=0.72)
    axes[0].set_yticks(yp); axes[0].set_yticklabels(epr["perturbed_gene"], fontsize=5.2)
    axes[0].set_xlabel("Robust edges (P > 0.8)")
    axes[0].set_title("All 30 regulators carry edges", loc="left")
    axes[0].set_ylim(-0.7, len(epr) - 0.3)
    from matplotlib.patches import Patch
    axes[0].legend(handles=[Patch(color=col[c], label=short[c]) for c in order],
                   fontsize=5, loc="lower right", frameon=False)

    # B: Jaccard heatmap
    J = np.array([[jac[a][b] for b in order] for a in order])
    im = axes[1].imshow(J, cmap="PuBu", vmin=0, vmax=0.3)
    labs = [short[c] for c in order]
    axes[1].set_xticks(range(len(order))); axes[1].set_xticklabels(labs, rotation=45, ha="right", fontsize=6)
    axes[1].set_yticks(range(len(order))); axes[1].set_yticklabels(labs, fontsize=6)
    for i in range(len(order)):
        for j in range(len(order)):
            v = J[i, j]
            axes[1].text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5,
                         color="white" if v > 0.15 else "0.25")
    axes[1].set_title("Classes converge on distinct programs", loc="left")
    cb = fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04)
    cb.set_label("Jaccard of target sets", fontsize=6); cb.ax.tick_params(labelsize=5.5)

    # C: ISG fold enrichment per class
    iso_o = iso.set_index("class").loc[order].reset_index()
    yp = np.arange(len(iso_o))[::-1]
    axes[2].barh(yp, iso_o["fold"], color=[col[c] for c in iso_o["class"]], height=0.66)
    axes[2].set_yticks(yp); axes[2].set_yticklabels([short[c] for c in iso_o["class"]], fontsize=6)
    axes[2].set_xlabel("Interferon-gene fold enrichment")
    axes[2].set_title("Interferon program is most concentrated under SAGA", loc="left")
    for y, f, p, k, n in zip(yp, iso_o["fold"], iso_o["p"], iso_o["ISGs"], iso_o["targets"]):
        sig = "***" if p < 1e-10 else ("**" if p < 1e-4 else ("*" if p < 0.05 else "ns"))
        axes[2].text(f + 0.4, y, f"{k}/{n} {sig}", va="center", fontsize=5.2, color="0.3")
    axes[2].set_xlim(0, iso_o["fold"].max() * 1.30)
    axes[2].axvline(1, color="0.6", lw=0.7, ls=":")
    for ax in axes[::2]:
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)

    fig.tight_layout(w_pad=2.2)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "27_regulator_class_programs.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
