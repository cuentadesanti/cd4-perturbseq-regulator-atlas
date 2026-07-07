#!/usr/bin/env python3
"""80/20 EDA of the CD4+ T cell Perturb-seq — supplementary tables only.

Unit of analysis: 1 row = perturbed gene × culture condition.
Loads no .h5ad (1.8 TB). Prints findings, saves figures to docs/figures/
and an actionable table to docs/tables/.

    python scripts/eda.py

NOTE on reproducibility: the cross-guide and cross-donor metrics
(single_guide_estimate, guide_correlation_all, donor_correlation_hits_mean)
are NOT in the supplementary CSV; they live in GWCD4i.DE_stats.h5ad (.obs). Here
we use a CROSS-CONDITION reproducibility proxy, which is computable.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "suppl_tables"
FIG = ROOT / "docs" / "figures"
TAB = ROOT / "docs" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

INK, ACCENT, UP, DOWN, AMBER, MUT = "#12161f", "#22b8c4", "#f0645a", "#4f86ea", "#e0a441", "#8b95a8"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]
plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": "#cfd6e0", "axes.labelcolor": "#333", "text.color": "#222",
    "xtick.color": "#555", "ytick.color": "#555", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#eceff4", "grid.linewidth": 1,
})

def h(t): print("\n" + "=" * 70 + f"\n  {t}\n" + "=" * 70)
def save(fig, name):
    fig.tight_layout(); fig.savefig(FIG / name, dpi=130, bbox_inches="tight"); plt.close(fig)
    print("  figure →", (FIG / name).relative_to(ROOT))

# ----------------------------------------------------------------------
de = pd.read_csv(DATA / "DE_stats.suppl_table.csv")
sg = pd.read_csv(DATA / "sgrna_library_metadata.suppl_table.csv")
sm = pd.read_csv(DATA / "sample_metadata.suppl_table.csv")
de["sig"] = de["ontarget_significant"].astype(bool)
de["kd_strength"] = -de["ontarget_effect_size"]   # more positive = stronger knockdown

h("SCOPE — what was used and what was NOT")
print("Used (supplementary tables, ~15 MB):")
print(f"  DE_stats     : {de.shape[0]:>7,} rows (perturbed gene × condition) × {de.shape[1]-2} cols")
print(f"  sgRNA library: {sg.shape[0]:>7,} guides")
print(f"  samples      : {sm.shape[0]:>7,} samples")
print("NOT loaded: the .h5ad/.h5mu (1.8 TB) — including the reproducibility metrics.")

# ----------------------------------------------------------------------
h("1 · EFFECTS ARE HEAVY-TAILED (use percentiles, not the mean)")
n = de["n_total_de_genes"]
print(f"median DEGs per perturbation     : {n.median():.0f}   (misleading mean: {n.mean():.1f})")
print(f"no effect (0 DEG)                : {(n==0).mean()*100:.1f}%")
print(f"percentile 90 / 99               : {n.quantile(.9):.0f} / {n.quantile(.99):.0f}")
print(f"hubs with >1000 DEGs             : {(n>1000).sum():,} ({(n>1000).mean()*100:.1f}%)")

# ----------------------------------------------------------------------
h("2 · KNOCKDOWN GATES MOST OF THE USABLE SIGNAL")
print(f"contrasts with significant on-target KD : {de['sig'].sum():,} / {len(de):,} ({de['sig'].mean()*100:.0f}%)")
print(f"offtarget_flag (possible off-target)    : {de['offtarget_flag'].sum():,} ({de['offtarget_flag'].mean()*100:.1f}%)")
med_sig = de.loc[de.sig, "n_downstream"].median()
med_non = de.loc[~de.sig, "n_downstream"].median()
print(f"median n_downstream — sig: {med_sig:.0f}  vs  non-sig: {med_non:.0f}")
frac = de.loc[de.sig, "n_downstream"].sum() / de["n_downstream"].sum()
print(f"→ significant contrasts concentrate {frac*100:.0f}% of all trans-effects")

# ----------------------------------------------------------------------
h("3 · STIMULATION BROADENS THE EFFECTS")
by = de.groupby("culture_condition")["n_total_de_genes"].agg(["mean", "median"]).reindex(COND_ORDER)
print(by.round(1).to_string())

# ----------------------------------------------------------------------
h("4 · TOP HUBS — plausible: TCR signaling and immunoregulation")
top = de.sort_values("n_downstream", ascending=False).head(12)
print(top[["target_contrast_gene_name", "culture_condition", "n_downstream", "n_up_genes", "n_down_genes"]].to_string(index=False))

# ----------------------------------------------------------------------
h("5 · LIBRARY — ~2 guides per gene coverage")
gpg = sg.groupby("target_gene_name").size()
print(f"unique target genes: {sg['target_gene_name'].nunique():,} · guides/gene (median): {gpg.median():.0f}")
print(gpg.value_counts().sort_index().head(4).to_string())

# ======================================================================
# ROBUST REGULATOR RANKING (actionable)
# ======================================================================
h("6 · ROBUST REGULATORS (actionable ranking)")
# Cross-condition reproducibility (proxy available without the h5ad):
#   for each gene tested in several conditions, how stable n_downstream is.
g = de.groupby("target_contrast_gene_name")
cond_counts = g["culture_condition"].nunique()
sig_counts = g["sig"].sum()
dn_min = g["n_downstream"].min()
dn_max = g["n_downstream"].max().replace(0, np.nan)
repro = (dn_min / dn_max).fillna(0)          # 1 = equal across all conditions
gene_lvl = pd.DataFrame({
    "n_conditions": cond_counts,
    "n_signif_conditions": sig_counts.astype(int),
    "xcond_reproducibility": repro.round(3),
})
de = de.merge(gene_lvl, left_on="target_contrast_gene_name", right_index=True, how="left")

# robust_score using ONLY available columns (documented):
#   log1p(downstream) · significant-KD gate · off-target penalty ·
#   cross-condition reproducibility term · fraction of significant conditions
de["robust_score"] = (
    np.log1p(de["n_downstream"])
    * de["sig"].astype(int)
    * np.where(de["offtarget_flag"], 0.6, 1.0)
    * (0.5 + 0.5 * de["xcond_reproducibility"])
    * (de["n_signif_conditions"] / 3.0)
).round(3)

cols = ["target_contrast_gene_name", "culture_condition", "n_downstream", "n_total_de_genes",
        "n_up_genes", "n_down_genes", "ontarget_effect_size", "ontarget_significant",
        "offtarget_flag", "n_signif_conditions", "xcond_reproducibility", "robust_score"]
top_robust = de.sort_values("robust_score", ascending=False)[cols].head(30)
out = TAB / "top_robust_regulators.csv"
top_robust.to_csv(out, index=False)
print(top_robust.head(15).to_string(index=False))
print("\n  table →", out.relative_to(ROOT), f"(top 30)")
print("  NOTE: robust_score uses only CSV columns. For a definitive ranking,")
print("  add single_guide_estimate + donor_correlation_hits_mean from DE_stats.h5ad.")

# ======================================================================
# FIGURES
# ======================================================================
h("GENERATING FIGURES")

# 01 — heavy-tailed distribution
fig, ax = plt.subplots(figsize=(7, 4))
labels = ["0", "1–10", "11–50", "51–200", "201–1000", ">1000"]
bins = [0, 1, 11, 51, 201, 1001, n.max() + 1]
counts = pd.cut(n, bins=bins, right=False, labels=labels).value_counts().reindex(labels)
ax.bar(labels, counts.values, color=ACCENT, edgecolor="white")
for i, v in enumerate(counts.values):
    ax.text(i, v + 300, f"{v:,}", ha="center", fontsize=9, color="#444")
ax.set_title("Heavy-tailed distribution of DEGs (summarize with percentiles, not the mean)", fontweight="bold", fontsize=11)
ax.set_ylabel("number of contrasts"); ax.set_xlabel("differentially expressed genes (DEGs)")
save(fig, "01_distribution_n_total_de_genes.png")

# 02 — DEGs by condition
fig, ax = plt.subplots(figsize=(6.5, 4))
data = [de.loc[de.culture_condition == c, "n_total_de_genes"].clip(upper=500) for c in COND_ORDER]
parts = ax.violinplot(data, showmedians=True, widths=0.8)
for pc in parts["bodies"]:
    pc.set_facecolor(ACCENT); pc.set_alpha(0.55); pc.set_edgecolor(ACCENT)
for k in ("cmedians", "cmaxes", "cmins", "cbars"): parts[k].set_color(MUT)
ax.set_xticks([1, 2, 3]); ax.set_xticklabels(COND_ORDER)
ax.set_title("Stimulated cells show broader effects", fontweight="bold", fontsize=11)
ax.set_ylabel("DEGs per perturbation (capped at 500)")
save(fig, "02_degs_by_condition.png")

# 03 — top hubs by condition
fig, ax = plt.subplots(figsize=(7, 5))
t = de.sort_values("n_downstream", ascending=False).head(15).iloc[::-1]
lab = t["target_contrast_gene_name"] + "  (" + t["culture_condition"].str.replace("Stim", "") + ")"
ax.barh(lab, t["n_up_genes"], color=UP, label="↑ up")
ax.barh(lab, t["n_down_genes"], left=t["n_up_genes"], color=DOWN, label="↓ down")
ax.set_title("Top hubs — enriched in T cell signaling", fontweight="bold", fontsize=11)
ax.set_xlabel("trans-affected genes"); ax.legend(loc="lower right", frameon=False); ax.grid(axis="y")
save(fig, "03_top_hubs_by_condition.png")

# 04 — KD strength vs trans-effects (kd_strength: further right = more KD)
fig, ax = plt.subplots(figsize=(6.5, 4.5))
s = de.dropna(subset=["kd_strength", "n_downstream"]).sample(min(6000, len(de)), random_state=1)
ax.scatter(s["kd_strength"], s["n_downstream"] + 1,
           c=s["sig"].map({True: ACCENT, False: "#c9d0da"}), s=8, alpha=0.5, edgecolors="none")
ax.set_yscale("log")
ax.set_xlabel("kd_strength = −ontarget_effect_size   (further right → stronger knockdown)")
ax.set_ylabel("trans-effects (n_downstream + 1, log)")
ax.set_title("Contrasts without a significant KD carry almost no trans-effects", fontweight="bold", fontsize=10.5)
ax.text(0.98, 0.06, "teal = significant KD · gray = non-sig.", transform=ax.transAxes,
        color=MUT, fontsize=9, ha="right")
save(fig, "04_ontarget_vs_downstream.png")

# 05 — guides per gene
fig, ax = plt.subplots(figsize=(6, 3.6))
vc = gpg.value_counts().sort_index(); vc = vc[vc.index <= 4]
ax.bar(vc.index.astype(str), vc.values, color=AMBER, edgecolor="white")
for i, (xx, v) in enumerate(vc.items()):
    ax.text(i, v + 150, f"{v:,}", ha="center", fontsize=9, color="#444")
ax.set_title("Library coverage: ~2 guides per gene", fontweight="bold", fontsize=11)
ax.set_xlabel("guides per target gene"); ax.set_ylabel("number of genes")
save(fig, "05_guides_per_gene.png")

# 06 — reproducibility (cross-condition) vs. effect magnitude
fig, ax = plt.subplots(figsize=(6.8, 4.6))
gg = de.drop_duplicates("target_contrast_gene_name")
gg = gg[(gg["n_conditions"] == 3) & (gg["n_signif_conditions"] >= 1)]
peak = de.groupby("target_contrast_gene_name")["n_downstream"].max()
gg = gg.assign(peak=gg["target_contrast_gene_name"].map(peak))
gg = gg[gg["peak"] > 0]
ax.scatter(gg["xcond_reproducibility"], gg["peak"] + 1, s=10, alpha=0.35,
           c=ACCENT, edgecolors="none")
ax.set_yscale("log")
ax.set_xlabel("cross-condition reproducibility  (min/max of n_downstream, 1 = stable)")
ax.set_ylabel("peak effect (max n_downstream, log)")
ax.set_title("Reproducibility vs. magnitude: top-right = robust regulators", fontweight="bold", fontsize=10.5)
ax.text(0.02, 0.96, "cross-CONDITION (proxy). Cross-DONOR requires DE_stats.h5ad.",
        transform=ax.transAxes, color=MUT, fontsize=8.5, va="top", style="italic")
save(fig, "06_reproducibility_vs_effects.png")

print("\n✓ EDA complete. Figures in docs/figures/ · table in docs/tables/")
