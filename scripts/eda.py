#!/usr/bin/env python3
"""EDA 80/20 del CD4+ T cell Perturb-seq — solo con las tablas suplementarias.

Unidad de análisis: 1 fila = gen perturbado × condición de cultivo.
No carga ningún .h5ad (1.8 TB). Imprime hallazgos, guarda figuras en docs/figures/
y una tabla accionable en docs/tables/.

    python scripts/eda.py

NOTA sobre reproducibilidad: las métricas cross-guide y cross-donante
(single_guide_estimate, guide_correlation_all, donor_correlation_hits_mean)
NO están en el CSV suplementario; viven en GWCD4i.DE_stats.h5ad (.obs). Aquí
usamos un proxy de reproducibilidad CROSS-CONDICIÓN, que sí es computable.
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
    print("  figura →", (FIG / name).relative_to(ROOT))

# ----------------------------------------------------------------------
de = pd.read_csv(DATA / "DE_stats.suppl_table.csv")
sg = pd.read_csv(DATA / "sgrna_library_metadata.suppl_table.csv")
sm = pd.read_csv(DATA / "sample_metadata.suppl_table.csv")
de["sig"] = de["ontarget_significant"].astype(bool)
de["kd_strength"] = -de["ontarget_effect_size"]   # más positivo = knockdown más fuerte

h("SCOPE — qué se usó y qué NO")
print("Usado (tablas suplementarias, ~15 MB):")
print(f"  DE_stats     : {de.shape[0]:>7,} filas (gen perturbado × condición) × {de.shape[1]-2} cols")
print(f"  sgRNA library: {sg.shape[0]:>7,} guías")
print(f"  samples      : {sm.shape[0]:>7,} muestras")
print("NO cargado: los .h5ad/.h5mu (1.8 TB) — incluidas las métricas de reproducibilidad.")

# ----------------------------------------------------------------------
h("1 · LOS EFECTOS SON HEAVY-TAILED (usa percentiles, no media)")
n = de["n_total_de_genes"]
print(f"mediana de DEGs por perturbación : {n.median():.0f}   (media engañosa: {n.mean():.1f})")
print(f"sin efecto (0 DEG)               : {(n==0).mean()*100:.1f}%")
print(f"percentil 90 / 99                : {n.quantile(.9):.0f} / {n.quantile(.99):.0f}")
print(f"hubs con >1000 DEGs              : {(n>1000).sum():,} ({(n>1000).mean()*100:.1f}%)")

# ----------------------------------------------------------------------
h("2 · EL KNOCKDOWN FILTRA LA MAYOR PARTE DE LA SEÑAL UTILIZABLE")
print(f"contrastes con KD on-target significativo: {de['sig'].sum():,} / {len(de):,} ({de['sig'].mean()*100:.0f}%)")
print(f"offtarget_flag (posible off-target)      : {de['offtarget_flag'].sum():,} ({de['offtarget_flag'].mean()*100:.1f}%)")
med_sig = de.loc[de.sig, "n_downstream"].median()
med_non = de.loc[~de.sig, "n_downstream"].median()
print(f"mediana n_downstream — sig: {med_sig:.0f}  vs  no-sig: {med_non:.0f}")
frac = de.loc[de.sig, "n_downstream"].sum() / de["n_downstream"].sum()
print(f"→ los contrastes significativos concentran el {frac*100:.0f}% de todos los trans-efectos")

# ----------------------------------------------------------------------
h("3 · LA ESTIMULACIÓN AMPLIA LOS EFECTOS")
by = de.groupby("culture_condition")["n_total_de_genes"].agg(["mean", "median"]).reindex(COND_ORDER)
print(by.round(1).to_string())

# ----------------------------------------------------------------------
h("4 · TOP HUBS — plausibles: señalización TCR e inmunoregulación")
top = de.sort_values("n_downstream", ascending=False).head(12)
print(top[["target_contrast_gene_name", "culture_condition", "n_downstream", "n_up_genes", "n_down_genes"]].to_string(index=False))

# ----------------------------------------------------------------------
h("5 · LIBRERÍA — cobertura ~2 guías por gen")
gpg = sg.groupby("target_gene_name").size()
print(f"genes diana únicos: {sg['target_gene_name'].nunique():,} · guías/gen (mediana): {gpg.median():.0f}")
print(gpg.value_counts().sort_index().head(4).to_string())

# ======================================================================
# RANKING DE REGULADORES ROBUSTOS (accionable)
# ======================================================================
h("6 · REGULADORES ROBUSTOS (ranking accionable)")
# Reproducibilidad cross-condición (proxy disponible sin el h5ad):
#   para cada gen tested en varias condiciones, cuán estable es n_downstream.
g = de.groupby("target_contrast_gene_name")
cond_counts = g["culture_condition"].nunique()
sig_counts = g["sig"].sum()
dn_min = g["n_downstream"].min()
dn_max = g["n_downstream"].max().replace(0, np.nan)
repro = (dn_min / dn_max).fillna(0)          # 1 = igual en todas las condiciones
gene_lvl = pd.DataFrame({
    "n_conditions": cond_counts,
    "n_signif_conditions": sig_counts.astype(int),
    "xcond_reproducibility": repro.round(3),
})
de = de.merge(gene_lvl, left_on="target_contrast_gene_name", right_index=True, how="left")

# robust_score SOLO con columnas disponibles (documentado):
#   log1p(downstream) · gate KD signif · penalización off-target ·
#   term de reproducibilidad cross-condición · fracción de condiciones significativas
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
print("\n  tabla →", out.relative_to(ROOT), f"(top 30)")
print("  NOTA: robust_score usa solo columnas del CSV. Para un ranking definitivo,")
print("  incorpora single_guide_estimate + donor_correlation_hits_mean desde DE_stats.h5ad.")

# ======================================================================
# FIGURAS
# ======================================================================
h("GENERANDO FIGURAS")

# 01 — distribución heavy-tailed
fig, ax = plt.subplots(figsize=(7, 4))
labels = ["0", "1–10", "11–50", "51–200", "201–1000", ">1000"]
bins = [0, 1, 11, 51, 201, 1001, n.max() + 1]
counts = pd.cut(n, bins=bins, right=False, labels=labels).value_counts().reindex(labels)
ax.bar(labels, counts.values, color=ACCENT, edgecolor="white")
for i, v in enumerate(counts.values):
    ax.text(i, v + 300, f"{v:,}", ha="center", fontsize=9, color="#444")
ax.set_title("Distribución heavy-tailed de DEGs (resumir con percentiles, no media)", fontweight="bold", fontsize=11)
ax.set_ylabel("nº de contrastes"); ax.set_xlabel("genes diferencialmente expresados (DEGs)")
save(fig, "01_distribution_n_total_de_genes.png")

# 02 — DEGs por condición
fig, ax = plt.subplots(figsize=(6.5, 4))
data = [de.loc[de.culture_condition == c, "n_total_de_genes"].clip(upper=500) for c in COND_ORDER]
parts = ax.violinplot(data, showmedians=True, widths=0.8)
for pc in parts["bodies"]:
    pc.set_facecolor(ACCENT); pc.set_alpha(0.55); pc.set_edgecolor(ACCENT)
for k in ("cmedians", "cmaxes", "cmins", "cbars"): parts[k].set_color(MUT)
ax.set_xticks([1, 2, 3]); ax.set_xticklabels(COND_ORDER)
ax.set_title("Las células estimuladas muestran efectos más amplios", fontweight="bold", fontsize=11)
ax.set_ylabel("DEGs por perturbación (cap 500)")
save(fig, "02_degs_by_condition.png")

# 03 — top hubs por condición
fig, ax = plt.subplots(figsize=(7, 5))
t = de.sort_values("n_downstream", ascending=False).head(15).iloc[::-1]
lab = t["target_contrast_gene_name"] + "  (" + t["culture_condition"].str.replace("Stim", "") + ")"
ax.barh(lab, t["n_up_genes"], color=UP, label="↑ up")
ax.barh(lab, t["n_down_genes"], left=t["n_up_genes"], color=DOWN, label="↓ down")
ax.set_title("Top hubs — enriquecidos en señalización de células T", fontweight="bold", fontsize=11)
ax.set_xlabel("genes trans-afectados"); ax.legend(loc="lower right", frameon=False); ax.grid(axis="y")
save(fig, "03_top_hubs_by_condition.png")

# 04 — KD strength vs trans-efectos (kd_strength: más a la derecha = más KD)
fig, ax = plt.subplots(figsize=(6.5, 4.5))
s = de.dropna(subset=["kd_strength", "n_downstream"]).sample(min(6000, len(de)), random_state=1)
ax.scatter(s["kd_strength"], s["n_downstream"] + 1,
           c=s["sig"].map({True: ACCENT, False: "#c9d0da"}), s=8, alpha=0.5, edgecolors="none")
ax.set_yscale("log")
ax.set_xlabel("kd_strength = −ontarget_effect_size   (más a la derecha → knockdown más fuerte)")
ax.set_ylabel("trans-efectos (n_downstream + 1, log)")
ax.set_title("Los contrastes sin KD significativo casi no concentran trans-efectos", fontweight="bold", fontsize=10.5)
ax.text(0.98, 0.06, "teal = KD significativo · gris = no signif.", transform=ax.transAxes,
        color=MUT, fontsize=9, ha="right")
save(fig, "04_ontarget_vs_downstream.png")

# 05 — guías por gen
fig, ax = plt.subplots(figsize=(6, 3.6))
vc = gpg.value_counts().sort_index(); vc = vc[vc.index <= 4]
ax.bar(vc.index.astype(str), vc.values, color=AMBER, edgecolor="white")
for i, (xx, v) in enumerate(vc.items()):
    ax.text(i, v + 150, f"{v:,}", ha="center", fontsize=9, color="#444")
ax.set_title("Cobertura de la librería: ~2 guías por gen", fontweight="bold", fontsize=11)
ax.set_xlabel("guías por gen diana"); ax.set_ylabel("nº de genes")
save(fig, "05_guides_per_gene.png")

# 06 — reproducibilidad (cross-condición) vs magnitud del efecto
fig, ax = plt.subplots(figsize=(6.8, 4.6))
gg = de.drop_duplicates("target_contrast_gene_name")
gg = gg[(gg["n_conditions"] == 3) & (gg["n_signif_conditions"] >= 1)]
peak = de.groupby("target_contrast_gene_name")["n_downstream"].max()
gg = gg.assign(peak=gg["target_contrast_gene_name"].map(peak))
gg = gg[gg["peak"] > 0]
ax.scatter(gg["xcond_reproducibility"], gg["peak"] + 1, s=10, alpha=0.35,
           c=ACCENT, edgecolors="none")
ax.set_yscale("log")
ax.set_xlabel("reproducibilidad cross-condición  (min/max de n_downstream, 1 = estable)")
ax.set_ylabel("efecto máximo (max n_downstream, log)")
ax.set_title("Reproducibilidad vs magnitud: arriba-derecha = reguladores robustos", fontweight="bold", fontsize=10.5)
ax.text(0.02, 0.96, "cross-CONDICIÓN (proxy). El cross-DONANTE requiere DE_stats.h5ad.",
        transform=ax.transAxes, color=MUT, fontsize=8.5, va="top", style="italic")
save(fig, "06_reproducibility_vs_effects.png")

print("\n✓ EDA completo. Figuras en docs/figures/ · tabla en docs/tables/")
