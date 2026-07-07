#!/usr/bin/env python3
"""Regulator ranking audit — valida y endurece el ranking de reguladores.

Sin descargas, sin PPLs, sin MCMC, sin dependencias nuevas. Solo hace el ranking
más defendible con: (1) comparación de baselines, (2) estabilidad por bootstrap,
(3) reguladores globales vs condition-specific.

    python scripts/audit_ranking.py

Outputs:
    docs/tables/ranking_baseline_comparison.csv
    docs/tables/hub_ranking_stability.csv
    docs/tables/top_global_regulators.csv
    docs/tables/top_condition_specific_regulators.csv
    docs/figures/08_kd_gate_changes_ranking.png
    docs/figures/10_global_vs_context_specific.png
    docs/figures/11_ranking_stability.png
    (+ añade stability_frequency y regulator_class a top_regulators_for_review.csv)
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from model_hubs import fit_fixed_effects   # reutiliza el GLM de efectos fijos

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "suppl_tables"
FIG = ROOT / "docs" / "figures"
TAB = ROOT / "docs" / "tables"
ACCENT, VIOLET, UP, MUT, AMBER = "#0a8f9c", "#6b5fc0", "#d6412f", "#586074", "#b47818"
B = 200
TOP_N = 30
SPEC_THRESHOLD = 0.6
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": "#cfd6e0", "font.size": 10, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": "#eceff4"})

SIGMA2E = None   # varianza de ruido global (se fija con todas las filas)


def eb_scores(sub):
    """Media posterior EB del efecto por gen (u_g) para un conjunto de filas."""
    g = sub.groupby("gene")
    d = g["work"].mean()
    n = g.size()
    s2 = SIGMA2E / n
    tau2 = max(float(d.var() - s2.mean()), 1e-6)
    shrink = tau2 / (tau2 + s2)
    return shrink * d


def main():
    global SIGMA2E
    de = pd.read_csv(DATA / "DE_stats.suppl_table.csv")
    de["ontarget_significant"] = de["ontarget_significant"].astype(bool)
    de = de.rename(columns={"target_contrast_gene_name": "gene"})

    print("== Regulator ranking audit ==")
    mu, alpha, _ = fit_fixed_effects(de.rename(columns={"gene": "target_contrast_gene_name"}))
    de["work"] = np.log(de["n_downstream"] + 0.5) - np.log(mu + 0.5)
    de["logmu"] = np.log(mu + 0.5)
    SIGMA2E = float(np.var(de["work"]))
    sig = de[de["ontarget_significant"]].copy()          # filas elegibles (KD gate)
    print(f"  filas totales={len(de):,} · elegibles (KD sig)={len(sig):,} · sigma2_e={SIGMA2E:.3f}")

    # ---------- 1. BASELINE COMPARISON ----------
    raw = de.groupby("gene")["n_downstream"].max()                     # ignora KD
    kdg = sig.groupby("gene")["n_downstream"].max()                    # peak entre sig
    eb = eb_scores(sig)                                                # EB con gate
    comp = pd.DataFrame({"raw_score": raw})
    comp["kd_gated_score"] = kdg
    comp["eb_score"] = eb.round(4)
    comp["raw_rank"] = comp["raw_score"].rank(ascending=False, method="min").astype("Int64")
    comp["kd_gated_rank"] = comp["kd_gated_score"].rank(ascending=False, method="min").astype("Int64")
    comp["eb_rank"] = comp["eb_score"].rank(ascending=False, method="min").astype("Int64")
    comp["dropped_by_kd_gate"] = comp["kd_gated_score"].isna()          # hub crudo sin KD válido
    comp["rank_shift_raw_to_eb"] = (comp["raw_rank"] - comp["eb_rank"]).astype("Float64")
    comp = comp.sort_values("raw_rank")
    comp.reset_index().rename(columns={"index": "gene"}).to_csv(
        TAB / "ranking_baseline_comparison.csv", index=False)
    n_drop = int(comp.head(TOP_N)["dropped_by_kd_gate"].sum())
    print(f"  [1] baseline comparison → ranking_baseline_comparison.csv "
          f"({n_drop}/{TOP_N} hubs crudos caen por el gate de KD)")

    # figura 08 — slopegraph de 3 columnas: raw -> KD-gated -> EB
    top_raw = comp.sort_values("raw_rank").head(20)
    CAP, SINK = 40, 45
    def ypos(rank, dropped=False):
        if dropped or pd.isna(rank):
            return SINK
        return min(int(rank), CAP)
    fig, ax = plt.subplots(figsize=(8, 6.8))
    for gene, r in top_raw.iterrows():
        y0 = ypos(r["raw_rank"])
        y1 = ypos(r["kd_gated_rank"], r["dropped_by_kd_gate"])
        y2 = ypos(r["eb_rank"], r["dropped_by_kd_gate"])
        if r["dropped_by_kd_gate"]:
            col = UP                                    # cae por el gate de KD
        elif int(r["eb_rank"]) <= TOP_N:
            col = ACCENT                                # sobrevive robusto
        else:
            col = AMBER                                 # demotado por shrinkage (condition-specific)
        ax.plot([0, 1, 2], [y0, y1, y2], "-o", color=col, lw=1.6, ms=4.5, alpha=.85)
        ax.text(-0.05, y0, gene, ha="right", va="center", fontsize=8, color="#333")
        if col is ACCENT:
            ax.text(2.06, y2, f"#{int(r['eb_rank'])}", ha="left", va="center", fontsize=7.5, color=col)
        elif r["dropped_by_kd_gate"]:
            ax.text(2.06, y2, "sin KD válido", ha="left", va="center", fontsize=7.5, color=col)
    ax.axhspan(SINK - 1.3, SINK + 1.3, color=UP, alpha=.06)
    ax.text(2.06, SINK, "" , fontsize=7)
    ax.text(1.0, SINK + 2.6, "amarillo = hub condition-specific demotado por el shrinkage EB",
            ha="center", fontsize=8, color=AMBER, style="italic")
    ax.set_xlim(-0.55, 2.7); ax.set_ylim(SINK + 4, 0)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(["raw\n(n_downstream)", "KD-gated\n(solo KD sig)", "EB\n(gate + shrinkage)"])
    ax.set_yticks([1, 10, 20, 30, 40]); ax.set_ylabel("posición en el ranking (top-40; ▼ = fuera)")
    ax.set_title("Por qué importan los gates de calidad", fontweight="bold", fontsize=12)
    ax.grid(False)
    fig.tight_layout(); fig.savefig(FIG / "08_kd_gate_changes_ranking.png", dpi=135, bbox_inches="tight")
    plt.close(fig)
    print("      figura → docs/figures/08_kd_gate_changes_ranking.png")

    # ---------- 2. RANKING STABILITY (bootstrap) ----------
    rng = np.random.default_rng(0)
    idx = np.arange(len(sig))
    sig_small = sig[["gene", "work"]].reset_index(drop=True)
    ranks = {}
    for b in range(B):
        samp = sig_small.iloc[rng.choice(idx, len(idx), replace=True)]
        r = eb_scores(samp).rank(ascending=False, method="first")
        for gene, rk in r.items():
            ranks.setdefault(gene, []).append(rk)
    rows = []
    for gene, rl in ranks.items():
        arr = np.array(rl)
        rows.append({
            "gene": gene,
            "selection_frequency_top30": round(float(np.mean(arr <= TOP_N)), 3),
            "median_rank": float(np.median(arr)),
            "rank_iqr": round(float(np.percentile(arr, 75) - np.percentile(arr, 25)), 1),
            "n_appearances": len(arr),
        })
    stab = pd.DataFrame(rows).sort_values(
        ["selection_frequency_top30", "median_rank"], ascending=[False, True])
    stab.to_csv(TAB / "hub_ranking_stability.csv", index=False)
    print(f"  [2] stability (bootstrap B={B}) → hub_ranking_stability.csv")

    # figura 11 — estabilidad de los top-30 por EB puntual
    eb_top = eb.sort_values(ascending=False).head(TOP_N).index
    s = stab.set_index("gene").reindex(eb_top)
    fig, ax = plt.subplots(figsize=(7, 7))
    y = np.arange(len(s))[::-1]
    freq = s["selection_frequency_top30"].fillna(0).values
    colors = [ACCENT if f >= 0.8 else (AMBER if f >= 0.5 else UP) for f in freq]
    ax.barh(y, freq, color=colors, edgecolor="white")
    ax.set_yticks(y); ax.set_yticklabels(s.index, fontsize=8)
    ax.set_xlim(0, 1.02); ax.axvline(0.8, color=MUT, ls="--", lw=1)
    ax.set_xlabel("frecuencia de selección en el top-30  (bootstrap B=200)")
    ax.set_title("Estabilidad del ranking — ¿el gen se queda en el top?", fontweight="bold", fontsize=11)
    ax.grid(axis="y")
    fig.tight_layout(); fig.savefig(FIG / "11_ranking_stability.png", dpi=135, bbox_inches="tight")
    plt.close(fig)
    print("      figura → docs/figures/11_ranking_stability.png")

    # ---------- 3. GLOBAL vs CONDITION-SPECIFIC ----------
    gsig = sig.groupby("gene")
    n_sig_cond = gsig["culture_condition"].nunique()
    dn_sum = gsig["n_downstream"].sum()
    dn_max = gsig["n_downstream"].max()
    cond_spec = (dn_max / dn_sum).round(3)
    peak_cond = (sig.sort_values("n_downstream", ascending=False)
                    .drop_duplicates("gene").set_index("gene")["culture_condition"])
    cls = pd.DataFrame({
        "regpower_eb_mean": eb.round(4),
        "n_signif_conditions": n_sig_cond,
        "condition_specificity": cond_spec,
        "peak_condition": peak_cond,
        "observed_n_downstream": dn_max,
    }).dropna(subset=["regpower_eb_mean"])
    cls["regulator_class"] = np.where(
        (cls["n_signif_conditions"] >= 2) & (cls["condition_specificity"] < SPEC_THRESHOLD),
        "global", "condition-specific")

    g_tab = (cls[cls.regulator_class == "global"].sort_values("regpower_eb_mean", ascending=False)
             .head(20).reset_index().rename(columns={"index": "gene"}))
    c_tab = (cls[cls.regulator_class == "condition-specific"].sort_values("regpower_eb_mean", ascending=False)
             .head(20).reset_index().rename(columns={"index": "gene"}))
    g_tab.to_csv(TAB / "top_global_regulators.csv", index=False)
    c_tab.to_csv(TAB / "top_condition_specific_regulators.csv", index=False)
    n_glob = int((cls.regulator_class == "global").sum())
    print(f"  [3] global vs condition-specific → 2 tablas "
          f"({n_glob} global / {len(cls)-n_glob} condition-specific)")

    # figura 10 — scatter especificidad vs poder regulatorio
    fig, ax = plt.subplots(figsize=(7.2, 5))
    for klass, col in [("global", ACCENT), ("condition-specific", VIOLET)]:
        d = cls[cls.regulator_class == klass]
        ax.scatter(d["condition_specificity"], d["regpower_eb_mean"], s=14, alpha=.35,
                   color=col, edgecolors="none", label=klass)
    lab = pd.concat([g_tab.head(5).set_index("gene"), c_tab.head(5).set_index("gene")])
    for k, (gene, r) in enumerate(lab.iterrows()):
        ax.annotate(gene, (r["condition_specificity"], r["regpower_eb_mean"]),
                    fontsize=7.5, color="#333", xytext=(4, 5 + (k % 3) * 8),
                    textcoords="offset points",
                    arrowprops=dict(arrowstyle="-", color="#bbb", lw=.5))
    ax.axvline(SPEC_THRESHOLD, color=MUT, ls="--", lw=1)
    ax.set_xlabel("condition_specificity = max / sum de n_downstream (KD sig)")
    ax.set_ylabel("poder regulatorio EB (u_g)")
    ax.set_title("Reguladores globales vs context-specific", fontweight="bold", fontsize=11)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout(); fig.savefig(FIG / "10_global_vs_context_specific.png", dpi=135, bbox_inches="tight")
    plt.close(fig)
    print("      figura → docs/figures/10_global_vs_context_specific.png")

    # ---------- 4. enriquecer top_regulators_for_review.csv ----------
    rev_path = TAB / "top_regulators_for_review.csv"
    if rev_path.exists():
        rev = pd.read_csv(rev_path)
        rev = rev.merge(stab[["gene", "selection_frequency_top30"]], on="gene", how="left")
        rev = rev.rename(columns={"selection_frequency_top30": "stability_frequency"})
        rev = rev.merge(cls[["regulator_class"]].reset_index().rename(columns={"index": "gene"}),
                        on="gene", how="left")
        rev.to_csv(rev_path, index=False)
        print("  [4] top_regulators_for_review.csv += stability_frequency, regulator_class")

    print("\n✓ Audit completo.")


if __name__ == "__main__":
    main()
