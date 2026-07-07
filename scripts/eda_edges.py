#!/usr/bin/env python3
"""EDA del análisis de edges (Modelo 1) — ¿resultado fuerte o bonus?

Usa docs/tables/robust_edges.csv si existe (generado por scripts/model_edges.py).
No descarga nada. Si no existe, informa y termina.

    python scripts/eda_edges.py

Outputs:
    docs/tables/edge_summary_by_regulator.csv
    docs/tables/top_downstream_genes.csv
    docs/figures/14_edges_per_regulator.png
    docs/figures/15_top_downstream_genes.png
    docs/figures/16_edge_direction_balance.png
    docs/figures/17_edges_by_condition.png
    docs/figures/18_edge_bipartite_network.png
    docs/EDGE_ANALYSIS.md
"""
from pathlib import Path
import sys
import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
DOCS = ROOT / "docs"
ACCENT, UP, DOWN, MUT, AMBER, VIOLET = "#0a8f9c", "#d6412f", "#2f5ed0", "#586074", "#b47818", "#6b5fc0"
plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white",
    "axes.edgecolor": "#cfd6e0", "font.size": 10, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": "#eceff4"})
SEP = "; "


def joinvals(s, n=5):
    return SEP.join(s[:n])


def main():
    edges_path = TAB / "robust_edges.csv"
    if not edges_path.exists():
        print("docs/tables/robust_edges.csv no existe — corre antes (opcional): "
              "python scripts/model_edges.py")
        sys.exit(0)
    e = pd.read_csv(edges_path)
    e["abs_theta"] = e["theta_post_mean"].abs()
    e["sign"] = np.where(e["theta_post_mean"] > 0, "pos", "neg")
    n_reg = e["perturbed_gene"].nunique()
    n_cond = e["condition"].nunique()
    print(f"== EDA edges · {len(e):,} edges · {n_reg} reguladores · {n_cond} condiciones ==")

    # ---------- 1. resumen por regulador ----------
    def reg_summary(g):
        gs = g.sort_values("abs_theta", ascending=False)
        return pd.Series({
            "n_edges": len(g),
            "n_positive": int((g["sign"] == "pos").sum()),
            "n_negative": int((g["sign"] == "neg").sum()),
            "mean_abs_theta": round(g["abs_theta"].mean(), 3),
            "median_abs_theta": round(g["abs_theta"].median(), 3),
            "top_downstream_genes": joinvals(list(gs["measured_gene"])),
        })
    by_reg = (e.groupby(["perturbed_gene", "condition"]).apply(reg_summary)
                .reset_index().sort_values("n_edges", ascending=False))
    by_reg.to_csv(TAB / "edge_summary_by_regulator.csv", index=False)
    print("  tabla → docs/tables/edge_summary_by_regulator.csv")

    # ---------- 2. top downstream genes ----------
    def down_summary(g):
        gs = g.sort_values("abs_theta", ascending=False)
        return pd.Series({
            "n_upstream_regulators": g["perturbed_gene"].nunique(),
            "n_positive_edges": int((g["sign"] == "pos").sum()),
            "n_negative_edges": int((g["sign"] == "neg").sum()),
            "mean_abs_theta": round(g["abs_theta"].mean(), 3),
            "conditions_seen": SEP.join(sorted(g["condition"].unique())),
            "top_upstream_regulators": joinvals(list(gs["perturbed_gene"].drop_duplicates())),
        })
    by_down = (e.groupby("measured_gene").apply(down_summary).reset_index()
                 .sort_values(["n_upstream_regulators", "mean_abs_theta"], ascending=False))
    by_down.to_csv(TAB / "top_downstream_genes.csv", index=False)
    print("  tabla → docs/tables/top_downstream_genes.csv")

    # ---------- figuras ----------
    # 14 — edges por regulador (stacked pos/neg)
    r = by_reg.groupby("perturbed_gene")[["n_positive", "n_negative"]].sum().sort_values("n_positive")
    fig, ax = plt.subplots(figsize=(6.8, 4))
    ax.barh(r.index, r["n_positive"], color=UP, label="activación (θ>0)")
    ax.barh(r.index, r["n_negative"], left=r["n_positive"], color=DOWN, label="represión (θ<0)")
    ax.set_xlabel("nº de edges robustos"); ax.set_title("Edges por regulador", fontweight="bold")
    ax.legend(frameon=False, loc="lower right"); ax.grid(axis="y")
    fig.tight_layout(); fig.savefig(FIG / "14_edges_per_regulator.png", dpi=135, bbox_inches="tight"); plt.close(fig)

    # 15 — top downstream por nº de reguladores upstream
    td = by_down.head(15).iloc[::-1]
    fig, ax = plt.subplots(figsize=(6.8, 5))
    ax.barh(td["measured_gene"], td["n_upstream_regulators"], color=ACCENT, edgecolor="white")
    ax.set_xlabel(f"nº de reguladores upstream (de {n_reg})")
    ax.set_title("Downstream más convergentes (targets compartidos)", fontweight="bold", fontsize=11)
    ax.grid(axis="y"); ax.set_xticks(range(0, n_reg + 1))
    fig.tight_layout(); fig.savefig(FIG / "15_top_downstream_genes.png", dpi=135, bbox_inches="tight"); plt.close(fig)

    # 16 — balance de dirección (100% apilado por regulador)
    fig, ax = plt.subplots(figsize=(6.8, 4))
    frac_pos = r["n_positive"] / (r["n_positive"] + r["n_negative"])
    ax.barh(r.index, frac_pos, color=UP, label="activación")
    ax.barh(r.index, 1 - frac_pos, left=frac_pos, color=DOWN, label="represión")
    ax.axvline(0.5, color=MUT, ls="--", lw=1)
    ax.set_xlim(0, 1); ax.set_xlabel("fracción de edges"); ax.grid(axis="y")
    ax.set_title("Balance de dirección — estos reguladores son mayormente activadores",
                 fontweight="bold", fontsize=10.5)
    ax.legend(frameon=False, loc="lower center", ncol=2)
    fig.tight_layout(); fig.savefig(FIG / "16_edge_direction_balance.png", dpi=135, bbox_inches="tight"); plt.close(fig)

    # 17 — edges por condición
    cc = e.groupby(["condition", "sign"]).size().unstack(fill_value=0)
    for c in ("pos", "neg"):
        if c not in cc: cc[c] = 0
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.bar(cc.index, cc["pos"], color=UP, label="activación")
    ax.bar(cc.index, cc["neg"], bottom=cc["pos"], color=DOWN, label="represión")
    ax.set_ylabel("nº de edges"); ax.set_title("Edges por condición", fontweight="bold")
    ax.legend(frameon=False); ax.grid(axis="x")
    fig.text(0.5, -0.02, "Nota: la demo eligió la condición pico por regulador → falta Rest",
             ha="center", fontsize=8, color=MUT, style="italic")
    fig.tight_layout(); fig.savefig(FIG / "17_edges_by_condition.png", dpi=135, bbox_inches="tight"); plt.close(fig)

    # 18 — red bipartita regulador → downstream convergentes
    core = by_down[by_down["n_upstream_regulators"] >= max(2, n_reg - 2)].head(18)["measured_gene"].tolist()
    sub = e[e["measured_gene"].isin(core)]
    regs = list(r.index)                                   # orden por nº de edges
    fig, ax = plt.subplots(figsize=(7.4, 8))
    ry = {g: i for i, g in enumerate(np.linspace(1, len(core), len(regs)))}
    reg_y = {g: v for g, v in zip(regs, np.linspace(1, len(core), len(regs)))}
    down_y = {g: i + 1 for i, g in enumerate(core[::-1])}
    for _, row in sub.iterrows():
        col = UP if row["theta_post_mean"] > 0 else DOWN
        ax.plot([0, 1], [reg_y[row["perturbed_gene"]], down_y[row["measured_gene"]]],
                color=col, lw=0.7, alpha=min(0.9, 0.25 + row["abs_theta"] / 4))
    for g, y in reg_y.items():
        ax.scatter(0, y, s=90, color=VIOLET, zorder=3)
        ax.text(-0.03, y, g, ha="right", va="center", fontsize=9, fontweight="bold")
    for g, y in down_y.items():
        ax.scatter(1, y, s=40, color="#888", zorder=3)
        ax.text(1.03, y, g, ha="left", va="center", fontsize=8)
    ax.set_xlim(-0.35, 1.4); ax.set_ylim(0, len(core) + 1); ax.axis("off")
    ax.text(0, len(core) + 0.8, "reguladores", ha="center", fontsize=9, color=VIOLET, fontweight="bold")
    ax.text(1, len(core) + 0.8, "downstream convergentes", ha="center", fontsize=9, color="#555", fontweight="bold")
    ax.set_title("Red bipartita — targets compartidos por el complejo SAGA",
                 fontweight="bold", fontsize=11)
    ax.text(0.5, -0.4, "rojo = activación · azul = represión", ha="center", fontsize=8, color=MUT)
    fig.tight_layout(); fig.savefig(FIG / "18_edge_bipartite_network.png", dpi=135, bbox_inches="tight"); plt.close(fig)
    print("  figuras → docs/figures/14–18_*.png")

    # ---------- EDGE_ANALYSIS.md ----------
    shared2 = int((e.groupby("measured_gene")["perturbed_gene"].nunique() >= 2).sum())
    pos_frac = (e["sign"] == "pos").mean()
    md = f"""# Análisis de edges (Modelo 1) — ¿resultado fuerte o bonus?

## Qué se analizó

La red de efectos con incertidumbre `docs/tables/robust_edges.csv`: **{len(e):,} edges** robustos
(`P(|efecto|>1.5×)>0.8`) de **{n_reg} reguladores** × **{n_cond} condiciones**, leídos por slice del
`DE_stats.h5ad` remoto. Se resumió por regulador y por gen downstream, con dirección y convergencia.

## ¿Se ven útiles?

**Sí como prueba de concepto, no como resultado fuerte.** Señales a favor:

- **Convergencia coherente**: {shared2} genes downstream son target de ≥2 reguladores. Los {n_reg}
  reguladores son co-miembros del **complejo SAGA** (TADA1/TADA2B/SGF29/SUPT20H) → que compartan
  targets es exactamente lo esperado, y **valida que el método recupera estructura biológica real**.
- **Dirección interpretable**: {pos_frac*100:.0f}% de los edges son de activación, consistente con
  SAGA como coactivador de la transcripción.
- La magnitud de los efectos es modesta y bien acotada (|θ| mediana ≈ {e['abs_theta'].median():.2f}).

## ¿Por qué NO es un resultado fuerte (todavía)

- **Cobertura mínima**: {n_reg} de 7,913 reguladores (0.08%), **seleccionados por el ranking EB**
  → muestra sesgada a un solo complejo, no representa el paisaje regulatorio.
- **Condiciones incompletas**: solo {', '.join(sorted(e['condition'].unique()))} (la demo tomó la
  condición pico por regulador; **falta Rest**).
- **Latencia remota**: escalar a ~150 reguladores son ~11 min de lectura del h5ad (4.5 s/fila medidos).
- **Semántica de la probabilidad**: `p_abs_effect_gt_1p5x` es **P(magnitud del efecto > 1.5×)**,
  NO P(existe una arista causal). No hay control de FDR a nivel de red ni de sparsity (spike-and-slab).

## Veredicto

Dejar como **bonus / proof-of-concept**: demuestra que el pipeline produce una red de efectos con
incertidumbre biológicamente sensata sin descargar 1.8 TB, pero la cobertura y la semántica no
sostienen todavía afirmaciones de red a escala. Para promoverlo a resultado fuerte: correr los ~150 top reguladores en
las 3 condiciones y añadir P(edge) con prior sparse.
"""
    (DOCS / "EDGE_ANALYSIS.md").write_text(md)
    print("  doc → docs/EDGE_ANALYSIS.md")
    print(f"\n  VEREDICTO: bonus / proof-of-concept "
          f"({n_reg} reguladores, {n_cond} cond, {shared2} targets compartidos)")
    print("✓ EDA de edges completo.")


if __name__ == "__main__":
    main()
