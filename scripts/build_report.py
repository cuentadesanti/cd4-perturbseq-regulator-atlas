#!/usr/bin/env python3
"""Genera el reporte judge-facing consolidado y la figura de overview del pipeline.

    python scripts/build_report.py

Produce:
    docs/figures/00_pipeline_overview.png
    docs/report.md   (consolida DATA_MODEL.md + EDA.md + MODELING.md + figuras + tabla top)
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

    box(1, 10, 20, 14, "CSV suplementarias", AMBER, "DE_stats · sgRNA · samples\n~15 MB · local")
    arrow(21, 29, 17, "features")
    box(29, 10, 18, 14, "EDA 80/20", ACCENT, "distribución, hubs,\nQC · scripts/eda.py")
    arrow(47, 55, 17, "ranking")
    box(55, 10, 20, 14, "Modelo 2 · EB", ACCENT, "reguladores robustos\nscripts/model_hubs.py")
    arrow(75, 83, 17, "candidatos")
    box(83, 4, 16, 26, "Modelo 1 (opcional)", VIOLET, "red probabilística\nstreaming h5ad\n17 GB · sin bajar")

    ax.text(50, 30.5, "Pipeline — de CSV local a reguladores con incertidumbre",
            ha="center", fontsize=12, fontweight="bold", color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "00_pipeline_overview.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("  figura → docs/figures/00_pipeline_overview.png")


def df_to_md(df):
    """Tabla markdown sin depender de `tabulate`."""
    cols = list(df.columns)
    head = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = "\n".join("| " + " | ".join(str(v) for v in row) + " |"
                     for row in df.itertuples(index=False))
    return "\n".join([head, sep, body])


def section(path):
    """Devuelve el cuerpo de un .md sin su primer título H1 (para anidar)."""
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
        n = len(pd.read_csv(ep))
        edges_line = (f"\nSe generó además una **red probabilística de edges** (bonus, Modelo 1) "
                      f"con **{n:,} edges robustos** (`P(|efecto|>1.5×)>0.8`) en "
                      f"`docs/tables/robust_edges.csv`.\n")

    # --- secciones de la auditoría (condicionales a que exista el audit) ---
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
## Naive hubs vs quality-aware regulators

Rankear por `n_downstream` crudo premia hubs que no sobreviven a los controles de calidad.
De los 30 hubs crudos top: **{n_drop} caen por el gate de KD** (sin knockdown on-target validado)
y **{len(demoted)} se demotan por el shrinkage EB** por ser condition-specific (la señal vive en
una sola condición). El ranking EB surface reguladores con efecto grande **y** estable.

![gate](figures/08_kd_gate_changes_ranking.png)

La estabilidad se auditó con bootstrap (B=200) sobre las filas elegibles: la frecuencia con que
cada gen cae en el top-30 (`stability_frequency`) está en `top_regulators_for_review.csv`. El
ranking es moderadamente estable — conviene leerlo como *conjunto* de reguladores robustos, no
como un orden exacto.

![stability](figures/11_ranking_stability.png)

## Global versus context-specific regulators

Separando por `condition_specificity = max/sum de n_downstream` entre condiciones con KD significativo:

- **Globales** (efecto estable en ≥2 condiciones): {glist}… — maquinaria de cromatina/transcripción.
- **Context-specific** (efecto concentrado en una condición): {clist}… — incluye señalización TCR
  (ZAP70, LCK), activa solo bajo estímulo.

Ambas clases son biología real; la distinción evita confundir un regulador universal con uno de contexto.
Tablas: `top_global_regulators.csv`, `top_condition_specific_regulators.csv`.

![globalvs](figures/10_global_vs_context_specific.png)
"""

    md = f"""# Reporte — Genome-scale CD4+ T cell Perturb-seq

*Reporte consolidado para revisión. Reproducible con `make all` (solo CSV locales).*

![pipeline](figures/00_pipeline_overview.png)

## Pregunta

¿Qué genes son **reguladores robustos** de los programas de células T CD4+, separando
señal real de ruido y priorizando por efecto **grande y reproducible**, no por conteos crudos?

## Resumen ejecutivo

- El efecto de las perturbaciones es **heavy-tailed**: mediana 2 DEGs, pero un 1.5% son hubs
  con >1000. Se resume con percentiles y rankings, no con la media.
- El **knockdown efectivo gatea la señal**: los contrastes con KD on-target significativo (62%)
  concentran el **85%** de todos los trans-efectos.
- Un modelo **empirical-Bayes** (pseudo-bayesiano) rankea reguladores por poder regulatorio
  latente con incertidumbre. El top robusto es maquinaria de **cromatina/transcripción**
  (complejo SAGA, Mediador, KDM1A, SETD2) — efecto grande **y** estable entre condiciones.
{edges_line}
## Top reguladores (para revisión)

{tbl}

Tabla completa (30, con todas las columnas): `docs/tables/top_regulators_for_review.csv`.

![ranking](figures/07_hub_posterior_ranking.png)
{audit_md}
## Hallazgos del EDA

![degs](figures/01_distribution_n_total_de_genes.png)
![hubs](figures/03_top_hubs_by_condition.png)
![kd](figures/04_ontarget_vs_downstream.png)
![repro](figures/06_reproducibility_vs_effects.png)

---

## Anexo A — Modelo de datos

{section("DATA_MODEL.md")}

---

## Anexo B — EDA

{section("EDA.md")}

---

## Anexo C — Modelado

{section("MODELING.md")}
"""
    (ROOT / "docs" / "report.md").write_text(md)
    print("  reporte → docs/report.md")


if __name__ == "__main__":
    print("== build_report ==")
    pipeline_overview()
    build_report()
    print("✓ Reporte generado.")
