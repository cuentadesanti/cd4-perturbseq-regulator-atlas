#!/usr/bin/env python3
"""Modelo 2 — ranking de reguladores por empirical-Bayes (pseudo-bayesiano).

IMPORTANTE: esto NO es un NB jerárquico completo. No hay PPL, ni random effects
formales, ni posterior conjunto muestreado. Es un NB de efectos fijos (statsmodels)
+ shrinkage empirical-Bayes normal-normal del efecto por gen. Se documenta como tal.

Corre solo con DE_stats.suppl_table.csv (local, ~15 MB). Sin descargas.

Salidas:
    docs/tables/hub_ranking_bayes.csv          (todos los genes, ranking completo)
    docs/tables/top_regulators_for_review.csv  (top 30, judge-facing)
    docs/figures/07_hub_posterior_ranking.png

Uso:
    python scripts/model_hubs.py
"""
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "suppl_tables"
FIG = ROOT / "docs" / "figures"
TAB = ROOT / "docs" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

ACCENT, MUT = "#22b8c4", "#8b95a8"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def fit_fixed_effects(de):
    """Media condicional (baseline) por GLM de efectos fijos: condición + calidad de KD.

    Poisson y NB comparten el MISMO modelo de media; como solo usamos mu (no inferencia
    sobre los coeficientes), ajustamos por IRLS estable:
      1) Poisson GLM  → mu, y dispersión de Pearson
      2) alpha (dispersión NB) por método de momentos
      3) NB GLM con alpha fijo (IRLS, converge de forma fiable) → mu final
    Devuelve (mu, alpha, engine).
    """
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    d = de.copy()
    d["sig_i"] = d["ontarget_significant"].astype(int)
    d["off_i"] = d["offtarget_flag"].astype(int)
    formula = "n_downstream ~ C(culture_condition) + sig_i + off_i"
    pois = smf.glm(formula, data=d, family=sm.families.Poisson()).fit()
    mu0 = np.asarray(pois.predict(d), float)
    y = d["n_downstream"].to_numpy(float)
    # método de momentos para la dispersión NB:  Var = mu + alpha*mu^2
    alpha = float(np.clip(np.mean(((y - mu0) ** 2 - mu0) / (mu0 ** 2)), 1e-3, 50))
    nb = smf.glm(formula, data=d, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
    mu = np.asarray(nb.predict(d), float)
    return mu, alpha, "NB GLM (alpha por método de momentos, IRLS)"


def empirical_bayes_gene(de, mu):
    """Shrinkage EB normal-normal del efecto por gen sobre la desviación log-rate.

    GATE DE CALIDAD: el efecto de un gen se estima SOLO con sus filas de knockdown
    on-target significativo. Un gen con mucho downstream pero sin KD validado es
    sospechoso (artefacto / off-target / diana no silenciada) y NO debe rankear alto.
    Genes sin ninguna condición significativa quedan fuera del ranking.

    work_i = log(y_i + .5) - log(mu_i + .5)   (desviación del baseline de efectos fijos)
    Por gen g (solo filas sig):  d_g = mean(work),  s2_g = sigma2_e / n_g
    Prior:      u_g ~ Normal(0, tau2)   con tau2 por método de momentos
    Posterior:  u_g | datos ~ Normal(shrink*d_g,  shrink*s2_g)
    """
    y = de["n_downstream"].to_numpy(float)
    work = np.log(y + 0.5) - np.log(mu + 0.5)
    df = pd.DataFrame({"gene": de["target_contrast_gene_name"].values, "work": work,
                       "logmu": np.log(mu + 0.5), "sig": de["ontarget_significant"].values})
    sigma2_e = float(np.var(work))               # ruido estimado sobre todas las filas
    df = df[df["sig"]].copy()                     # gate: solo KD significativo
    g = df.groupby("gene")
    d_g = g["work"].mean()
    n_g = g["work"].size()
    logmu_g = g["logmu"].mean()

    s2_g = sigma2_e / n_g
    tau2 = max(float(np.var(d_g) - s2_g.mean()), 1e-6)   # método de momentos
    shrink = tau2 / (tau2 + s2_g)
    u = shrink * d_g                             # media posterior (log-rate)
    post_var = shrink * s2_g
    sd = np.sqrt(post_var)

    thr = np.quantile(u, 0.99)                    # umbral top 1% sobre estimaciones puntuales
    p_top = pd.Series(norm.sf(thr, loc=u, scale=sd), index=u.index)
    expected = np.exp(logmu_g + u) - 0.5          # downstream esperado del gen

    out = pd.DataFrame({
        "regpower_eb_mean": u.round(4),
        "regpower_eb_sd": sd.round(4),
        "p_top_1pct": p_top.round(4),
        "n_conditions": n_g,
        "expected_downstream": expected.clip(lower=0).round(1),
        "shrinkage": shrink.round(3),
    })
    return out.sort_values("regpower_eb_mean", ascending=False), tau2, sigma2_e


def main():
    de = pd.read_csv(DATA / "DE_stats.suppl_table.csv")
    de["ontarget_significant"] = de["ontarget_significant"].astype(bool)
    de["offtarget_flag"] = de["offtarget_flag"].astype(bool)

    # features exploratorias por gen (reutiliza la lógica del EDA)
    gg = de.groupby("target_contrast_gene_name")
    feat = pd.DataFrame({
        "n_signif_conditions": gg["ontarget_significant"].sum().astype(int),
        "xcond_reproducibility": (gg["n_downstream"].min() /
                                  gg["n_downstream"].max().replace(0, np.nan)).fillna(0).round(3),
        "any_offtarget": gg["offtarget_flag"].any(),
    })

    print("== Modelo 2 · empirical-Bayes (pseudo-bayesiano) ==")
    mu, alpha, engine = fit_fixed_effects(de)
    print(f"  efectos fijos: {engine} · alpha_NB={alpha:.3f}")
    rank, tau2, sig2 = empirical_bayes_gene(de, mu)
    rank = rank.join(feat)
    rank.insert(0, "rank", np.arange(1, len(rank) + 1))
    print(f"  tau2 (prior)={tau2:.4f} · sigma2_e (ruido)={sig2:.4f} · genes={len(rank):,}")

    out_full = TAB / "hub_ranking_bayes.csv"
    rank.reset_index().rename(columns={"gene": "target_contrast_gene_name"}).to_csv(out_full, index=False)
    print("  tabla →", out_full.relative_to(ROOT))

    # ---- top_regulators_for_review.csv (judge-facing) ----
    # fila = gen en su condición de mayor efecto observado ENTRE las significativas
    peak = (de[de["ontarget_significant"]].sort_values("n_downstream", ascending=False)
              .drop_duplicates("target_contrast_gene_name")
              .set_index("target_contrast_gene_name"))
    top = rank.head(30).copy()
    rows = []
    for gene, r in top.iterrows():
        pk = peak.loc[gene]
        note = (f"KD sig {int(r.n_signif_conditions)}/3 cond; "
                f"repro cross-cond {r.xcond_reproducibility:.2f}; "
                f"{'posible off-target' if r.any_offtarget else 'sin off-target'}")
        rows.append({
            "rank": int(r["rank"]),
            "gene": gene,
            "condition": pk["culture_condition"],
            "regpower_eb_mean": r.regpower_eb_mean,
            "regpower_eb_sd": r.regpower_eb_sd,
            "p_top_1pct": r.p_top_1pct,
            "observed_n_downstream": int(pk["n_downstream"]),
            "expected_downstream": r.expected_downstream,
            "ontarget_significant": bool(pk["ontarget_significant"]),
            "single_guide_estimate": "NA (requiere DE_stats.h5ad)",
            "n_guides": "NA (requiere DE_stats.h5ad)",
            "offtarget_flag": bool(pk["offtarget_flag"]),
            "interpretation_note": note,
        })
    review = pd.DataFrame(rows)
    out_rev = TAB / "top_regulators_for_review.csv"
    review.to_csv(out_rev, index=False)
    print("  tabla →", out_rev.relative_to(ROOT))
    print("\n  TOP 12 reguladores (por regpower EB):")
    print(review.head(12)[["rank", "gene", "condition", "regpower_eb_mean",
                            "p_top_1pct", "observed_n_downstream"]].to_string(index=False))

    # ---- figura 07 ----
    plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white",
        "axes.edgecolor": "#cfd6e0", "font.size": 10, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": "#eceff4"})
    t = rank.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.errorbar(t["regpower_eb_mean"], range(len(t)), xerr=t["regpower_eb_sd"],
                fmt="o", color=ACCENT, ecolor=MUT, elinewidth=1.4, capsize=3, markersize=6)
    ax.set_yticks(range(len(t))); ax.set_yticklabels(t.index)
    ax.axvline(0, color="#cfd6e0", lw=1)
    ax.set_xlabel("poder regulatorio (log-rate, media posterior EB ± sd)")
    ax.set_title("Top 20 reguladores — ranking empirical-Bayes con incertidumbre",
                 fontweight="bold", fontsize=11)
    fig.tight_layout(); fig.savefig(FIG / "07_hub_posterior_ranking.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("  figura → docs/figures/07_hub_posterior_ranking.png")
    print("\n✓ Modelo 2 completo.")


if __name__ == "__main__":
    main()
