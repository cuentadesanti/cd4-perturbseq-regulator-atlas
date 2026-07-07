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

    reproducibility_sensitivity_audit(de, rank)   # auditoría opcional (no reemplaza el core)
    print("\n✓ Modelo 2 completo.")


def reproducibility_sensitivity_audit(de, core_rank):
    """Auditoría de SENSIBILIDAD guide/donor-aware.

    NO reestima el posterior EB ni es un modelo nuevo: repondera el score EB con
    reproducibilidad REAL del .obs de DE_stats.h5ad (cross-guide amplia, cross-donor
    parcial) y observa qué reguladores sobreviven. Solo corre si existe el metadata;
    el ranking core NO depende de esto."""
    meta_path = TAB / "de_obs_reproducibility_metadata.csv"
    if not meta_path.exists():
        print("  [sensitivity audit] omitido (no existe de_obs_reproducibility_metadata.csv). "
              "El core no depende de esto.")
        return
    AMBER = "#e0a441"
    meta = pd.read_csv(meta_path).drop(columns=["target_contrast_gene_name"], errors="ignore")
    sig = de[de["ontarget_significant"]][
        ["target_contrast_gene_name", "target_contrast", "culture_condition", "n_downstream"]].copy()
    m = sig.merge(meta, on=["target_contrast", "culture_condition"], how="left")

    # --- coverage del join (explícito: cross-donor es parcial) ---
    cov = pd.DataFrame([{
        "n_rows_de": len(de),
        "n_rows_sig": len(sig),
        "n_rows_joined": int(m["single_guide_estimate"].notna().sum()),
        "pct_single_guide_available": round(m["single_guide_estimate"].notna().mean() * 100, 1),
        "pct_guide_corr_available": round(m["guide_correlation_all"].notna().mean() * 100, 1),
        "pct_donor_corr_available": round(m["donor_correlation_hits_mean"].notna().mean() * 100, 1),
    }])
    cov.to_csv(TAB / "reproducibility_coverage.csv", index=False)
    print(f"  [sensitivity audit] coverage: guide_corr {cov.pct_guide_corr_available[0]:.0f}% · "
          f"donor_corr {cov.pct_donor_corr_available[0]:.0f}% (parcial) → peso neutral donde falta")

    def agg(g):
        return pd.Series({
            "guide_corr": g["guide_correlation_all"].mean(),
            "donor_corr": g["donor_correlation_hits_mean"].mean(),
            "single_guide_frac": g["single_guide_estimate"].astype(float).mean(),
            "n_guides_med": g["n_guides"].median(),
            "peak_condition": g.loc[g["n_downstream"].idxmax(), "culture_condition"],
        })
    pg = m.groupby("target_contrast_gene_name").apply(agg, include_groups=False)

    pg["guide_score"] = pg["guide_corr"].clip(0, 1)
    pg["donor_score"] = pg["donor_corr"].clip(0, 1)
    pg["single_guide_penalty"] = (1 - 0.4 * pg["single_guide_frac"].fillna(1)).round(3)
    # reweighting (neutral 0.75 donde falta la métrica → no penaliza por dato ausente)
    gt = (0.5 + 0.5 * pg["guide_score"]).fillna(0.75)
    dt = (0.5 + 0.5 * pg["donor_score"]).fillna(0.75)
    pg["repro_weight"] = (gt * dt * pg["single_guide_penalty"]).round(3)
    pg["donor_metadata"] = np.where(pg["donor_corr"].notna(), "presente", "ausente (peso neutral)")

    out = core_rank[["rank", "regpower_eb_mean"]].rename(columns={"rank": "eb_rank"}).join(pg, how="left")
    out["repro_weight"] = out["repro_weight"].fillna(0.75)          # genes sin metadata → neutral
    out["donor_metadata"] = out["donor_metadata"].fillna("ausente (peso neutral)")
    out["reweighted_score"] = (out["regpower_eb_mean"] * out["repro_weight"]).round(4)
    out = out.sort_values("reweighted_score", ascending=False)
    out["new_rank"] = np.arange(1, len(out) + 1)
    out = out.rename(columns={"eb_rank": "old_rank"})
    out["rank_shift"] = out["old_rank"] - out["new_rank"]           # + = subió con reproducibilidad
    out["guide_corr"] = out["guide_corr"].round(3)
    out["donor_corr"] = out["donor_corr"].round(3)
    out["single_guide_frac"] = out["single_guide_frac"].round(2)

    full = out.reset_index()
    full = full.rename(columns={full.columns[0]: "gene"})
    full[["gene", "old_rank", "new_rank", "rank_shift", "regpower_eb_mean", "guide_corr",
          "donor_corr", "single_guide_frac", "repro_weight", "reweighted_score",
          "donor_metadata", "n_guides_med", "peak_condition"]].to_csv(
        TAB / "hub_ranking_bayes_reproducibility_aware.csv", index=False)

    # --- audit interpretable: unión de top-30 de ambos rankings, con razón textual ---
    uni = full[(full["old_rank"] <= 30) | (full["new_rank"] <= 30)].copy()
    def classify(r):
        if r["old_rank"] <= 30 and r["new_rank"] <= 30:
            st = "sobrevive"
        elif r["old_rank"] <= 30 and r["new_rank"] > 30:
            st = "demotado"
        else:
            st = "promovido"
        if st == "promovido":
            rs = "promoted: reproducibilidad guide/donor alta (subvalorada por el core)"
        elif st == "demotado":
            if r["single_guide_frac"] >= 0.5:
                rs = "demoted: single-guide (sin chequeo cross-guide)"
            elif pd.notna(r["guide_corr"]) and r["guide_corr"] < 0.3:
                rs = "demoted: baja correlación cross-guide"
            else:
                rs = "demoted: reproducibilidad más baja que sus pares"
        else:  # sobrevive
            rs = ("survives: EB alto; donor metadata ausente (peso neutral)"
                  if r["donor_metadata"].startswith("ausente")
                  else "survives: EB alto y reproducible (guide + donor)")
        return pd.Series({"status": st, "reason": rs})
    uni[["status", "reason"]] = uni.apply(classify, axis=1)
    audit = uni[["gene", "old_rank", "new_rank", "status", "rank_shift", "guide_corr",
                 "donor_corr", "single_guide_frac", "donor_metadata", "reason"]]
    audit.sort_values("new_rank").to_csv(TAB / "reproducibility_audit.csv", index=False)

    # top-30 guide/donor-aware (judge-facing) con las columnas del join
    top = full.sort_values("new_rank").head(30)[
        ["gene", "peak_condition", "old_rank", "new_rank", "rank_shift", "guide_corr",
         "donor_corr", "single_guide_frac", "repro_weight", "reweighted_score", "donor_metadata"]]
    top.to_csv(TAB / "top_regulators_reproducibility_aware.csv", index=False)

    n_dem = int((uni["status"] == "demotado").sum())
    n_pro = int((uni["status"] == "promovido").sum())
    print(f"  [sensitivity audit] join OK · {n_dem} demotados / {n_pro} promovidos del top-30")
    print("      tablas → hub_ranking_bayes_reproducibility_aware.csv, "
          "top_regulators_reproducibility_aware.csv, reproducibility_audit.csv, reproducibility_coverage.csv")

    # figura 19 — slopegraph EB core -> reweighted (guide/donor-aware) para el top-20 EB
    top20 = full.sort_values("old_rank").head(20)
    CAP = 40
    fig, ax = plt.subplots(figsize=(7.4, 6.6))
    for _, r in top20.iterrows():
        y0, y1 = min(int(r["old_rank"]), CAP), min(int(r["new_rank"]), CAP)
        col = ACCENT if r["new_rank"] <= 30 else AMBER
        ax.plot([0, 1], [y0, y1], "-o", color=col, lw=1.6, ms=5, alpha=.85)
        ax.text(-0.04, y0, r["gene"], ha="right", va="center", fontsize=8, color="#333")
        if r["new_rank"] > 30:
            ax.text(1.04, y1, f"#{int(r['new_rank'])}", ha="left", va="center", fontsize=7.5, color=col)
    ax.set_xlim(-0.45, 1.5); ax.set_ylim(CAP + 3, 0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["ranking EB\n(core)", "guide/donor-aware\n(EB score reponderado)"])
    ax.set_yticks([1, 10, 20, 30, 40]); ax.set_ylabel("posición en el ranking (top-40)")
    ax.set_title("Auditoría de sensibilidad — ¿quién sobrevive a la reproducibilidad real?",
                 fontweight="bold", fontsize=11)
    ax.grid(False)
    ax.text(0.5, CAP + 2.2, "amarillo = demotado al reponderar por reproducibilidad real (no es un modelo nuevo)",
            ha="center", fontsize=8, color=AMBER, style="italic")
    fig.tight_layout(); fig.savefig(FIG / "19_reproducibility_aware_ranking_shift.png", dpi=135, bbox_inches="tight")
    plt.close(fig)
    print("      figura → docs/figures/19_reproducibility_aware_ranking_shift.png")


if __name__ == "__main__":
    main()
