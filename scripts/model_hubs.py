#!/usr/bin/env python3
"""Model 2 — regulator ranking by empirical Bayes (pseudo-Bayesian).

IMPORTANT: this is NOT a full hierarchical NB. There is no PPL, no formal random
effects, no jointly sampled posterior. It is a fixed-effects NB (statsmodels)
+ a normal-normal empirical-Bayes shrinkage of the per-gene effect. Documented as such.

Runs from DE_stats.suppl_table.csv alone (local, ~15 MB). No downloads.

Outputs:
    docs/tables/hub_ranking_bayes.csv          (all genes, full ranking)
    docs/tables/top_regulators_for_review.csv  (top 30, judge-facing)
    docs/figures/07_hub_posterior_ranking.png

Usage:
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
    """Conditional mean (baseline) via a fixed-effects GLM: condition + KD quality.

    Poisson and NB share the SAME mean model; since we only use mu (no inference
    on the coefficients), we fit by stable IRLS:
      1) Poisson GLM  → mu, and Pearson dispersion
      2) alpha (NB dispersion) by method of moments
      3) NB GLM with fixed alpha (IRLS, converges reliably) → final mu
    Returns (mu, alpha, engine).
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
    # method of moments for the NB dispersion:  Var = mu + alpha*mu^2
    alpha = float(np.clip(np.mean(((y - mu0) ** 2 - mu0) / (mu0 ** 2)), 1e-3, 50))
    nb = smf.glm(formula, data=d, family=sm.families.NegativeBinomial(alpha=alpha)).fit()
    mu = np.asarray(nb.predict(d), float)
    return mu, alpha, "NB GLM (alpha by method of moments, IRLS)"


def empirical_bayes_gene(de, mu):
    """Normal-normal EB shrinkage of the per-gene effect on the log-rate deviation.

    QUALITY GATE: a gene's effect is estimated ONLY from its rows with a significant
    on-target knockdown. A gene with lots of downstream but no validated KD is
    suspicious (artifact / off-target / target not silenced) and should NOT rank high.
    Genes with no significant condition are excluded from the ranking.

    work_i = log(y_i + .5) - log(mu_i + .5)   (deviation from the fixed-effects baseline)
    Per gene g (sig rows only):  d_g = mean(work),  s2_g = sigma2_e / n_g
    Prior:      u_g ~ Normal(0, tau2)   with tau2 by method of moments
    Posterior:  u_g | data ~ Normal(shrink*d_g,  shrink*s2_g)
    """
    y = de["n_downstream"].to_numpy(float)
    work = np.log(y + 0.5) - np.log(mu + 0.5)
    df = pd.DataFrame({"gene": de["target_contrast_gene_name"].values, "work": work,
                       "logmu": np.log(mu + 0.5), "sig": de["ontarget_significant"].values})
    sigma2_e = float(np.var(work))               # noise estimated over all rows
    df = df[df["sig"]].copy()                     # gate: significant KD only
    g = df.groupby("gene")
    d_g = g["work"].mean()
    n_g = g["work"].size()
    logmu_g = g["logmu"].mean()

    s2_g = sigma2_e / n_g
    tau2 = max(float(np.var(d_g) - s2_g.mean()), 1e-6)   # method of moments
    shrink = tau2 / (tau2 + s2_g)
    u = shrink * d_g                             # posterior mean (log-rate)
    post_var = shrink * s2_g
    sd = np.sqrt(post_var)

    thr = np.quantile(u, 0.99)                    # top-1% threshold on point estimates
    p_top = pd.Series(norm.sf(thr, loc=u, scale=sd), index=u.index)
    expected = np.exp(logmu_g + u) - 0.5          # expected downstream for the gene

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

    # exploratory per-gene features (reuses the EDA logic)
    gg = de.groupby("target_contrast_gene_name")
    feat = pd.DataFrame({
        "n_signif_conditions": gg["ontarget_significant"].sum().astype(int),
        "xcond_reproducibility": (gg["n_downstream"].min() /
                                  gg["n_downstream"].max().replace(0, np.nan)).fillna(0).round(3),
        "any_offtarget": gg["offtarget_flag"].any(),
    })

    print("== Model 2 · empirical Bayes (pseudo-Bayesian) ==")
    mu, alpha, engine = fit_fixed_effects(de)
    print(f"  fixed effects: {engine} · alpha_NB={alpha:.3f}")
    rank, tau2, sig2 = empirical_bayes_gene(de, mu)
    rank = rank.join(feat)
    rank.insert(0, "rank", np.arange(1, len(rank) + 1))
    print(f"  tau2 (prior)={tau2:.4f} · sigma2_e (noise)={sig2:.4f} · genes={len(rank):,}")

    out_full = TAB / "hub_ranking_bayes.csv"
    rank.reset_index().rename(columns={"gene": "target_contrast_gene_name"}).to_csv(out_full, index=False)
    print("  table →", out_full.relative_to(ROOT))

    # ---- top_regulators_for_review.csv (judge-facing) ----
    # row = gene in its highest observed-effect condition AMONG the significant ones
    peak = (de[de["ontarget_significant"]].sort_values("n_downstream", ascending=False)
              .drop_duplicates("target_contrast_gene_name")
              .set_index("target_contrast_gene_name"))
    top = rank.head(30).copy()
    rows = []
    for gene, r in top.iterrows():
        pk = peak.loc[gene]
        note = (f"KD sig {int(r.n_signif_conditions)}/3 cond; "
                f"cross-cond repro {r.xcond_reproducibility:.2f}; "
                f"{'possible off-target' if r.any_offtarget else 'no off-target'}")
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
            "single_guide_estimate": "NA (requires DE_stats.h5ad)",
            "n_guides": "NA (requires DE_stats.h5ad)",
            "offtarget_flag": bool(pk["offtarget_flag"]),
            "interpretation_note": note,
        })
    review = pd.DataFrame(rows)
    out_rev = TAB / "top_regulators_for_review.csv"
    review.to_csv(out_rev, index=False)
    print("  table →", out_rev.relative_to(ROOT))
    print("\n  TOP 12 regulators (by EB regpower):")
    print(review.head(12)[["rank", "gene", "condition", "regpower_eb_mean",
                            "p_top_1pct", "observed_n_downstream"]].to_string(index=False))

    # ---- figure 07 ----
    plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white",
        "axes.edgecolor": "#cfd6e0", "font.size": 10, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": "#eceff4"})
    t = rank.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.errorbar(t["regpower_eb_mean"], range(len(t)), xerr=t["regpower_eb_sd"],
                fmt="o", color=ACCENT, ecolor=MUT, elinewidth=1.4, capsize=3, markersize=6)
    ax.set_yticks(range(len(t))); ax.set_yticklabels(t.index)
    ax.axvline(0, color="#cfd6e0", lw=1)
    ax.set_xlabel("regulatory power (log-rate, EB posterior mean ± sd)")
    ax.set_title("Top 20 regulators — empirical-Bayes ranking with uncertainty",
                 fontweight="bold", fontsize=11)
    fig.tight_layout(); fig.savefig(FIG / "07_hub_posterior_ranking.png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("  figure → docs/figures/07_hub_posterior_ranking.png")

    reproducibility_sensitivity_audit(de, rank)   # optional audit (does not replace the core)
    print("\n✓ Model 2 complete.")


def reproducibility_sensitivity_audit(de, core_rank):
    """Guide/donor-aware SENSITIVITY audit.

    Does NOT re-estimate the EB posterior and is not a new model: it reweights the EB
    score with REAL reproducibility from the .obs of DE_stats.h5ad (broad cross-guide,
    partial cross-donor) and observes which regulators survive. Runs only if the metadata
    exists; the core ranking does NOT depend on this."""
    meta_path = TAB / "de_obs_reproducibility_metadata.csv"
    if not meta_path.exists():
        print("  [sensitivity audit] skipped (de_obs_reproducibility_metadata.csv not found). "
              "The core does not depend on this.")
        return
    AMBER = "#e0a441"
    meta = pd.read_csv(meta_path).drop(columns=["target_contrast_gene_name"], errors="ignore")
    sig = de[de["ontarget_significant"]][
        ["target_contrast_gene_name", "target_contrast", "culture_condition", "n_downstream"]].copy()
    m = sig.merge(meta, on=["target_contrast", "culture_condition"], how="left")

    # --- join coverage (explicit: cross-donor is partial) ---
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
          f"donor_corr {cov.pct_donor_corr_available[0]:.0f}% (partial) → neutral weight where missing")

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
    # reweighting (neutral 0.75 where the metric is missing → no penalty for absent data)
    gt = (0.5 + 0.5 * pg["guide_score"]).fillna(0.75)
    dt = (0.5 + 0.5 * pg["donor_score"]).fillna(0.75)
    pg["repro_weight"] = (gt * dt * pg["single_guide_penalty"]).round(3)
    pg["donor_metadata"] = np.where(pg["donor_corr"].notna(), "present", "absent (neutral weight)")

    out = core_rank[["rank", "regpower_eb_mean"]].rename(columns={"rank": "eb_rank"}).join(pg, how="left")
    out["repro_weight"] = out["repro_weight"].fillna(0.75)          # genes without metadata → neutral
    out["donor_metadata"] = out["donor_metadata"].fillna("absent (neutral weight)")
    out["reweighted_score"] = (out["regpower_eb_mean"] * out["repro_weight"]).round(4)
    out = out.sort_values("reweighted_score", ascending=False)
    out["new_rank"] = np.arange(1, len(out) + 1)
    out = out.rename(columns={"eb_rank": "old_rank"})
    out["rank_shift"] = out["old_rank"] - out["new_rank"]           # + = rose with reproducibility
    out["guide_corr"] = out["guide_corr"].round(3)
    out["donor_corr"] = out["donor_corr"].round(3)
    out["single_guide_frac"] = out["single_guide_frac"].round(2)

    full = out.reset_index()
    full = full.rename(columns={full.columns[0]: "gene"})
    full[["gene", "old_rank", "new_rank", "rank_shift", "regpower_eb_mean", "guide_corr",
          "donor_corr", "single_guide_frac", "repro_weight", "reweighted_score",
          "donor_metadata", "n_guides_med", "peak_condition"]].to_csv(
        TAB / "hub_ranking_bayes_reproducibility_aware.csv", index=False)

    # --- interpretable audit: union of top-30 from both rankings, with a text reason ---
    uni = full[(full["old_rank"] <= 30) | (full["new_rank"] <= 30)].copy()
    def classify(r):
        if r["old_rank"] <= 30 and r["new_rank"] <= 30:
            st = "survives"
        elif r["old_rank"] <= 30 and r["new_rank"] > 30:
            st = "demoted"
        else:
            st = "promoted"
        if st == "promoted":
            rs = "promoted: high guide/donor reproducibility (undervalued by the core)"
        elif st == "demoted":
            if r["single_guide_frac"] >= 0.5:
                rs = "demoted: single-guide (no cross-guide check)"
            elif pd.notna(r["guide_corr"]) and r["guide_corr"] < 0.3:
                rs = "demoted: low cross-guide correlation"
            else:
                rs = "demoted: lower reproducibility than its peers"
        else:  # survives
            rs = ("survives: high EB; donor metadata absent (neutral weight)"
                  if r["donor_metadata"].startswith("absent")
                  else "survives: high EB and reproducible (guide + donor)")
        return pd.Series({"status": st, "reason": rs})
    uni[["status", "reason"]] = uni.apply(classify, axis=1)
    audit = uni[["gene", "old_rank", "new_rank", "status", "rank_shift", "guide_corr",
                 "donor_corr", "single_guide_frac", "donor_metadata", "reason"]]
    audit.sort_values("new_rank").to_csv(TAB / "reproducibility_audit.csv", index=False)

    # top-30 guide/donor-aware (judge-facing) with the join columns
    top = full.sort_values("new_rank").head(30)[
        ["gene", "peak_condition", "old_rank", "new_rank", "rank_shift", "guide_corr",
         "donor_corr", "single_guide_frac", "repro_weight", "reweighted_score", "donor_metadata"]]
    top.to_csv(TAB / "top_regulators_reproducibility_aware.csv", index=False)

    n_dem = int((uni["status"] == "demoted").sum())
    n_pro = int((uni["status"] == "promoted").sum())
    print(f"  [sensitivity audit] join OK · {n_dem} demoted / {n_pro} promoted from the top-30")
    print("      tables → hub_ranking_bayes_reproducibility_aware.csv, "
          "top_regulators_reproducibility_aware.csv, reproducibility_audit.csv, reproducibility_coverage.csv")

    # figure 19 — slopegraph EB core -> reweighted (guide/donor-aware) for the top-20 EB
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
    ax.set_xticklabels(["EB ranking\n(core)", "guide/donor-aware\n(reweighted EB score)"])
    ax.set_yticks([1, 10, 20, 30, 40]); ax.set_ylabel("rank position (top-40)")
    ax.set_title("Sensitivity audit — who survives real reproducibility?",
                 fontweight="bold", fontsize=11)
    ax.grid(False)
    ax.text(0.5, CAP + 2.2, "yellow = demoted when reweighting by real reproducibility (not a new model)",
            ha="center", fontsize=8, color=AMBER, style="italic")
    fig.tight_layout(); fig.savefig(FIG / "19_reproducibility_aware_ranking_shift.png", dpi=135, bbox_inches="tight")
    plt.close(fig)
    print("      figure → docs/figures/19_reproducibility_aware_ranking_shift.png")


if __name__ == "__main__":
    main()
