#!/usr/bin/env python3
"""Analysis 2 — the decisive specificity control (fully offline).

Question: is interferon (ISG) de-repression a SAGA specialty, or does *any*
transcription/chromatin perturbation de-repress interferon (generic stress)?

For three effect-size-comparable groups of regulators we build the convergent target set
(genes hit by >= half the group at |log2FC| > log2(1.5) on the cached log_fc, self-edges
removed — exactly as analyze_class_programs.py) and run the identical ISG hypergeometric test
(48-gene ISG core in the 10,282-gene measured universe):

  SAGA (positive class)          — SGF29, TADA1, TADA2B, SUPT20H, TAF6L, USP22, SUPT7L, ATXN7L3
  other chromatin/transcription  — non-SAGA remodelers / Polycomb / Pol II elongation / GTFs
  random (effect-size-matched)   — B random regulator groups matched to SAGA's n_downstream

Decisive comparison: SAGA fold vs other-chromatin fold vs the random null.

    python scripts/analyze_chromatin_stress_control.py

Outputs: docs/tables/chromatin_stress_control.csv · docs/figures/29_chromatin_stress_control.png
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
LOG2_1P5 = np.log2(1.5)
RNG = np.random.RandomState(0)
B_RANDOM = 500

SAGA = ["SGF29", "TADA1", "TADA2B", "SUPT20H", "TAF6L", "USP22", "SUPT7L", "ATXN7L3"]
# non-SAGA chromatin / general-transcription controls (MED1 excluded: Mediator is its own class)
OTHER_CHROM = ["SMARCA4", "SMARCB1", "ARID1A", "SMARCC1", "PBRM1",  # BAF / SWI-SNF
               "EZH2", "EED",                                       # Polycomb
               "CDK9", "CCNT1", "SUPT4H1",                          # Pol II elongation
               "TAF1", "TBP",                                       # general transcription
               "CHD4", "KMT2A", "BRD4"]                             # remodeler / MLL / bromodomain

# Curated type-I interferon core (same set as analyze_class_programs.py).
ISG_REF = {
    "IFI44", "IFI44L", "IFI6", "IFI35", "IFI27", "OAS1", "OAS2", "OAS3", "OASL", "MX1", "MX2",
    "BST2", "RIGI", "DDX58", "XAF1", "GBP1", "GBP2", "GBP4", "CMPK2", "HELZ2", "ISG15", "IRF7",
    "STAT1", "STAT2", "HERC5", "HERC6", "USP18", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3",
    "RSAD2", "SAMD9", "SAMD9L", "EPSTI1", "PARP9", "PARP12", "DTX3L", "LY6E", "PLSCR1", "APOL1",
    "APOL2", "APOL3", "APOL4", "APOL6", "TRIM22", "SP100", "IFI16",
}


def peak_rows(obs, genes, gcol="target_contrast_gene_name"):
    """Row index (into log_fc / obs) of each gene's peak-n_downstream *significant* KD."""
    sub = obs[obs[gcol].isin(genes) & obs["ontarget_significant"].astype(bool)]
    pk = sub.sort_values("n_downstream", ascending=False).drop_duplicates(gcol)
    return pk.set_index(gcol)


def convergent_targets(rows_df, lf, var_names, self_names):
    """Genes hit by >= half the group's members at |log2FC|>log2(1.5); self-edges removed.
    Returns (convergent_gene_set, frac_up_on_ISG_edges, n_members_used, median_n_downstream)."""
    members = list(rows_df.index)
    if len(members) < 2:
        return set(), None, len(members), None
    hit_counts = np.zeros(lf.shape[1], dtype=int)
    isg_idx = np.array([i for i, g in enumerate(var_names) if g in ISG_REF])
    isg_up, isg_tot = 0, 0
    for g in members:
        r = int(rows_df.loc[g, "row"])
        y = np.asarray(lf[r], dtype=float)
        hit = np.abs(y) > LOG2_1P5
        hit_counts += hit
        # direction on ISG edges this member hits
        m = hit[isg_idx]
        isg_up += int((y[isg_idx][m] > 0).sum())
        isg_tot += int(m.sum())
    thr = max(2, int(np.ceil(len(members) / 2)))
    conv_mask = hit_counts >= thr
    conv = {var_names[i] for i in np.where(conv_mask)[0]} - set(self_names)
    frac_up = round(isg_up / isg_tot, 3) if isg_tot else None
    med_nd = float(np.median(rows_df["n_downstream"]))
    return conv, frac_up, len(members), med_nd


def isg_test(conv, universe, isg_u):
    N, K = len(universe), len(isg_u)
    t = conv & universe
    n = len(t)
    k = len(t & isg_u)
    p = float(hypergeom.sf(k - 1, N, K, n)) if n else 1.0
    fold = (k / n) / (K / N) if n else 0.0
    return n, k, round(fold, 2), p


def main():
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv")
    var = pd.read_csv(CACHE / "fingerprint_var.csv")
    var_names = var["gene_name"].values
    lf = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
    obs = obs.reset_index(drop=True)
    obs["row"] = np.arange(len(obs))
    universe = set(var_names)
    isg_u = ISG_REF & universe
    gcol = "target_contrast_gene_name"

    saga_rows = peak_rows(obs, SAGA)
    saga_rows = saga_rows.assign(row=saga_rows["row"])[["row", "n_downstream"]]
    chrom_present = [g for g in OTHER_CHROM if g in set(obs[gcol])]
    chrom_rows = peak_rows(obs, chrom_present)[["row", "n_downstream"]]
    all_test = set(SAGA) | set(chrom_present)

    print(f"universe={len(universe)} · ISG in universe={len(isg_u)}")
    print(f"SAGA members used={len(saga_rows)} (median n_downstream={saga_rows.n_downstream.median():.0f})")
    print(f"other-chromatin used={len(chrom_rows)}: {chrom_present}")

    rows = []
    # SAGA
    conv, up, nm, nd = convergent_targets(saga_rows, lf, var_names, all_test)
    n, k, fold, p = isg_test(conv, universe, isg_u)
    saga_fold = fold
    rows.append({"group": "SAGA/chromatin", "n_members": nm, "median_n_downstream": nd,
                 "n_convergent_targets": n, "n_ISG": k, "isg_fold": fold, "isg_p": p,
                 "frac_isg_up_on_KD": up})
    print(f"  SAGA           : targets={n} ISG={k} fold={fold}x p={p:.2e} frac_up={up}")
    # other-chromatin
    conv, up, nm, nd = convergent_targets(chrom_rows, lf, var_names, all_test)
    n, k, fold, p = isg_test(conv, universe, isg_u)
    chrom_fold = fold
    rows.append({"group": "other chromatin/transcription", "n_members": nm, "median_n_downstream": nd,
                 "n_convergent_targets": n, "n_ISG": k, "isg_fold": fold, "isg_p": p,
                 "frac_isg_up_on_KD": up})
    print(f"  other-chromatin: targets={n} ISG={k} fold={fold}x p={p:.2e} frac_up={up}")

    # random null — effect-size matched to SAGA's n_downstream range, same group size
    lo, hi = saga_rows.n_downstream.min(), saga_rows.n_downstream.max()
    sig = obs[obs["ontarget_significant"].astype(bool)]
    pool = (sig.sort_values("n_downstream", ascending=False).drop_duplicates(gcol))
    pool = pool[(pool["n_downstream"] >= lo * 0.8) & (pool["n_downstream"] <= hi * 1.2)]
    pool = pool[~pool[gcol].isin(all_test)]
    pool_idx = pool[[gcol, "row", "n_downstream"]].reset_index(drop=True)
    print(f"  random pool (matched n_downstream {lo}-{hi}): {len(pool_idx)} genes · {B_RANDOM} draws")
    null_folds, null_isgN, null_up = [], [], []
    gsz = len(saga_rows)
    for b in range(B_RANDOM):
        pick = pool_idx.iloc[RNG.choice(len(pool_idx), gsz, replace=False)]
        rdf = pick.set_index(gcol)[["row", "n_downstream"]]
        conv, up, _, _ = convergent_targets(rdf, lf, var_names, all_test | set(pick[gcol]))
        n, k, fold, _ = isg_test(conv, universe, isg_u)
        null_folds.append(fold); null_isgN.append(k)
        if up is not None:
            null_up.append(up)
    null_folds = np.array(null_folds)
    null_up_mean = round(float(np.mean(null_up)), 3) if null_up else None
    nm_mean, nm_sd = float(null_folds.mean()), float(null_folds.std())
    saga_z = (saga_fold - nm_mean) / (nm_sd + 1e-9)
    saga_p = float((np.sum(null_folds >= saga_fold) + 1) / (B_RANDOM + 1))
    chrom_p = float((np.sum(null_folds >= chrom_fold) + 1) / (B_RANDOM + 1))
    rows.append({"group": "random (effect-size matched)", "n_members": gsz,
                 "median_n_downstream": float(pool_idx.n_downstream.median()),
                 "n_convergent_targets": None, "n_ISG": round(float(np.mean(null_isgN)), 1),
                 "isg_fold": round(nm_mean, 2), "isg_p": None, "frac_isg_up_on_KD": null_up_mean,
                 "null_sd": round(nm_sd, 2), "null_p95": round(float(np.percentile(null_folds, 95)), 2)})
    print(f"  random null    : fold mean={nm_mean:.2f}±{nm_sd:.2f} (95th pct {np.percentile(null_folds,95):.2f})")
    print(f"  SAGA vs null   : z={saga_z:+.1f} · p_perm={saga_p:.4f}")
    print(f"  chromatin vs null: p_perm={chrom_p:.4f}")

    df = pd.DataFrame(rows)
    df.to_csv(TAB / "chromatin_stress_control.csv", index=False)
    saga_up, chrom_up = rows[0]["frac_isg_up_on_KD"], rows[1]["frac_isg_up_on_KD"]
    # verdict: separate MAGNITUDE (fold) from DIRECTION (de-repression) specificity
    mag_specific = saga_p < 0.05 and saga_fold >= 2 * nm_mean
    dir_gap = (saga_up or 0) - (null_up_mean or 0.5)
    dir_specific = (saga_up or 0) >= 0.85 and dir_gap >= 0.25
    if mag_specific and dir_specific:
        verdict = "SAGA-specific in both magnitude and direction"
    elif saga_p < 0.05 or dir_gap >= 0.15:
        verdict = ("Largely general — the ISG-enrichment MAGNITUDE is mostly a strong-perturbation-"
                   "under-stimulation effect (matched random ~%.1fx, other-chromatin %.1fx). SAGA is only "
                   "marginally above random in magnitude (%.1fx, p=%.3f) but the MOST consistently "
                   "de-repressive (%.0f%% ISG-up on KD vs %.0f%% random, %.0f%% other-chromatin). "
                   "Requalify: nominate SAGA as a de-repressive control point, do not claim strong specificity."
                   % (nm_mean, chrom_fold, saga_fold, saga_p, saga_up * 100,
                      (null_up_mean or 0.5) * 100, (chrom_up or 0.5) * 100))
    elif chrom_p < 0.05 and chrom_fold >= 2 * nm_mean:
        verdict = "general chromatin/transcription stress (not SAGA-specific)"
    else:
        verdict = "no enrichment beyond matched random"
    summary = {"saga_fold": saga_fold, "other_chromatin_fold": chrom_fold,
               "random_null_mean_fold": round(nm_mean, 2), "random_null_sd": round(nm_sd, 2),
               "saga_z_vs_null": round(float(saga_z), 2), "saga_p_vs_null": saga_p,
               "chromatin_p_vs_null": chrom_p, "n_random_draws": B_RANDOM,
               "frac_isg_up_saga": saga_up, "frac_isg_up_chromatin": chrom_up, "frac_isg_up_random": null_up_mean,
               "magnitude_specific": bool(mag_specific), "direction_specific": bool(dir_specific),
               "verdict": verdict}
    json.dump(summary, open(TAB / "chromatin_stress_control_summary.json", "w"), indent=2)
    print(f"\n  VERDICT: {verdict}")

    _figure(null_folds, saga_fold, chrom_fold, saga_up, chrom_up, null_up_mean)
    print("  figure → docs/figures/29_chromatin_stress_control.png")


def _figure(null_folds, saga_fold, chrom_fold, saga_up, chrom_up, null_up_mean):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    labels = ["SAGA/\nchromatin", "other chromatin/\ntranscription", "random\n(matched)"]
    colors = ["#8e44ad", "#e0a441", "#8b95a8"]
    x = np.arange(3)
    nm = float(null_folds.mean())
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

    # A — magnitude (ISG fold): SAGA barely above the matched-random null
    folds = [saga_fold, chrom_fold, nm]
    ax1.bar(x, folds, color=colors, width=0.62, zorder=2)
    p95 = float(np.percentile(null_folds, 95)); p5 = float(np.percentile(null_folds, 5))
    ax1.errorbar(2, nm, yerr=[[nm - p5], [p95 - nm]], fmt="none", ecolor="#5d6779", capsize=5, lw=1.4, zorder=3)
    ax1.axhline(p95, color="#5d6779", ls=":", lw=1)
    ax1.text(2.45, p95, "random 95th", fontsize=7.5, color="#5d6779", va="center")
    for xi, f in zip(x, folds):
        ax1.text(xi, f + max(folds) * 0.02, f"{f:.1f}×", ha="center", fontsize=10, fontweight="bold")
    ax1.set_xticks(x); ax1.set_xticklabels(labels, fontsize=9)
    ax1.set_ylabel("ISG fold-enrichment of convergent targets")
    ax1.set_title("Magnitude — mostly a general strong-perturbation effect", fontsize=9.5, fontweight="bold", loc="left")
    ax1.set_ylim(0, max(folds) * 1.2)

    # B — direction (fraction of ISG edges that are de-repressive on KD): the SAGA-specific signal
    ups = [saga_up or 0, chrom_up or 0, null_up_mean or 0.5]
    ax2.bar(x, ups, color=colors, width=0.62, zorder=2)
    ax2.axhline(0.5, color="#999", ls="--", lw=1); ax2.text(2.45, 0.5, "no bias", fontsize=7.5, color="#999", va="center")
    for xi, u in zip(x, ups):
        ax2.text(xi, u + 0.02, f"{u*100:.0f}%", ha="center", fontsize=10, fontweight="bold")
    ax2.set_xticks(x); ax2.set_xticklabels(labels, fontsize=9)
    ax2.set_ylabel("Fraction of ISG hits UP on knockdown (de-repressive)")
    ax2.set_title("Direction — SAGA is the MOST consistent de-repressor", fontsize=9.5, fontweight="bold", loc="left")
    ax2.set_ylim(0, 1.05)

    fig.suptitle("Specificity control — the interferon enrichment is largely generic (magnitude);\n"
                 "SAGA's distinction is the consistency of de-repression, not a unique magnitude",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "29_chromatin_stress_control.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
