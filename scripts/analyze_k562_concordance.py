#!/usr/bin/env python3
"""Step 2 — cross-cell-type concordance: K562 (Replogle 2022) vs CD4 primary T cells.

Splits regulators into UNIVERSAL (concordant across cell types) vs T-SPECIFIC
(divergent), and cross-validates the Step-3 CP program names against an independent
dataset.

DATA (not in repo; ~375 MB, fits on disk — the pre-aggregated bulk product, NOT the
8.8 GB cell-level file):
    Replogle et al. 2022, "Genome-wide Perturb-seq", processed datasets deposit
    (figshare plus 10.25452/figshare.plus.20029387). Perturbation x gene Z-normalized
    signatures — exactly the concordance unit, so no pseudobulking of raw cells.

    curl -L -o data/cache/K562_gwps_normalized_bulk.h5ad \
        https://ndownloader.figshare.com/files/35773217   # K562_gwps_normalized_bulk_01.h5ad

METHOD (harmonization is load-bearing — K562 CRISPRi cancer line vs CD4 CRISPRi primary):
  - JOIN genes by Ensembl (intersect, not force); JOIN regulators by Ensembl at the
    GENE level (collapse Replogle's promoter/guide rows by unweighted mean — X has no
    NaN, so no imputation) and DROP non-targeting / core_control rows.
  - CD4 side: Rest-preferred logFC (K562 is a resting proliferating line, so Rest is the
    apples-to-apples state; fall back to the peak-|effect| condition where no Rest row).
  - Concordance = per-regulator Pearson over CO-MEASURED genes on PER-REGULATOR
    z-scored effect vectors (no impute-to-0: 0=="no change" != NaN=="not measured").
  - NULL: permute the regulator<->regulator correspondence ACROSS datasets (genes held
    aligned) — asks "does this pair concord more than two random regulators?".
  - Controls reported, not assumed: (2) power = spearman(concordance, min-KD-rank) both
    sides + a well-powered stratified split; (3) coverage = spearman(concordance,
    n_co_measured). NB the power control is CONSERVATIVE: an immune gene lowly expressed
    in K562 has weak K562 KD *because* it is T-specific, so the well-powered split is a
    FLOOR on true T-specificity, not the total.

Outputs: docs/tables/k562_concordance.csv, docs/figures/36_k562_universal_vs_specific.png
"""
from pathlib import Path
import h5py, numpy as np, pandas as pd
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"; FIG = ROOT / "docs" / "figures"
BULK = CACHE / "K562_gwps_normalized_bulk.h5ad"


def build_k562():
    h = h5py.File(BULK, "r")
    X = np.asarray(h["X"], np.float64)                       # (11258, 8248), 0 NaN
    gt = h["obs/gene_transcript"][:].astype(str)             # "{n}_{SYM}_{promoter}_{ENSG}"
    vg = h["var/gene_id"][:].astype(str)                     # Ensembl gene axis
    core = np.asarray(h["obs/core_control"][:], bool) if "obs/core_control" in h else np.zeros(len(gt), bool)
    ncf = np.nan_to_num(np.asarray(h["obs/num_cells_filtered"][:], float), nan=0.0) if "obs/num_cells_filtered" in h else np.zeros(len(gt))

    def ens(s):
        for p in s.split("_"):
            if p.startswith("ENSG"): return p
        return None
    reg = np.array([ens(s) for s in gt], dtype=object)
    ctrl = np.array([r is None for r in reg]) | core          # drop non-targeting / core controls
    keep = ~ctrl
    Xr, er, nc = X[keep], reg[keep].astype(str), ncf[keep]
    uniq, inv = np.unique(er, return_inverse=True)            # collapse to gene level (unweighted mean)
    s = np.zeros((len(uniq), Xr.shape[1])); c = np.zeros(len(uniq))
    np.add.at(s, inv, Xr); np.add.at(c, inv, 1.0)
    Keff = s / c[:, None]
    kcells = {r: float(nc[er == r].sum()) for r in uniq}
    return Keff, uniq.astype(str), vg, kcells, int(ctrl.sum())


def main():
    Keff, Kreg, vg, kcells, n_ctrl = build_k562()
    gpos = {g: i for i, g in enumerate(vg)}; kri = {r: i for i, r in enumerate(Kreg)}

    cd4var = pd.read_csv(CACHE / "fingerprint_var.csv"); cd4obs = pd.read_csv(CACHE / "fingerprint_obs.csv")
    lfc = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
    cd4obs["sig"] = cd4obs["ontarget_significant"].astype(str).str.strip().eq("True")
    cd4obs["tc"] = cd4obs["target_contrast"].astype(str)
    gid2col = {g: i for i, g in enumerate(cd4var["gene_id"].astype(str))}

    shared_genes = [g for g in vg if g in gid2col]            # Ensembl intersect
    shared_reg = [r for r in Kreg if r in set(cd4obs.loc[cd4obs.sig, "tc"])]
    print(f"controls dropped={n_ctrl} | shared genes={len(shared_genes)} | shared regulators={len(shared_reg)}")
    kcol = np.array([gpos[g] for g in shared_genes]); ccol = np.array([gid2col[g] for g in shared_genes])

    sig = cd4obs[cd4obs.sig & cd4obs.tc.isin(set(shared_reg))].copy(); sig["absz"] = sig["ontarget_effect_size"].abs()
    cd4row, kdc = {}, {}
    for tc, g in sig.groupby("tc"):
        src = g[g.culture_condition == "Rest"] if (g.culture_condition == "Rest").any() else g
        cd4row[tc] = src["absz"].idxmax(); kdc[tc] = float(src["absz"].max())
    regs = [r for r in shared_reg if r in cd4row]
    Ck = np.array([np.asarray(lfc[cd4row[r]], np.float64)[ccol] for r in regs])   # NO imputation
    Kk = np.array([Keff[kri[r]][kcol] for r in regs])
    finC, finK = np.isfinite(Ck), np.isfinite(Kk); co = finC & finK
    usable = co.sum(1) >= 100
    regs = list(np.array(regs)[usable]); Ck, Kk, co, finC, finK = Ck[usable], Kk[usable], co[usable], finC[usable], finK[usable]
    n_co = co.sum(1); n = len(regs); G = len(shared_genes)

    def zc(v, m): w = v[m]; return (w - w.mean()) / (w.std() + 1e-9)
    obs_r = np.array([(zc(Ck[i], co[i]) * zc(Kk[i], co[i])).mean() for i in range(n)])
    rng = np.random.default_rng(0); B = 300; perm = np.empty((n, B))
    for b in range(B):
        j = rng.permutation(n); bad = j == np.arange(n); j[bad] = (j[bad] + 1) % n
        for i in range(n):
            m = finC[i] & finK[j[i]]
            perm[i, b] = (zc(Ck[i], m) * zc(Kk[j[i]], m)).mean() if m.sum() >= 100 else 0.0
    nf = perm.ravel(); q95 = np.quantile(nf, .95); med = np.quantile(nf, .5); emp_p = (perm >= obs_r[:, None]).mean(1)
    cls = np.where((obs_r > q95) & (emp_p < .05), "universal", np.where(obs_r <= med, "T-specific", "intermediate"))

    KDc = np.array([kdc[r] for r in regs]); KDk = np.array([abs(Keff[kri[r]][gpos[r]]) if r in gpos else np.nan for r in regs])
    rk = lambda x: pd.Series(x).rank(pct=True).to_numpy()
    mk = np.fmin(rk(KDc), rk(KDk)); okm = np.isfinite(mk)
    wp = (rk(KDc) >= .5) & (rk(KDk) >= .5)
    print(f"NULL median={med:.3f} q95={q95:.3f} | SPLIT universal={int((cls=='universal').sum())} "
          f"T-specific={int((cls=='T-specific').sum())} intermediate={int((cls=='intermediate').sum())}")
    print(f"POWER spearman(concordance,min-KD)={spearmanr(obs_r[okm],mk[okm]).statistic:.3f} "
          f"| COVERAGE spearman(concordance,n_co)={spearmanr(obs_r,n_co).statistic:.3f}")
    print(f"well-powered n={int(wp.sum())}: universal={int((wp&(cls=='universal')).sum())} T-specific={int((wp&(cls=='T-specific')).sum())}")

    e2s = dict(zip(cd4obs.tc, cd4obs.target_contrast_gene_name.astype(str)))
    df = pd.DataFrame({"ensembl": regs, "symbol": [e2s.get(r, r) for r in regs], "pearson_z": obs_r.round(4),
                       "empirical_p": emp_p.round(4), "n_co_measured": n_co, "class": cls,
                       "kd_cd4_abs": np.round(KDc, 3), "kd_k562_ontarget_abs": np.round(KDk, 3),
                       "well_powered": wp, "k562_cells": [int(kcells.get(r, 0)) for r in regs]})
    hb = pd.read_csv(TAB / "hub_ranking_bayes.csv")
    df["donor_robust"] = df.symbol.isin(set(hb.loc[hb.get("donor_robust", False) == True, "target_contrast_gene_name"].astype(str)))
    lab = pd.read_csv(TAB / "operator_cp_factors_3106_labeled.csv")
    f1 = set(lab[lab.factor == 1].top_regulators.iloc[0].split(";")); f6 = set(lab[lab.factor == 6].top_regulators.iloc[0].split(";"))
    df["program"] = np.where(df.symbol.isin(f6), "factor6(activation+SAGA)", np.where(df.symbol.isin(f1), "factor1(proteostasis/WRC)", ""))
    df = df.sort_values("pearson_z", ascending=False)
    TAB.mkdir(parents=True, exist_ok=True); df.to_csv(TAB / "k562_concordance.csv", index=False)
    make_figure(df, q95); print("wrote k562_concordance.csv + figure 36")


def make_figure(df, q95):
    ACCENT, VIOLET, MUT, GREY = "#0a8f9c", "#6b5fc0", "#8b95a8", "#c8cdd6"
    plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "axes.edgecolor": "#cfd6e0",
        "font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#eceff4"})
    np.random.seed(0)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 5.4), gridspec_kw={"width_ratios": [1, 1.15]})
    x = df.pearson_z.values
    a.hist(x, bins=90, color=GREY, edgecolor="white", linewidth=.3)
    a.axvline(0, color=MUT, ls="--", lw=1.2); a.axvline(q95, color=ACCENT, lw=1.4)
    a.axvspan(q95, x.max(), color=ACCENT, alpha=.07); a.axvspan(x.min(), 0, color=VIOLET, alpha=.07)
    ym = a.get_ylim()[1]
    a.text(q95 + .006, ym * .92, "universal\n(> null q95)", color=ACCENT, fontsize=8.5, va="top")
    a.text(-.006, ym * .92, "T-specific\n(≤ null median)", color=VIOLET, fontsize=8.5, va="top", ha="right")
    a.set_xlabel("per-regulator concordance  (z-scored Pearson, K562 vs CD4-Rest)"); a.set_ylabel("regulators")
    a.set_xlim(-.15, .25); a.set_title("Cross-cell-type concordance of 6,407 regulators", fontweight="bold", fontsize=11)
    wpr = df.well_powered & df.donor_robust
    a.text(.02, .98, f"well-powered ∩ donor-robust (conservative floor):\n"
           f"{int((wpr&(df['class']=='universal')).sum())} universal  ·  {int((wpr&(df['class']=='T-specific')).sum())} T-specific",
           transform=a.transAxes, va="top", fontsize=8, color="#444",
           bbox=dict(boxstyle="round,pad=0.4", fc="#f6f8fb", ec="#dfe5ec"))
    prog = {"factor6(activation+SAGA)": ("Factor 6\nSAGA + early activation", 1),
            "factor1(proteostasis/WRC)": ("Factor 1\nproteostasis + immune WRC", 0)}
    colc = {"universal": ACCENT, "T-specific": VIOLET, "intermediate": MUT}
    b.axvline(0, color=MUT, ls="--", lw=1); b.axvline(q95, color=ACCENT, lw=1.2)
    for key, (lab_, y) in prog.items():
        sub = df[df.program == key].sort_values("pearson_z").reset_index(drop=True)
        for i, r in sub.iterrows():
            yj = y + (0.20 if i % 2 else -0.20)
            b.scatter(r.pearson_z, yj, s=150 if r.donor_robust else 70, color=colc[r["class"]],
                      edgecolors="black" if r.well_powered else "none", linewidths=1.6, alpha=.9, zorder=3)
            b.annotate(r.symbol, (r.pearson_z, yj), fontsize=7.6, xytext=(0, 11 if i % 2 else -15),
                       textcoords="offset points", ha="center", color="#333", zorder=4)
    b.set_yticks([0, 1]); b.set_yticklabels([prog["factor1(proteostasis/WRC)"][0], prog["factor6(activation+SAGA)"][0]], fontsize=9)
    b.set_ylim(-.7, 1.7); b.set_xlim(-.09, .23); b.grid(axis="y")
    b.set_xlabel("concordance (K562 vs CD4)")
    b.set_title("Naming ↔ K562 cross-validation: SAGA universal, immune-WRC T-specific", fontweight="bold", fontsize=10.5)
    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], marker='o', color='w', markerfacecolor=ACCENT, markersize=9, label='universal'),
           Line2D([0], [0], marker='o', color='w', markerfacecolor=VIOLET, markersize=9, label='T-specific'),
           Line2D([0], [0], marker='o', color='w', markerfacecolor=MUT, markersize=9, label='intermediate'),
           Line2D([0], [0], marker='o', color='w', markerfacecolor='#ccc', markeredgecolor='black', markersize=9, label='well-powered (edge)'),
           Line2D([0], [0], marker='o', color='w', markerfacecolor='#ccc', markersize=12, label='donor-robust (large)')]
    b.legend(handles=leg, frameon=False, fontsize=7.5, loc="lower right")
    FIG.mkdir(parents=True, exist_ok=True); fig.tight_layout()
    fig.savefig(FIG / "36_k562_universal_vs_specific.png", dpi=140, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
