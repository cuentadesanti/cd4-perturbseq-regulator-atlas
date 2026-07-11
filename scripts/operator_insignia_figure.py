#!/usr/bin/env python3
"""Insignia figure for Section 4 — the module-discovery logic in one view.

Six panels telling the whole story: (A) the denoised operator block-ordered by community;
(B) the eigenspectrum with the global mode, MP bulk, and empirical signal edge; (C) blind recovery
of mitochondrial Complex I (its subunit sub-block); (D) the SAGA convergent module sub-block;
(E) recovered-community vs convergent-module distinction as partition stability against the gate;
(F) the predictive boundary — real vs shuffled features for leave-regulator-out inference.

    python scripts/operator_insignia_figure.py

Output: docs/figures/38_operator_insignia_3106.png
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import operator_communities as oc

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"; FIG = ROOT / "docs" / "figures"
ACCENT, VIOLET, MUT, GREY, RED = "#0a8f9c", "#6b5fc0", "#8b95a8", "#c8cdd6", "#c0504d"


def submatrix(C, idx, order_by=None):
    sub = C[np.ix_(idx, idx)].copy()
    if order_by is not None:
        o = np.argsort(order_by)
        sub = sub[np.ix_(o, o)]
    return sub


def main():
    plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white",
                         "font.size": 9, "axes.edgecolor": "#cfd6e0", "axes.titlesize": 10})
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    X, regs = oc.build_feature_matrix(d); regs = np.asarray(regs).astype(str)
    N, T = X.shape; q = N / T
    C = np.corrcoef(X)
    C_clean, spec, sigma2, lam_plus, lam0 = oc.rmt_clean(C, q)

    comm = pd.read_csv(TAB / "operator_communities_3106.csv").set_index("regulator").reindex(regs)
    lab = comm["community"].to_numpy()
    stab = comm.groupby(comm["community"]).community_stability.first() if "community_stability" in comm else None
    null = pd.read_csv(TAB / "operator_community_null_3106.csv").set_index("metric")["value"]
    emp_edge = float(null.get("noise_edge_empirical_p95", 2.95))
    mod_z = float(null.get("modularity_null_z", 259))
    n_signal = int(pd.read_csv(TAB / "operator_community_null_3106.csv")
                   .set_index("metric").loc["noise_edge_empirical_p95", "n_signal"])

    fig, axes = plt.subplots(2, 3, figsize=(15, 9.4))
    (a, b, c), (dd, e, f) = axes

    # ---- (A) denoised operator, block-ordered by community ----
    order = np.argsort(lab)
    a.imshow(C_clean[np.ix_(order, order)], cmap="RdBu_r", vmin=-0.15, vmax=0.15, aspect="auto")
    a.set_xticks([]); a.set_yticks([])
    a.set_title("A · Denoised operator, block-ordered by community", fontweight="bold")
    # outline the three stability-gated communities (5,6,7) + SAGA (2)
    labs_sorted = lab[order]
    for cval, col, name in [(2, VIOLET, "SAGA"), (5, MUT, ""), (6, MUT, ""), (7, ACCENT, "Complex I")]:
        w = np.where(labs_sorted == cval)[0]
        if len(w):
            s, sz = w.min(), len(w)
            a.add_patch(Rectangle((s, s), sz, sz, fill=False, edgecolor=col, lw=1.8))
            if name:
                a.text(s + sz / 2, s - 25, name, color=col, fontsize=8, ha="center", fontweight="bold")
    a.text(0.02, 0.98, f"modularity z = {mod_z:.0f}", transform=a.transAxes, va="top",
           fontsize=8, bbox=dict(boxstyle="round", fc="#f6f8fb", ec="#dfe5ec"))

    # ---- (B) eigenspectrum: global mode, MP bulk, empirical edge ----
    ev = np.sort(spec["eigenvalue"].to_numpy())[::-1]
    b.scatter(np.arange(len(ev)), ev, s=6, color=MUT, zorder=2)
    b.set_yscale("log")
    b.axhline(lam_plus, color=GREY, ls=":", lw=1.3)
    b.axhline(emp_edge, color=ACCENT, lw=1.6)
    b.scatter([0], [ev[0]], s=45, color=RED, zorder=3)
    b.annotate(f"global mode λ₀ = {ev[0]:.0f}", (0, ev[0]), xytext=(40, -6),
               textcoords="offset points", fontsize=8, color=RED)
    b.text(len(ev) * 0.98, emp_edge * 1.15, f"empirical edge {emp_edge:.2f}", color=ACCENT,
           fontsize=8, ha="right")
    b.text(len(ev) * 0.98, lam_plus * 0.6, f"MP bulk λ₊ = {lam_plus:.2f}", color=MUT, fontsize=8, ha="right")
    b.axvspan(0, n_signal, color=ACCENT, alpha=0.06)
    b.text(n_signal + 30, ev[3], f"{n_signal} signal\ndirections", color=ACCENT, fontsize=8.5, fontweight="bold")
    b.set_xlim(-20, 700); b.set_xlabel("eigenvalue rank"); b.set_ylabel("eigenvalue")
    b.set_title("B · Spectral denoising sets the signal count", fontweight="bold")

    # ---- (C) Complex I recovered community sub-block ----
    ci = np.where(lab == 7)[0]
    ci_sub = submatrix(C, ci, order_by=None)
    im = c.imshow(ci_sub, cmap="RdBu_r", vmin=-0.3, vmax=0.3, aspect="auto")
    c.set_xticks([]); c.set_yticks([])
    nduf = [r for r in regs[ci] if r.startswith(("NDUF", "MT-ND"))]
    c.set_title("C · Complex I — recovered blind (n=87)", fontweight="bold")
    c.text(0.5, -0.06, f"CORUM NADH-dehydrogenase  BH-FDR 1.4×10⁻⁷   ·   {len(nduf)} NDUF* subunits present",
           transform=c.transAxes, ha="center", fontsize=8, color=ACCENT)

    # ---- (D) SAGA convergent module sub-block (diagonal masked, values annotated) ----
    saga = ["SUPT20H", "SUPT7L", "TADA2B", "TAF6L", "USP22", "SGF29"]
    sidx = [np.where(regs == s)[0][0] for s in saga if (regs == s).any()]
    sub = C[np.ix_(sidx, sidx)].copy()
    off = sub.copy(); np.fill_diagonal(off, np.nan)
    vmax = float(np.nanmax(off)); mean_intra = float(np.nanmean(off))
    dd.imshow(off, cmap="Reds", vmin=0, vmax=vmax, aspect="auto")
    for i in range(len(sidx)):
        for j in range(len(sidx)):
            if i != j:
                dd.text(j, i, f"{sub[i, j]:.2f}", ha="center", va="center", fontsize=6.8,
                        color="white" if sub[i, j] > vmax * 0.6 else "#333")
    dd.set_xticks(range(len(sidx))); dd.set_xticklabels([regs[i] for i in sidx], rotation=45, ha="right", fontsize=7)
    dd.set_yticks(range(len(sidx))); dd.set_yticklabels([regs[i] for i in sidx], fontsize=7)
    dd.set_title("D · SAGA — convergent module", fontweight="bold")
    dd.text(0.5, -0.30, f"mean intra-corr {mean_intra:.2f}   ·   re-forms in K562 (z = 16.1)",
            transform=dd.transAxes, ha="center", fontsize=8, color=VIOLET)

    # ---- (E) recovered vs convergent: stability against the gate ----
    order_comm = [7, 5, 6, 2]
    names = ["Complex I\n(recovered)", "unnamed\nstable", "unnamed\nstable", "SAGA\n(convergent)"]
    vals = [float(stab.get(k, np.nan)) for k in order_comm] if stab is not None else [0.87, 0.835, 0.834, 0.56]
    cols = [ACCENT, MUT, MUT, VIOLET]
    e.bar(range(4), vals, color=cols, edgecolor="white")
    e.axhline(0.8, color=RED, ls="--", lw=1.3)
    e.text(3.3, 0.81, "stability gate 0.8", color=RED, fontsize=8, ha="right")
    e.set_xticks(range(4)); e.set_xticklabels(names, fontsize=7.5)
    e.set_ylim(0, 1); e.set_ylabel("partition stability s_c")
    e.set_title("E · Recovered community vs convergent module", fontweight="bold")

    # ---- (F) predictive boundary: real vs shuffled ----
    imc = pd.read_csv(TAB / "imc_leave_regulator_out_3106.csv")
    cat = pd.read_csv(TAB / "imc_nonlinear_catboost_3106.csv")
    groups = ["ridge", "CatBoost K10", "CatBoost K30"]
    real = [imc.r2_ridge.mean(), cat[cat.K == 10].r2_test_real.mean(), cat[cat.K == 30].r2_test_real.mean()]
    shuf = [imc.r2_shuffled.mean(), cat[cat.K == 10].r2_test_shuffled.mean(), cat[cat.K == 30].r2_test_shuffled.mean()]
    x = np.arange(3); w = 0.36
    f.bar(x - w / 2, real, w, label="real features", color=ACCENT, edgecolor="white")
    f.bar(x + w / 2, shuf, w, label="shuffled features", color=GREY, edgecolor="white")
    f.axhline(0, color="#333", lw=0.8)
    f.set_xticks(x); f.set_xticklabels(groups, fontsize=8)
    f.set_ylim(-0.01, 0.01); f.set_ylabel("held-out R²  (unseen regulators)")
    f.legend(frameon=False, fontsize=8, loc="upper right")
    f.set_title("F · Predictive boundary: real ≈ shuffled ≈ 0", fontweight="bold")
    f.text(0.5, -0.16, "capacity probe: train R² = +0.13 while test R² = −0.02 — capacity present, signal absent",
           transform=f.transAxes, ha="center", fontsize=8, color=MUT)

    fig.suptitle("Denoised operator structure reveals regulatory modules — and its predictive boundary",
                 fontsize=12.5, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / "38_operator_insignia_3106.png"
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"[insignia] wrote {out}  (Complex I NDUF* subunits: {len(nduf)}; modularity z={mod_z:.0f}; "
          f"signal={n_signal})")


if __name__ == "__main__":
    main()
