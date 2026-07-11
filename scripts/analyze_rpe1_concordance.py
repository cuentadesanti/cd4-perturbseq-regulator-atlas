#!/usr/bin/env python3
"""Third cell type — RPE1 (Replogle 2022) concordance vs CD4 primary T cells.

Robustness of the "universal (SAGA) vs T-specific (DOCK2)" split: today it rests on n=2 cell
types (CD4 vs K562), both extremes. RPE1 (epithelial, non-cancerous, hTERT, near-euploid, p53+)
is a genuinely third context. If SAGA reappears universal in RPE1 and T-specific regulators stay
flat, specificity goes from a pairwise comparison to a property sustained across 3 diverse contexts.

METHOD is IDENTICAL to scripts/analyze_k562_concordance.py (this is its RPE1 twin — same build,
same per-regulator z-scored Pearson over co-measured genes, same cross-dataset permutation null q95,
same universal/T-specific classification, same power + coverage controls, same well-powered floor).
Only the input file changes: rpe1_normalized_bulk_01.h5ad (figshare plus 20029387, file 35775512),
the essential-scale sibling of the K562 bulk product.

CHECKPOINT A gate (pre-registered, same rule as K562): <50 shared regulators -> STOP & document;
>=100 -> proceed; 50-100 -> proceed with a reduced-power caveat. RPE1 is ESSENTIAL-scale (~2.4k
targets), so also verify the anchors (6 SAGA subunits + DOCK2) are actually perturbed in RPE1 —
if not, the specificity test loses part of its point even when the global overlap passes.

    curl -L -o data/cache/rpe1_normalized_bulk.h5ad https://ndownloader.figshare.com/files/35775512
    python scripts/analyze_rpe1_concordance.py

Outputs: docs/tables/rpe1_concordance.csv, docs/tables/rpe1_checkpointA.txt
"""
from pathlib import Path
import h5py, numpy as np, pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"
BULK = CACHE / "rpe1_normalized_bulk.h5ad"
SAGA = ["SUPT20H", "TAF6L", "TADA2B", "SUPT7L", "USP22", "SGF29"]
TSPEC_ANCHOR = ["DOCK2"]


def _ens(s):
    for p in s.split("_"):
        if p.startswith("ENSG"):
            return p
    return None


def build_rpe1():
    """Collapse Replogle promoter/guide rows to GENE level (unweighted mean), drop controls.
    Schema-identical to build_k562()."""
    h = h5py.File(BULK, "r")
    X = np.asarray(h["X"], np.float64)                       # (2679, 8749), 0 NaN
    gt = h["obs/gene_transcript"][:].astype(str)             # "{n}_{SYM}_{promoter}_{ENSG}"
    vg = h["var/gene_id"][:].astype(str)                     # Ensembl gene axis
    core = np.asarray(h["obs/core_control"][:], bool) if "obs/core_control" in h else np.zeros(len(gt), bool)
    ncf = np.nan_to_num(np.asarray(h["obs/num_cells_filtered"][:], float), nan=0.0) if "obs/num_cells_filtered" in h else np.zeros(len(gt))
    reg = np.array([_ens(s) for s in gt], dtype=object)
    sym = np.array([s.split("_")[1] if len(s.split("_")) > 1 else s for s in gt], dtype=object)
    ctrl = np.array([r is None for r in reg]) | core          # drop non-targeting / core controls
    keep = ~ctrl
    Xr, er, nc = X[keep], reg[keep].astype(str), ncf[keep]
    uniq, inv = np.unique(er, return_inverse=True)            # collapse to gene level (unweighted mean)
    s = np.zeros((len(uniq), Xr.shape[1])); c = np.zeros(len(uniq))
    np.add.at(s, inv, Xr); np.add.at(c, inv, 1.0)
    Keff = s / c[:, None]
    kcells = {r: float(nc[er == r].sum()) for r in uniq}
    rpe1_syms = set(sym[keep].astype(str))
    return Keff, uniq.astype(str), vg, kcells, int(ctrl.sum()), rpe1_syms


def main():
    Keff, Kreg, vg, kcells, n_ctrl, rpe1_syms = build_rpe1()
    gpos = {g: i for i, g in enumerate(vg)}; kri = {r: i for i, r in enumerate(Kreg)}

    cd4var = pd.read_csv(CACHE / "fingerprint_var.csv"); cd4obs = pd.read_csv(CACHE / "fingerprint_obs.csv")
    lfc = np.load(CACHE / "log_fc.f32.npy", mmap_mode="r")
    cd4obs["sig"] = cd4obs["ontarget_significant"].astype(str).str.strip().eq("True")
    cd4obs["tc"] = cd4obs["target_contrast"].astype(str)
    gid2col = {g: i for i, g in enumerate(cd4var["gene_id"].astype(str))}
    e2s = dict(zip(cd4obs.tc, cd4obs.target_contrast_gene_name.astype(str)))

    shared_genes = [g for g in vg if g in gid2col]            # Ensembl intersect
    shared_reg = [r for r in Kreg if r in set(cd4obs.loc[cd4obs.sig, "tc"])]
    shared_syms = {e2s.get(r, r) for r in shared_reg}

    # ---------- CHECKPOINT A gate (write before building anything heavier) ----------
    n_shared = len(shared_reg)
    saga_rpe1 = [a for a in SAGA if a in rpe1_syms]; saga_shared = [a for a in SAGA if a in shared_syms]
    dock2_rpe1 = [a for a in TSPEC_ANCHOR if a in rpe1_syms]
    decision = ("STOP — overlap insuficiente (<50)" if n_shared < 50 else
                "PROCEED WITH CAVEAT (50-100, reduced power)" if n_shared < 100 else "PROCEED (>=100)")
    TAB.mkdir(parents=True, exist_ok=True)
    with open(TAB / "rpe1_checkpointA.txt", "w") as f:
        f.write("RPE1 Checkpoint A — regulator overlap gate (pre-registered)\n")
        f.write(f"RPE1 targeting genes (controls dropped {n_ctrl}): {len(Kreg)}\n")
        f.write(f"CD4 significant regulators (unique): {cd4obs.loc[cd4obs.sig, 'tc'].nunique()}\n")
        f.write(f"SHARED (RPE1 perturbed ∩ CD4 significant, by Ensembl): {n_shared}\n")
        f.write(f"shared genes (Ensembl intersect): {len(shared_genes)}\n\n")
        f.write("Anchor presence (the specificity-test anchors):\n")
        f.write(f"  SAGA subunits perturbed in RPE1: {saga_rpe1} ({len(saga_rpe1)}/6)\n")
        f.write(f"  SAGA subunits in shared set:     {saga_shared} ({len(saga_shared)}/6)\n")
        f.write(f"  DOCK2 (T-specific anchor) in RPE1: {'YES' if dock2_rpe1 else 'NO — not in essential-scale set'}\n\n")
        f.write(f"DECISION: {decision}\n")
        if not dock2_rpe1:
            f.write("CAVEAT: DOCK2 is not targeted in RPE1's essential-scale screen (it is a hematopoietic-\n"
                    "specific GEF, not core-essential), so its T-specific arm cannot be measured directly here.\n"
                    "Compensated below by a CLASS-level test: are K562-classified T-specific regulators that ARE\n"
                    "shared with RPE1 also non-universal in RPE1?\n")
    print(f"[Checkpoint A] shared={n_shared} | SAGA in RPE1 {len(saga_rpe1)}/6 {saga_rpe1} | "
          f"DOCK2 in RPE1: {'yes' if dock2_rpe1 else 'NO'} | -> {decision}")
    if n_shared < 50:
        print("[GATE FAILED] documented in rpe1_checkpointA.txt; STOP (not forced).")
        return

    # ---------- concordance pipeline (identical to K562) ----------
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
    n_co = co.sum(1); n = len(regs)

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
          f"T-specific={int((cls=='T-specific').sum())} intermediate={int((cls=='intermediate').sum())} (n={n})")
    print(f"POWER spearman(concordance,min-KD)={spearmanr(obs_r[okm],mk[okm]).statistic:.3f} "
          f"| COVERAGE spearman(concordance,n_co)={spearmanr(obs_r,n_co).statistic:.3f}")
    print(f"well-powered n={int(wp.sum())}: universal={int((wp&(cls=='universal')).sum())} "
          f"T-specific={int((wp&(cls=='T-specific')).sum())}")

    df = pd.DataFrame({"ensembl": regs, "symbol": [e2s.get(r, r) for r in regs], "pearson_z": obs_r.round(4),
                       "empirical_p": emp_p.round(4), "n_co_measured": n_co, "class": cls,
                       "kd_cd4_abs": np.round(KDc, 3), "kd_rpe1_ontarget_abs": np.round(KDk, 3),
                       "min_kd_rank": np.round(mk, 4), "well_powered": wp,
                       "rpe1_cells": [int(kcells.get(r, 0)) for r in regs]})
    hb = pd.read_csv(TAB / "hub_ranking_bayes.csv")
    df["donor_robust"] = df.symbol.isin(set(hb.loc[hb.get("donor_robust", False) == True, "target_contrast_gene_name"].astype(str)))
    lab = pd.read_csv(TAB / "operator_cp_factors_3106_labeled.csv")
    f1 = set(lab[lab.factor == 1].top_regulators.iloc[0].split(";")); f6 = set(lab[lab.factor == 6].top_regulators.iloc[0].split(";"))
    df["program"] = np.where(df.symbol.isin(f6), "factor6(activation+SAGA)", np.where(df.symbol.isin(f1), "factor1(proteostasis/WRC)", ""))
    df = df.sort_values("pearson_z", ascending=False)
    df.to_csv(TAB / "rpe1_concordance.csv", index=False)

    # ---------- the anchor tests ----------
    look = df.set_index("symbol")
    print("\n[ANCHOR · SAGA universality in RPE1] (strong result = >=2 of the available subunits universal):")
    nsaga_uni = 0
    for a in SAGA:
        if a in look.index:
            r = look.loc[a]; nsaga_uni += int(r["class"] == "universal")
            print(f"  {a:9s} pearson_z={r.pearson_z:+.4f}  class={r['class']}  well_powered={r.well_powered}")
        else:
            print(f"  {a:9s} — not perturbed in RPE1 essential-scale set")
    print(f"  => {nsaga_uni} of {sum(a in look.index for a in SAGA)} available SAGA subunits are UNIVERSAL in RPE1")

    print("\n[ANCHOR · DOCK2 T-specific in RPE1]:")
    if "DOCK2" in look.index:
        r = look.loc["DOCK2"]; print(f"  DOCK2 pearson_z={r.pearson_z:+.4f} class={r['class']}")
    else:
        print("  DOCK2 not in RPE1 essential-scale set — cannot measure directly (see class-level test).")

    # class-level compensation: are K562-T-specific regulators also non-universal in RPE1?
    try:
        k = pd.read_csv(TAB / "k562_concordance.csv")
        k_tspec = set(k.loc[k["class"] == "T-specific", "symbol"].astype(str))
        both = df[df.symbol.isin(k_tspec)]
        if len(both):
            frac_nonuni = float((both["class"] != "universal").mean())
            print(f"\n[CLASS-LEVEL T-specific robustness] K562-T-specific regulators shared with RPE1: n={len(both)} "
                  f"| mean RPE1 pearson_z={both.pearson_z.mean():+.4f} | fraction NON-universal in RPE1={frac_nonuni:.2f} "
                  f"(universal={int((both['class']=='universal').sum())} T-specific={int((both['class']=='T-specific').sum())})")
    except FileNotFoundError:
        pass

    print("\n[DONE] -> rpe1_concordance.csv, rpe1_checkpointA.txt")


if __name__ == "__main__":
    main()
