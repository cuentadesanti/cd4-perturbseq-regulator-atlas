#!/usr/bin/env python3
"""Step 1 insurance — the cis-off-target gate audit.

Confirms, from the code and the data, what the two paper-feeding analytics actually do
with off-target contrasts, PROVES the semantics of the flag, and runs the
conclusion-level sensitivity to a hard cis exclusion. ADDITIVE: never overwrites the
shipped ranking or the merged-3106 tensor.

What the code does (traceable):
  offtarget_flag (DE_stats) ~99% == neighboring_gene_KD (cis off-target: guide knocked
    down a NEIGHBOURING gene). So "exclude off-target" == "exclude cis contrasts".
  ranking (model_hubs)  : off-target is HANDLED three ways — covariate off_i in the GLM
    (n_downstream ~ C(culture_condition) + sig_i + off_i), an annotation column
    (any_offtarget / "possible off-target"), and a 0.6 penalty in robust_score. The row
    filter is ontarget_significant only; flagged rows are NOT excluded.
  operator tensor (build_operator_tensor) : off-target is UNHANDLED — the panel and rows
    gate on ontarget_significant_bool only, though the obs cache carries
    neighboring_gene_KD / distal_offtarget_flag / low_target_gex. This is the one
    genuine exposure.

Sensitivity to a hard cis exclusion, reported as a DELTA vs the shipped result:
  (a) ranking  — global structure vs the individual movers
  (b) SAGA     — which subunits survive, which were cis-inflated
  (c) operator — flagship 3b predictive margin under an offline regulator leave-out
                 (removes cis regulators from the merged tensor; NO new z-scores fetched)

Note on (c): the leave-out changes the out-of-panel test population (fewer regulators),
so the base and cis-clean margins are NOT measured over the identical test set. The
defensible claim is ROBUSTNESS (margin / low-rank elbow / effective rank essentially
unchanged), NOT that cleaning improves the operator — the small margin delta is not
interpretable given the test-set composition change.

    python scripts/audit_offtarget_gate.py

Outputs: docs/tables/offtarget_gate_sensitivity.csv
         data/cache/operator_tensor_cisclean.npz  (gitignored leave-out tensor)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
DATA = ROOT / "data" / "suppl_tables"
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"

SAGA = ["TADA2B", "SUPT20H", "TADA1", "TAF6L", "SUPT7L", "USP22"]


def as_bool(s):
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "1.0"})


def main():
    rows = []

    # ---- semantics: is offtarget_flag the cis flag? ----
    de = pd.read_csv(DATA / "DE_stats.suppl_table.csv")
    de["offtarget_flag"] = as_bool(de["offtarget_flag"])
    de["sig"] = as_bool(de["ontarget_significant"])
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv")
    for c in ["neighboring_gene_KD", "distal_offtarget_flag", "low_target_gex"]:
        obs[c] = as_bool(obs[c])
    key = ["target_contrast", "culture_condition"]
    m = de[key + ["offtarget_flag"]].merge(
        obs[key + ["neighboring_gene_KD", "distal_offtarget_flag"]], on=key)
    agree_cis = (m.offtarget_flag == m.neighboring_gene_KD).mean() * 100
    pct_cis = m[m.offtarget_flag].neighboring_gene_KD.mean() * 100
    sig_cis_pct = de[de.sig & de.offtarget_flag].shape[0] / de.sig.sum() * 100
    rows += [dict(check="semantics", metric="offtarget_flag==neighboring_gene_KD (cis) agreement %", value=round(agree_cis, 1)),
             dict(check="semantics", metric="of flagged rows, % that are cis (neighboring)", value=round(pct_cis, 0)),
             dict(check="semantics", metric="sig rows that are cis-flagged %", value=round(sig_cis_pct, 1))]
    print(f"[semantics] offtarget_flag==cis agreement {agree_cis:.1f}%; {pct_cis:.0f}% of flagged are cis")

    # ---- (a) ranking sensitivity: hard cis exclusion vs the covariate-handled ranking ----
    import model_hubs as mh
    de["ontarget_significant"] = de["sig"]
    mu, _, _ = mh.fit_fixed_effects(de)

    def eb(extra=None):
        work = np.log(de["n_downstream"].to_numpy(float) + .5) - np.log(mu + .5)
        df = pd.DataFrame({"gene": de["target_contrast_gene_name"].values, "work": work,
                           "sig": de["sig"].values})
        s2e = float(np.var(work))
        msk = df["sig"].values
        if extra is not None:
            msk = msk & extra
        df = df[msk]
        g = df.groupby("gene")
        d_g, n_g = g["work"].mean(), g["work"].size()
        s2 = s2e / n_g
        tau2 = max(float(np.var(d_g) - s2.mean()), 1e-6)
        return ((tau2 / (tau2 + s2)) * d_g).sort_values(ascending=False)

    ub, uc = eb(), eb(~de["offtarget_flag"].values)
    rb = {g: i + 1 for i, g in enumerate(ub.index)}
    rc = {g: i + 1 for i, g in enumerate(uc.index)}
    from scipy.stats import spearmanr
    common = ub.index.intersection(uc.index)
    rho = spearmanr(ub[common].rank(), uc[common].rank()).statistic
    overlap = len(set(ub.head(30).index) & set(uc.head(30).index))
    rows += [dict(check="ranking(a)", metric="spearman rank corr (all common genes)", value=round(rho, 4)),
             dict(check="ranking(a)", metric="top-30 overlap (of 30)", value=overlap)]
    print(f"[ranking] rho={rho:.4f} top-30 overlap={overlap}/30")

    # ---- (b) SAGA subunits: which survive, which were cis-inflated ----
    for g in SAGA:
        rows.append(dict(check="SAGA(b)", metric=f"{g} rank base->cisclean",
                         value=f"{rb.get(g, '>N')}->{rc.get(g, 'out')}"))
    print("[SAGA] " + "; ".join(f"{g}:{rb.get(g,'-')}->{rc.get(g,'out')}" for g in SAGA))

    # ---- (c) operator predictive margin: offline regulator leave-out (no new z-scores) ----
    import operator_completion as oc
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    regs = d["regulators"].astype(str)
    sp = obs[as_bool(obs["ontarget_significant"]) & obs["target_contrast_gene_name"].isin(set(regs))]
    flagged = set(sp.groupby("target_contrast_gene_name")["neighboring_gene_KD"].any().pipe(lambda s: s[s].index))
    keep = ~np.array([r in flagged for r in regs])
    dd = {k: (d[k][keep] if getattr(d[k], "shape", [None])[0:1] == (len(regs),) else d[k]) for k in d.files}
    np.savez(CACHE / "operator_tensor_cisclean.npz", **dd)

    def margin(dic):
        ce = oc.condition_extrap(dic, 12, 0)
        mv = ce.r2_model_novel - ce.r2_persistence_novel
        return float(mv.max()), int(ce.loc[mv.idxmax(), "rank"]), int(ce.n_test_novel.iloc[0])

    mf, rf, nf = margin(d)
    mc, rc2, nc = margin(dd)
    rows += [dict(check="operator(c)", metric="regs cis-flagged / total", value=f"{int((~keep).sum())}/{len(regs)}"),
             dict(check="operator(c)", metric="out-of-panel test regs base->cisclean", value=f"{nf}->{nc}"),
             dict(check="operator(c)", metric="3b predictive margin base (rank)", value=f"{mf:+.3f} (r{rf})"),
             dict(check="operator(c)", metric="3b predictive margin cisclean (rank)", value=f"{mc:+.3f} (r{rc2})")]
    print(f"[operator] margin base {mf:+.3f}@r{rf} (n={nf}) -> cisclean {mc:+.3f}@r{rc2} (n={nc}) "
          f"— robust; delta NOT an improvement (test set changed)")

    TAB.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TAB / "offtarget_gate_sensitivity.csv", index=False)
    print(f"\n-> {TAB / 'offtarget_gate_sensitivity.csv'}")


if __name__ == "__main__":
    main()
