#!/usr/bin/env python3
"""Step 6, 2nd pass · Task B — the SAGA module across cell types (K562).

The unit is the SAGA SUBUNIT SET, not the whole (mostly non-SAGA) community 2. Two tests:
 (1) CD4<->K562 concordance per subunit (from docs/tables/k562_concordance.csv): how many of
     the SAGA subunits are universal in K562 — the n>1 result replaces the single SUPT20H anchor.
 (2) Non-circular module test: build the regulator-regulator correlation IN K562 (Replogle bulk,
     same pipeline as CD4) and ask whether the SAGA subunits form a tighter cluster (mean intra
     correlation) than random same-size regulator sets -> empirical p / z.

    python scripts/k562_saga_module.py

Outputs: docs/tables/k562_saga_module_3106.csv (+ a printed module-modularity statement)
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_k562_concordance as k562

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"
SAGA = ["SUPT20H", "TAF6L", "TADA2B", "SUPT7L", "USP22", "SGF29", "KAT2B"]


def main():
    rng = np.random.default_rng(0)
    # ---- (1) per-subunit K562 concordance (from the shipped table) ----
    kc = pd.read_csv(TAB / "k562_concordance.csv")
    sub = kc[kc.symbol.isin(SAGA)][["symbol", "pearson_z", "class", "well_powered", "donor_robust"]].copy()
    present = set(sub.symbol)
    rows = [dict(subunit=s, en_k562_shared=s in present,
                 pearson_z=float(sub.loc[sub.symbol == s, "pearson_z"].iloc[0]) if s in present else np.nan,
                 klass=sub.loc[sub.symbol == s, "class"].iloc[0] if s in present else "not_in_shared",
                 well_powered=bool(sub.loc[sub.symbol == s, "well_powered"].iloc[0]) if s in present else False,
                 donor_robust=bool(sub.loc[sub.symbol == s, "donor_robust"].iloc[0]) if s in present else False)
            for s in SAGA]
    df = pd.DataFrame(rows)
    n_univ = int((df.klass == "universal").sum())
    print(f"[B1] SAGA subunits in K562 shared set: {int(df.en_k562_shared.sum())}/{len(SAGA)} | "
          f"universal: {n_univ} | mean pearson_z: {df.pearson_z.mean():.3f}", flush=True)

    # ---- (2) non-circular: SAGA module cohesion in the K562 reg-reg correlation ----
    Keff, Kreg, vg, kcells, n_ctrl = k562.build_k562()             # (regulators x 8248 genes), z-normalized
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv")[["target_contrast", "target_contrast_gene_name"]].drop_duplicates()
    s2e = dict(zip(obs.target_contrast_gene_name.astype(str), obs.target_contrast.astype(str)))
    kri = {r: i for i, r in enumerate(Kreg)}
    saga_ens = [s2e.get(s) for s in SAGA if s2e.get(s) in kri]
    idx = [kri[e] for e in saga_ens]
    Z = Keff                                                      # rows = regulators, cols = genes
    finite = np.isfinite(Z).all(1)                                # keep fully-measured regulators for the null pool
    Ck = np.corrcoef(Z[idx])                                      # SAGA x SAGA correlation in K562
    iu = np.triu_indices(len(idx), 1)
    obs_mean = float(np.nanmean(Ck[iu]))
    pool = np.where(finite)[0]
    B = 2000; nullmeans = np.empty(B)
    for b in range(B):
        pick = rng.choice(pool, len(idx), replace=False)
        Cn = np.corrcoef(Z[pick]); nullmeans[b] = np.nanmean(Cn[np.triu_indices(len(idx), 1)])
    p = float((nullmeans >= obs_mean).mean())
    z = float((obs_mean - nullmeans.mean()) / (nullmeans.std(ddof=1) + 1e-12))
    print(f"[B2] SAGA module cohesion in K562: mean intra-corr={obs_mean:.3f} "
          f"vs random same-size {nullmeans.mean():.3f}±{nullmeans.std(ddof=1):.3f} | z={z:.1f} p={p:.4f} "
          f"(n_subunits_in_K562={len(idx)})", flush=True)

    df["k562_module_intra_corr"] = round(obs_mean, 4)
    df["k562_module_null_mean"] = round(float(nullmeans.mean()), 4)
    df["k562_module_z"] = round(z, 2)
    df["k562_module_p"] = round(p, 4)
    TAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(TAB / "k562_saga_module_3106.csv", index=False)
    print(f"[B DONE] universal={n_univ}/{int(df.en_k562_shared.sum())} shared; module z={z:.1f} p={p:.4f} "
          f"-> k562_saga_module_3106.csv", flush=True)


if __name__ == "__main__":
    main()
