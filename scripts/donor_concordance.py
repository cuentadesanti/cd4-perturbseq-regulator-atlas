#!/usr/bin/env python3
"""Real donor-level robustness (critique #4): the honest replication unit.

The original "donor-aware" audit reweighted the EB score with donor-correlation metadata that was
present on only ~19% of contrasts (neutral-weighted, i.e. inert, for the other 81%). This uses the
dedicated per-donor DE object instead:

    s3://.../GWCD4i.DE_stats.by_donors.h5mu

whose `obs` carries `donor_correlation_hits_*` on **100%** of contrasts (3,993 KD-gated). These are
the mean / worst of the **6 pairwise donor correlations** (4 donors) of a perturbation's effect
vector on its hit genes — i.e. how well the transcriptomic response replicates across donors.

We compute, per regulator (KD-gated, peak condition):
  * `donor_corr_hits_mean` — mean of the 6 pairwise donor correlations
  * `donor_corr_hits_min`  — worst pair (strict analog of "concordant in >=3/4 donors")
  * a hard flag `donor_robust = donor_corr_hits_min >= MIN_THRESH`
and test whether donor reproducibility is (a) independent of effect size and (b) higher for the top
regulators than baseline -- and we flag large-effect hubs that FAIL to replicate across donors.

    python scripts/donor_concordance.py
Outputs: docs/tables/donor_concordance.csv (per regulator) + a printed summary.
Streams the obs only (a few MB); no 15.7 GB download.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
CACHE = Path("/Users/cuentadesanti/code/hackaton/data/cache")
TAB = ROOT / "docs" / "tables"
S3_URL = "s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.by_donors.h5mu"
MIN_THRESH = 0.5   # worst-pair donor correlation for the hard "donor_robust" flag


def _col(g, name, h5py):
    n = g[name]
    if isinstance(n, h5py.Group):  # AnnData categorical
        cats = [c.decode() if isinstance(c, bytes) else c for c in n["categories"][:]]
        codes = n["codes"][:]
        return np.array([cats[i] if i >= 0 else None for i in codes], dtype=object)
    a = n[:]
    if a.dtype.kind in ("S", "O"):
        return np.array([x.decode() if isinstance(x, bytes) else x for x in a], dtype=object)
    return a


def load_donor_obs():
    cpath = CACHE / "donor_obs.csv"
    if cpath.exists():
        return pd.read_csv(cpath)
    import h5py, fsspec
    CACHE.mkdir(parents=True, exist_ok=True)
    f = fsspec.open(S3_URL, anon=True, default_cache_type="readahead").open()
    obs = h5py.File(f, "r")["obs"]
    cols = ["target_contrast_gene_name", "culture_condition", "n_downstream",
            "ontarget_significant", "n_cells_target", "donor_correlation_hits_mean",
            "donor_correlation_hits_min", "donor_correlation_all_mean"]
    df = pd.DataFrame({c: _col(obs, c, h5py) for c in cols})
    df = df.rename(columns={"target_contrast_gene_name": "gene", "culture_condition": "condition",
                            "n_cells_target": "n_cells", "donor_correlation_hits_mean": "donor_corr_hits_mean",
                            "donor_correlation_hits_min": "donor_corr_hits_min",
                            "donor_correlation_all_mean": "donor_corr_all_mean"})
    df.to_csv(cpath, index=False)
    return df


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    d = load_donor_obs()
    sig = d[d.ontarget_significant == True].copy()
    print(f"KD-gated contrasts with donor data: {len(sig)} "
          f"(coverage {100*sig.donor_corr_hits_mean.notna().mean():.0f}%)")

    # independence of reproducibility from effect size / power (contrast level)
    rho_nd = spearmanr(sig.donor_corr_hits_mean, sig.n_downstream).correlation
    rho_nc = spearmanr(sig.donor_corr_hits_mean, sig.n_cells).correlation
    print(f"donor reproducibility vs effect size: rho(n_downstream)={rho_nd:.3f}  rho(n_cells)={rho_nc:.3f} "
          f"(near 0 -> an independent axis)")

    # per regulator: peak (max n_downstream) contrast, with its donor reproducibility
    per = (sig.sort_values("n_downstream", ascending=False)
              .groupby("gene", sort=False).first().reset_index())
    per["rank_n_downstream"] = per["n_downstream"].rank(ascending=False).astype(int)
    per["donor_robust"] = per["donor_corr_hits_min"] >= MIN_THRESH
    per = per.sort_values("n_downstream", ascending=False)
    per.to_csv(TAB / "donor_concordance.csv", index=False)

    # are the top regulators more donor-reproducible than the rest?
    top = per.head(50); rest = per.iloc[50:]
    u_p = mannwhitneyu(top.donor_corr_hits_mean, rest.donor_corr_hits_mean, alternative="greater").pvalue
    print(f"\ntop-50 donor_corr_hits_mean median={top.donor_corr_hits_mean.median():.3f}  "
          f"rest median={rest.donor_corr_hits_mean.median():.3f}  (Mann-Whitney p={u_p:.1e})")
    print(f"top-30 donor_robust (worst-pair >= {MIN_THRESH}): "
          f"{int(per.head(30).donor_robust.sum())}/30")

    # the payload: large-effect hubs that FAIL cross-donor replication
    fails = per.head(60)[~per.head(60).donor_robust].sort_values("n_downstream", ascending=False)
    print(f"\nlarge-effect hubs (top 60) that are NOT donor-robust (demote candidates):")
    print(fails[["gene", "condition", "n_downstream", "donor_corr_hits_mean",
                 "donor_corr_hits_min"]].head(12).round(3).to_string(index=False))
    print(f"\nwrote {TAB/'donor_concordance.csv'}  ({len(per)} regulators)")


if __name__ == "__main__":
    main()
