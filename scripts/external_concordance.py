#!/usr/bin/env python3
"""External ground truth (critique #6): correlate our regulator ranking against
independent references bundled by the original Marson-lab analysis repo
(emdann/GWT_perturbseq_analysis_2025), instead of only internal consistency.

Two references, two roles:
  * Schmidt & Steinhart 2022 CRISPRi screen (metadata/SchmidtSteinhart2022_...csv):
    an INDEPENDENT functional screen (cytokine phenotypes) in CD4+ T cells. Genes with a
    strong screen phenotype are known functional regulators.
  * Polarization regulator coefficients (metadata/suppl_tables/polarization_prediction_...csv):
    the PAPER'S OWN regulator-importance ranking from a different modelling task on the same cells.

We also use these as an ARBITER between our two candidate metrics: whichever of
n_downstream (thresholded count) vs mag_sig (effect magnitude over the FDR-significant set)
agrees better with an orthogonal screen is the better metric -- a decisive, external tie-break
for critique #1.

    python scripts/external_concordance.py
Outputs: docs/tables/external_concordance.csv
Reference CSVs are expected in data/external/ (download once; see REFS below).
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, hypergeom

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "docs" / "tables"
EXT = ROOT / "data" / "external"
RAW = "https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/master"
REFS = {
    "schmidt2022": f"{RAW}/metadata/SchmidtSteinhart2022_CRISPRi_screen_gene_phenotypes.csv",
    "polarization": f"{RAW}/metadata/suppl_tables/polarization_prediction_condition_comparison_regulator_coefficients.csv",
}


def ensure_refs():
    EXT.mkdir(parents=True, exist_ok=True)
    import urllib.request
    for k, url in REFS.items():
        p = EXT / f"{k}.csv"
        if not p.exists():
            print(f"  downloading {k} ...")
            urllib.request.urlretrieve(url, p)
    return {k: EXT / f"{k}.csv" for k in REFS}


def schmidt_gene_score(path):
    """Per gene: strongest functional evidence = -log10(min FDR across phenotypes/directions)."""
    df = pd.read_csv(path)
    fdr = df[["neg|fdr", "pos|fdr"]].min(axis=1).clip(lower=1e-300)
    g = pd.DataFrame({"gene": df["id"], "fdr": fdr}).groupby("gene").fdr.min()
    return (-np.log10(g)).rename("schmidt_score")


def polar_gene_score(path):
    """Per regulator: max |coef_mean| across signatures = paper's own importance."""
    df = pd.read_csv(path)
    s = df.assign(a=df["coef_mean"].abs()).groupby("regulator").a.max()
    return s.rename("polar_score")


def concordance(rank, score, name, ref_name, k=100):
    """rank: our per-gene series (higher=stronger). score: external per-gene series."""
    j = pd.concat([rank.rename("ours"), score.rename("ext")], axis=1, join="inner").dropna()
    if len(j) < 20:
        return None
    rho, p = spearmanr(j.ours, j.ext)
    # top-k overlap enrichment: our top-k vs external top-k, hypergeometric
    N = len(j)
    kk = min(k, N // 3)
    ours_top = set(j.ours.sort_values(ascending=False).head(kk).index)
    ext_top = set(j.ext.sort_values(ascending=False).head(kk).index)
    ov = len(ours_top & ext_top)
    # P(overlap >= ov) under hypergeometric
    p_enr = float(hypergeom.sf(ov - 1, N, kk, kk))
    return {"our_metric": name, "reference": ref_name, "n_shared_genes": N,
            "spearman": round(float(rho), 3), "spearman_p": f"{p:.1e}",
            "topk": kk, "topk_overlap": ov, "topk_overlap_p": f"{p_enr:.1e}"}


def main():
    refs = ensure_refs()
    per = pd.read_csv(TAB / "effect_size_ranking.csv", index_col=0)
    schmidt = schmidt_gene_score(refs["schmidt2022"])
    polar = polar_gene_score(refs["polarization"])

    our_metrics = {
        "n_downstream": per["n_downstream_peak"],
        "mag_sig": per["mag_sig_peak"],
        "mean_sig(intensity)": per["mean_sig_peak"] if "mean_sig_peak" in per else None,
        "l2_z": per["l2_z_peak"] if "l2_z_peak" in per else None,
    }
    rows = []
    for nm, r in our_metrics.items():
        if r is None:
            continue
        for ref_name, score in (("schmidt2022_screen", schmidt), ("paper_polarization_coef", polar)):
            res = concordance(r, score, nm, ref_name)
            if res:
                rows.append(res)
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "external_concordance.csv", index=False)
    print(out.to_string(index=False))
    print("\n--- arbiter (critique #1): which of our metrics agrees more with external truth ---")
    for ref_name in out.reference.unique():
        sub = out[out.reference == ref_name].sort_values("spearman", ascending=False)
        best = sub.iloc[0]
        print(f"  {ref_name}: best = {best.our_metric} (spearman={best.spearman})")
    print(f"\nwrote {TAB/'external_concordance.csv'}")


if __name__ == "__main__":
    main()
