#!/usr/bin/env python3
"""Fold the donor-robustness flag (Tier 1) into the main ranking and the fingerprint programs.

Policy (per review): donor_robust is added as a COLUMN + caveat, never a re-sort. Keeping a
large-effect hub visible at its rank *with* a "fails donor concordance" flag IS the finding
("rank 4, not donor-reproducible"); re-sorting would hide it.

Two propagations:
  1. Ranking tables get donor_corr_hits_{mean,min} + donor_robust columns (left join by gene).
  2. Fingerprint programs: every ASSIGNED NEIGHBOR is checked against the donor axis so the program
     table is internally consistent with the ranking. Assigned neighbors that fail the donor check
     (e.g. ATF7IP2, EIF1AX in TCR) are a program-level false positive and are flagged explicitly.

donor_robust is {True, False, NaN} — NaN means "no KD-gated donor data", which is NOT the same as
failing (absence of evidence != evidence of fragility).

    python scripts/annotate_donor_robust.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "docs" / "tables"
DON = ["donor_corr_hits_mean", "donor_corr_hits_min", "donor_robust"]

RANKING_TABLES = {   # file -> gene-key column
    "hub_ranking_bayes.csv": "target_contrast_gene_name",
    "top_robust_regulators.csv": "target_contrast_gene_name",
    "top_regulators_for_review.csv": "gene",
    "top_regulators_reproducibility_aware.csv": "gene",
}


def donor_map():
    d = pd.read_csv(TAB / "donor_concordance.csv")
    # one row per gene (peak condition already selected upstream); guard against dups
    d = d.drop_duplicates("gene", keep="first").set_index("gene")
    return d[["donor_corr_hits_mean", "donor_corr_hits_min", "donor_robust"]]


def annotate_ranking(dmap):
    for fname, key in RANKING_TABLES.items():
        p = TAB / fname
        if not p.exists():
            continue
        df = pd.read_csv(p)
        if key not in df.columns:
            print(f"  ! {fname}: no key column {key}, skipped")
            continue
        df = df.drop(columns=[c for c in DON if c in df.columns])  # idempotent
        j = df.join(dmap, on=key)  # left join, preserves row order (NO re-sort)
        for c in DON:
            df[c] = j[c].values
        df.to_csv(p, index=False)
        cov = df["donor_robust"].notna().mean()
        fragile = df[df["donor_robust"] == False]
        print(f"  {fname}: annotated {len(df)} rows (donor coverage {100*cov:.0f}%), "
              f"{len(fragile)} flagged donor-fragile")


def annotate_programs(dmap):
    # 1) per-gene program findings get the donor flag
    fp = TAB / "fingerprint_findings.csv"
    if fp.exists():
        df = pd.read_csv(fp)
        df = df.drop(columns=[c for c in ["donor_corr_hits_min", "donor_robust"] if c in df.columns])
        j = df.join(dmap[["donor_corr_hits_min", "donor_robust"]], on="gene")
        df["donor_corr_hits_min"] = j["donor_corr_hits_min"].values
        df["donor_robust"] = j["donor_robust"].values
        df.to_csv(fp, index=False)

    # 2) program evidence: flag assigned neighbors that fail the donor check
    ep = TAB / "program_label_evidence.csv"
    ev = pd.read_csv(ep)
    robust_of = dmap["donor_robust"].to_dict()
    minc_of = dmap["donor_corr_hits_min"].to_dict()

    def fragile_list(neigh):
        if not isinstance(neigh, str) or not neigh:
            return "", 0, 0
        genes = [g for g in neigh.split(";") if g]
        frag = [g for g in genes if robust_of.get(g) == False]
        return ";".join(frag), len(frag), len(genes)

    rows = ev["assigned_neighbors"].apply(fragile_list)
    ev["donor_fragile_neighbors"] = [r[0] for r in rows]
    ev["n_donor_fragile_neighbors"] = [r[1] for r in rows]
    ev["n_assigned_neighbors"] = [r[2] for r in rows]
    ev.to_csv(ep, index=False)

    print("\n  program-level donor check (assigned neighbors):")
    for r in ev.itertuples():
        if isinstance(r.assigned_neighbors, str) and r.assigned_neighbors:
            frag = r.donor_fragile_neighbors or "-"
            print(f"    {r.program_label}: {r.n_donor_fragile_neighbors}/{r.n_assigned_neighbors} "
                  f"assigned neighbors fail donor concordance -> {frag}")
    # detailed line for the flagged ones
    flagged = []
    for r in ev.itertuples():
        for g in (r.donor_fragile_neighbors.split(";") if isinstance(r.donor_fragile_neighbors, str) and r.donor_fragile_neighbors else []):
            flagged.append((g, r.program_label, minc_of.get(g)))
    if flagged:
        print("\n  program false positives (assigned but donor-fragile):")
        for g, prog, mn in flagged:
            print(f"    {g}  ({prog})  worst-pair donor rho = {mn:.3f}")


def main():
    dmap = donor_map()
    print(f"donor_concordance: {len(dmap)} genes, "
          f"{int((dmap.donor_robust == False).sum())} donor-fragile")
    print("\nranking tables:")
    annotate_ranking(dmap)
    annotate_programs(dmap)
    print("\ndone.")


if __name__ == "__main__":
    main()
