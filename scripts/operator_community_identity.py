#!/usr/bin/env python3
"""Step 6, 2nd pass · Tasks C+D — community identity (enrichment) + CP convergence.

Task C — hypergeometric enrichment (BH-FDR) of each community vs known complexes, giving the
stable communities a biological identity. Curated regulator complexes (SAGA/Mediator/TCR) are
always tested; pass --corum PATH to allComplexes.txt to add the broad CORUM database (the CORUM
5.x site is a JS SPA with no working direct download, so this is a manual file until provided —
we do NOT fabricate memberships).
Task D — cross-tab each CP factor's top regulators to communities (consistency cross-check).

    python scripts/operator_community_identity.py [--corum data/genesets/allComplexes.txt]

Outputs: docs/tables/operator_community_enrichment_3106.csv, operator_community_cp_overlap_3106.csv
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import _opkernels as op
from analyze_fingerprints import COMPLEXES

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "docs" / "tables"
TARGETS = [2, 5, 6, 7]   # SAGA community + the 3 stable (5/6/7)


def load_corum(path, background):
    """Parse CORUM allComplexes.txt -> {complex_name: [gene symbols]} for human complexes."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    org = next((c for c in df.columns if "organism" in c.lower()), None)
    if org is not None:
        df = df[df[org].str.contains("Human", case=False, na=False)]
    namec = next(c for c in df.columns if "complex" in c.lower() and "name" in c.lower())
    subc = next(c for c in df.columns if "subunits" in c.lower() and "gene" in c.lower())
    sets = {}
    for _, r in df.iterrows():
        members = [g.strip() for g in str(r[subc]).replace(";", " ").split() if g.strip()]
        if 2 <= len(set(members) & set(background)):        # only complexes with >=2 members present
            sets[str(r[namec])[:60]] = members
    return sets


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--corum", default=None); args = ap.parse_args()
    comm = pd.read_csv(TAB / "operator_communities_3106.csv")
    bg = comm.regulator.tolist()
    gene_sets = dict(COMPLEXES); source = "curated(SAGA/Mediator/TCR)"
    if args.corum and Path(args.corum).exists():
        gene_sets.update(load_corum(args.corum, bg)); source = "curated+CORUM"
        print(f"[C] CORUM loaded: {len(gene_sets)} complex sets total", flush=True)
    else:
        print("[C] CORUM not provided -> curated sets only (broad identity BLOCKED; "
              "pass --corum allComplexes.txt)", flush=True)

    # ---- Task C: enrichment per community ----
    rows = []
    for c in TARGETS:
        regs = comm[comm.community == c].regulator.tolist()
        e = op.hypergeometric_enrichment(regs, gene_sets, bg)
        e = e[e.set_size > 0].copy()
        e.insert(0, "community", c); e.insert(1, "stable", bool(comm[comm.community == c].is_stable.iloc[0]))
        rows.append(e)
    enr = pd.concat(rows, ignore_index=True)
    enr["source"] = source
    enr = enr[["community", "stable", "set_name", "n_overlap", "set_size", "pvalue", "fdr", "overlap_genes", "source"]]
    enr.to_csv(TAB / "operator_community_enrichment_3106.csv", index=False)
    sig = enr[enr.fdr < 0.05]
    print(f"[C] enrichment: {len(sig)} community×complex hits at FDR<0.05", flush=True)
    for _, r in sig.iterrows():
        print(f"    community {r.community} (stable={r.stable}): {r.set_name}  FDR={r.fdr:.1e}  "
              f"({r.n_overlap}/{r.set_size})", flush=True)

    # ---- Task D: CP factor -> community overlap ----
    c2c = dict(zip(comm.regulator, comm.community)); stable = set(comm[comm.is_stable].community)
    cp = pd.read_csv(TAB / "operator_cp_factors_3106.csv"); drows = []
    for _, r in cp.iterrows():
        tops = [g for g in str(r.top_regulators).split(";") if g in c2c]
        if not tops:
            continue
        vc = pd.Series([c2c[g] for g in tops]).value_counts()
        drows.append(dict(factor=int(r.factor), gating=r.gating_shape, n_top_in_operator=len(tops),
                          dominant_community=int(vc.index[0]), dom_is_stable=int(vc.index[0]) in stable,
                          concentration=round(vc.iloc[0] / len(tops), 2), community_spread=len(vc),
                          top_regs=";".join(tops)))
    pd.DataFrame(drows).to_csv(TAB / "operator_community_cp_overlap_3106.csv", index=False)
    print(f"[D] CP->community overlap: {len(drows)} factors -> operator_community_cp_overlap_3106.csv", flush=True)


if __name__ == "__main__":
    main()
