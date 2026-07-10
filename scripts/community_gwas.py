#!/usr/bin/env python3
"""Task 3 — autoimmune-GWAS signal of the regulator communities, MATCHED-null only.

Asks whether the stable, unnamed communities (comm 5 n=177, comm 6 n=132; with comm 2=SAGA and
comm 7=Complex I as expected positive/negative controls) carry more autoimmune genetic-risk signal
than chance. A naive overlap vs the genome is INADMISSIBLE (immune genes are long, SNP-dense, near
the MHC). Two mandatory null layers:

  Layer 1 — universe = the 3106 panel regulators (all T-expressed/perturbable), NOT the genome.
  Layer 2 — resample MATCHED by gene-length decile, and report WITH and WITHOUT the MHC block
            (chr6 ~28-34 Mb) removed from module AND universe.

GWAS: Open Targets genetic_association scores for 5 autoimmune diseases (SLE/RA/MS/Sjogren/T1D),
GENETIC_MIN=0.01. Load statistic = sum of per-gene max genetic score over the community.

HONESTY: this is a genetics corollary / NOMINATION ("these modules touch risk genes"), NOT a claim
that a module causes disease. Directionality is not claimed.

    python scripts/community_gwas.py --nperm 10000

Outputs: docs/tables/community_gwas_matchednull_3106.csv, docs/tables/community_gwas_genes_3106.csv
Blocks (does NOT fabricate) if Open Targets or Ensembl are unreachable.
"""
import argparse, json, subprocess, time, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyze_module_gwas import disease_genetic_scores, DISEASES, GENETIC_MIN

ROOT = Path(__file__).resolve().parent.parent
TAB = ROOT / "docs" / "tables"; GS = ROOT / "data" / "genesets"
MHC = ("6", 28_000_000, 34_000_000)   # chr6 classical MHC block


def ensembl_coords(symbols):
    """{symbol: (chr, start, end)} via Ensembl REST batch lookup; cached to gene_lengths.csv."""
    cache = GS / "gene_lengths.csv"
    if cache.exists():
        c = pd.read_csv(cache, dtype={"chr": str})
        have = dict(zip(c.symbol, zip(c.chr.astype(str), c.start, c.end)))
        if set(symbols) <= set(have):
            return have
    coords = {}
    for i in range(0, len(symbols), 900):
        batch = symbols[i:i + 900]
        body = json.dumps({"symbols": batch})
        r = subprocess.run(["curl", "-s", "-m", "60", "-X", "POST",
                            "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
                            "-H", "Content-Type: application/json", "-H", "Accept: application/json",
                            "-d", body], capture_output=True, text=True)
        if r.returncode != 0 or not r.stdout.strip():
            sys.exit("[BLOCKED] Ensembl REST unreachable — cannot build the length-matched null. "
                     "Provide network or data/genesets/gene_lengths.csv (symbol,chr,start,end). "
                     "NOT fabricating lengths.")
        d = json.loads(r.stdout)
        for s, v in d.items():
            if isinstance(v, dict) and "start" in v:
                coords[s] = (str(v.get("seq_region_name")), int(v["start"]), int(v["end"]))
        print(f"  ensembl coords {min(i+900,len(symbols))}/{len(symbols)}", flush=True)
        time.sleep(0.5)
    GS.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"symbol": s, "chr": c, "start": a, "end": b} for s, (c, a, b) in coords.items()]
                 ).to_csv(cache, index=False)
    return coords


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0); args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    comm = pd.read_csv(TAB / "operator_communities_3106.csv")
    regs = comm.regulator.astype(str).tolist()

    # ---- GWAS scores (5 diseases) -> per-regulator max score + per-disease + diseases hit ----
    dz_scores = {}
    for name, efo in DISEASES.items():
        dz_scores[name] = disease_genetic_scores(efo)
        print(f"  OT {name}: {len(dz_scores[name])} associated targets", flush=True)
    per = {r: {dz: dz_scores[dz].get(r, 0.0) for dz in DISEASES} for r in regs}
    max_score = {r: max(per[r].values()) for r in regs}
    diseases_hit = {r: [dz for dz in DISEASES if per[r][dz] >= GENETIC_MIN] for r in regs}

    # per-gene detail table
    gt = pd.DataFrame([dict(community=comm.set_index("regulator").loc[r, "community"], gene=r,
                            max_genetic_score=round(max_score[r], 3),
                            diseases_hit=";".join(diseases_hit[r]),
                            **{f"score_{dz}": round(per[r][dz], 3) for dz in DISEASES})
                       for r in regs if max_score[r] >= GENETIC_MIN])
    TAB.mkdir(parents=True, exist_ok=True)
    gt.sort_values(["community", "max_genetic_score"], ascending=[True, False]).to_csv(
        TAB / "community_gwas_genes_3106.csv", index=False)

    # ---- gene lengths + MHC flag ----
    coords = ensembl_coords(regs)
    length = np.array([coords[r][2] - coords[r][1] + 1 if r in coords else np.nan for r in regs], float)
    in_mhc = np.array([r in coords and coords[r][0] == MHC[0] and coords[r][1] < MHC[2] and coords[r][2] > MHC[1]
                       for r in regs])
    score = np.array([max_score[r] for r in regs])
    print(f"  coords resolved {int(np.isfinite(length).sum())}/{len(regs)} | MHC genes in panel: {int(in_mhc.sum())}", flush=True)

    regs = np.array(regs)
    def load(mask):                     # community load = sum of per-gene max genetic score
        return float(score[mask].sum())

    def emp_p(obs, comm_mask, universe_mask, matched):
        """empirical p of the community load vs `nperm` resamples of the same size from universe."""
        uni = np.where(universe_mask)[0]; n = int(comm_mask.sum())
        if matched:
            # length-decile-stratified resample: match the community's per-decile counts
            fin = np.isfinite(length) & universe_mask
            dec = np.full(len(regs), -1)
            valid = np.where(fin)[0]
            dec[valid] = pd.qcut(length[valid], 10, labels=False, duplicates="drop")
            comm_dec = dec[comm_mask & fin]
            need = pd.Series(comm_dec).value_counts().to_dict()
            pools = {d: valid[dec[valid] == d] for d in need}
            if any(len(pools[d]) < need[d] for d in need):
                return np.nan, np.nan
            null = np.empty(args.nperm)
            for b in range(args.nperm):
                pick = np.concatenate([rng.choice(pools[d], need[d], replace=False) for d in need])
                null[b] = score[pick].sum()
        else:
            null = np.array([score[rng.choice(uni, n, replace=False)].sum() for _ in range(args.nperm)])
        return float((null >= obs).mean()), float(null.mean())

    rows = []
    for c in sorted(comm.community.unique()):
        cm = (comm.community.values == c)
        n = int(cm.sum()); obs = load(cm)
        # Layer 1 (uniform panel null)
        p_uni, mu_uni = emp_p(obs, cm, np.ones(len(regs), bool), matched=False)
        # Layer 2 (length-matched)
        p_mat, mu_mat = emp_p(obs, cm, np.ones(len(regs), bool), matched=True)
        # Layer 2 + MHC excluded (from module AND universe)
        keep = ~in_mhc
        obs_x = load(cm & keep)
        p_mhc, mu_mhc = emp_p(obs_x, cm & keep, keep, matched=True)
        # top disease
        cg = gt[gt.community == c]
        top_dz = (cg[[f"score_{dz}" for dz in DISEASES]].sum().idxmax().replace("score_", "")
                  if len(cg) else "none")
        n_hits = int((score[cm] >= GENETIC_MIN).sum())
        rows.append(dict(community=c, n=n, is_stable=bool(comm[comm.community == c].is_stable.iloc[0]),
                         n_gwas_hits=n_hits, gwas_load_obs=round(obs, 3),
                         null_mean_matched=round(mu_mat, 3) if np.isfinite(mu_mat) else np.nan,
                         fold=round(obs / mu_mat, 2) if mu_mat else np.nan,
                         p_emp_uniform=round(p_uni, 4), p_emp=round(p_mat, 4) if np.isfinite(p_mat) else np.nan,
                         p_emp_MHCexcluded=round(p_mhc, 4) if np.isfinite(p_mhc) else np.nan,
                         top_disease=top_dz))
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "community_gwas_matchednull_3106.csv", index=False)
    print("\n[Task 3] community × autoimmune-GWAS (length-matched null, MHC-robust):")
    print(out[["community", "n", "is_stable", "n_gwas_hits", "fold", "p_emp_uniform", "p_emp",
               "p_emp_MHCexcluded", "top_disease"]].to_string(index=False), flush=True)
    for c in (5, 6):
        r = out[out.community == c].iloc[0]
        verdict = ("POSITIVE" if (r.p_emp < 0.05 and r.p_emp_MHCexcluded < 0.05) else "negative")
        print(f"  comm {c} (stable, unnamed): {verdict} (p_matched={r.p_emp}, p_MHCexcl={r.p_emp_MHCexcluded})", flush=True)


if __name__ == "__main__":
    main()
