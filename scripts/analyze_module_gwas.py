#!/usr/bin/env python3
"""Analysis 1 (GWAS layer) — autoimmune genetic overlap via Open Targets (needs network).

For each gene in the convergent interferon module, ask whether it is also a genetic-association
(GWAS) risk gene for an autoimmune disease (SLE, RA, MS, Sjogren's, T1D) in Open Targets. If a
module gene is a disease risk gene AND the SAGA/Mediator coactivators restrain it (de-repressive
direction, from Analysis 1), that is a mechanistic hypothesis. Also checks the SAGA regulators
themselves for reported autoimmune genetic association.

HONESTY: this is a genetics corollary / nomination, not a claim of novel disease association.
Overlap of an interferon module with autoimmune genetics is expected.

    python scripts/analyze_module_gwas.py

Outputs: docs/tables/module_gwas_hits.csv · docs/tables/module_gwas_summary.json
"""
import json
import subprocess
import time
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
OT = "https://api.platform.opentargets.org/api/v4/graphql"

DISEASES = {"SLE": "MONDO_0007915", "RA": "MONDO_0008383", "MS": "MONDO_0005301",
            "Sjogren": "MONDO_0010030", "T1D": "MONDO_0005147"}
SAGA_REGS = ["SGF29", "TADA1", "TADA2B", "SUPT20H", "TAF6L", "USP22", "SUPT7L", "ATXN7L3"]
GENETIC_MIN = 0.01   # minimum genetic_association datatype score to count as a GWAS/genetic hit

Q = """query($efo:String!,$idx:Int!){ disease(efoId:$efo){
  associatedTargets(page:{index:$idx,size:500}){ count
    rows{ target{ approvedSymbol } score datatypeScores{ id score } } } } }"""


def gql(query, variables, tries=3):
    """POST a GraphQL query via curl (system CA certs work; Python's bundled ones don't here)."""
    body = json.dumps({"query": query, "variables": variables})
    for t in range(tries):
        r = subprocess.run(["curl", "-s", "-m", "40", "-X", "POST", OT,
                            "-H", "Content-Type: application/json", "-d", body],
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip():
            try:
                return json.loads(r.stdout)
            except Exception:
                pass
        if t == tries - 1:
            raise RuntimeError(f"Open Targets query failed: {r.stderr[:200] or r.stdout[:200]}")
        time.sleep(2)


def disease_genetic_scores(efo):
    """{symbol: genetic_association score} for all associated targets of a disease."""
    out = {}
    idx = 0
    while True:
        d = gql(Q, {"efo": efo, "idx": idx})["data"]["disease"]["associatedTargets"]
        for row in d["rows"]:
            sym = row["target"]["approvedSymbol"]
            gscore = next((s["score"] for s in row["datatypeScores"] if s["id"] == "genetic_association"), 0.0)
            if gscore > 0:
                out[sym] = max(out.get(sym, 0.0), round(gscore, 3))
        idx += 1
        if idx * 500 >= d["count"] or idx >= 8:   # cap at 4000 targets/disease
            break
    return out


def main():
    mod = pd.read_csv(TAB / "convergent_module_genes.csv")
    module = dict(zip(mod["measured_gene"], mod["is_ISG"]))

    print("Querying Open Targets genetic associations for 5 autoimmune diseases…")
    dz_scores = {}
    for name, efo in DISEASES.items():
        s = disease_genetic_scores(efo)
        dz_scores[name] = s
        print(f"  {name:8s} ({efo}): {len(s)} genetically-associated targets")

    # module gene × disease genetic hits
    rows = []
    for g, isg in module.items():
        hits = {name: dz_scores[name].get(g, 0.0) for name in DISEASES}
        n_dz = sum(1 for v in hits.values() if v >= GENETIC_MIN)
        if n_dz:
            rows.append({"gene": g, "is_ISG": bool(isg), "n_autoimmune_diseases": n_dz,
                         **{f"gwas_{k}": hits[k] for k in DISEASES},
                         "max_genetic_score": round(max(hits.values()), 3)})
    hits_df = pd.DataFrame(rows).sort_values(["n_autoimmune_diseases", "max_genetic_score"], ascending=False)
    hits_df.to_csv(TAB / "module_gwas_hits.csv", index=False)

    # SAGA regulators themselves
    saga_rows = []
    for g in SAGA_REGS:
        hits = {name: dz_scores[name].get(g, 0.0) for name in DISEASES}
        saga_rows.append({"regulator": g, **{f"gwas_{k}": hits[k] for k in DISEASES},
                          "max_genetic_score": round(max(hits.values()), 3)})
    saga_df = pd.DataFrame(saga_rows)

    n_module = len(module)
    n_hit = len(hits_df)
    n_isg_hit = int(hits_df["is_ISG"].sum()) if n_hit else 0
    print(f"\n  module genes that are autoimmune genetic (GWAS) risk genes: {n_hit}/{n_module} "
          f"({n_isg_hit} of them ISGs)")
    if n_hit:
        print(hits_df.head(12)[["gene", "is_ISG", "n_autoimmune_diseases", "max_genetic_score"]].to_string(index=False))
    print("\n  SAGA regulators' own autoimmune genetic association (max score across the 5 diseases):")
    print(saga_df[["regulator", "max_genetic_score"]].to_string(index=False))

    summary = {
        "module_size": n_module, "diseases": DISEASES,
        "genetic_score_min": GENETIC_MIN,
        "module_genes_with_autoimmune_gwas": n_hit,
        "of_which_ISGs": n_isg_hit,
        "top_hits": hits_df.head(15)["gene"].tolist() if n_hit else [],
        "saga_regulators_max_score": {r["regulator"]: r["max_genetic_score"] for r in saga_rows},
        "framing": ("Genetics corollary / nomination, not a novel disease association. An interferon "
                    "module overlapping autoimmune GWAS genes is expected; the point is that these are "
                    "restrained by the coactivators (de-repressive direction), giving a mechanistic "
                    "hypothesis. Enrichment magnitude is not SAGA-specific (see chromatin_stress_control)."),
    }
    json.dump(summary, open(TAB / "module_gwas_summary.json", "w"), indent=2)
    _figure(hits_df, saga_df)
    print("\n  tables → module_gwas_hits.csv · module_gwas_summary.json · figure → 31_module_gwas_autoimmune.png")


def _figure(hits_df, saga_df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from pathlib import Path as _P
    FIG = ROOT / "docs" / "figures"
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    top = hits_df.head(16).iloc[::-1]
    y = np.arange(len(top))
    colors = ["#c0392b" if isg else "#8e44ad" for isg in top["is_ISG"]]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [2.1, 1]})
    ax1.barh(y, top["max_genetic_score"], color=colors, height=0.68)
    ax1.set_yticks(y); ax1.set_yticklabels(top["gene"], fontsize=8.5)
    for yi, s, n in zip(y, top["max_genetic_score"], top["n_autoimmune_diseases"]):
        ax1.text(s + 0.01, yi, f"  {n} dz", va="center", fontsize=7.5, color="0.35")
    ax1.set_xlabel("max Open Targets genetic-association score (autoimmune)")
    ax1.set_title(f"{len(hits_df)}/163 module genes are autoimmune GWAS risk genes",
                  fontsize=10, fontweight="bold", loc="left")
    ax1.legend(handles=[Patch(color="#c0392b", label="ISG"), Patch(color="#8e44ad", label="non-ISG")],
               fontsize=8, loc="lower right", frameon=False)
    ax1.set_xlim(0, 1.0)

    sd = saga_df.sort_values("max_genetic_score").set_index("regulator")["max_genetic_score"]
    yy = np.arange(len(sd))
    ax2.barh(yy, sd.values, color="#2980b9", height=0.6)
    ax2.set_yticks(yy); ax2.set_yticklabels(sd.index, fontsize=8.5)
    ax2.set_xlabel("genetic-association score")
    ax2.set_title("SAGA regulators' own\nautoimmune genetics", fontsize=9.5, fontweight="bold", loc="left")
    ax2.set_xlim(0, max(0.5, sd.max() * 1.15))

    fig.suptitle("Genetics corollary — module genes restrained by SAGA/Mediator overlap autoimmune GWAS "
                 "(e.g. STAT4); nomination, not novel association", fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _P(FIG).mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "31_module_gwas_autoimmune.png", dpi=200, bbox_inches="tight")
    plt.close(fig)



if __name__ == "__main__":
    main()
