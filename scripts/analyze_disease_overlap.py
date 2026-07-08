#!/usr/bin/env python3
"""Analysis 1 — disease-signature overlap (the medical bridge), offline part.

Tests whether the convergent interferon module (docs/tables/convergent_module_genes.csv)
overlaps clinically-tracked type-I interferon signatures more than expected, against the
correct background (the 10,282 MEASURED genes, not the genome).

HONESTY: the overlap is expected/confirmatory almost by construction (the module's ISG core
is essentially the clinical signature). The reportable payload is NOT "module is disease-
relevant" but:
  (1) DIRECTION — the module edges are de-repressive (knockdown RAISES these ISGs), so the
      SAGA/Mediator coactivators RESTRAIN the clinically-tracked interferon axis → they are
      *nominated* (not proven) as upstream control points; and
  (2) SPECIFICITY of the module for interferon, shown by an exhaustion contrast set that does
      NOT overlap. See analyze_chromatin_stress_control.py: the enrichment MAGNITUDE is largely
      a general strong-perturbation effect, so this is nomination, not a claim of SAGA-specific
      disease control.

    python scripts/analyze_disease_overlap.py

Outputs: docs/tables/module_disease_overlap.csv · docs/figures/30_module_disease_overlap.png
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"

# Published, small, public signatures (curated cores — cited in the docs).
SIGNATURES = {
    # Lupus / SLE type-I IFN signature (Baechler/Kirou/Feng core).
    "Lupus type-I IFN (21-gene core)": [
        "IFI27", "IFI44", "IFI44L", "IFI6", "MX1", "OAS1", "OAS3", "OASL", "RSAD2", "SIGLEC1",
        "USP18", "HERC5", "HERC6", "LY6E", "ISG15", "IFIT1", "IFIT3", "EPSTI1", "SPATS2L",
        "PLSCR1", "LAMP3"],
    # Interferonopathy IFN score (Aicardi-Goutières / SAVI 6-7-gene panel + common extras).
    "Interferonopathy IFN score (AGS/SAVI)": [
        "IFI27", "IFI44", "IFI44L", "IFIT1", "ISG15", "RSAD2", "SIGLEC1", "USP18", "CXCL10"],
    # Canonical type-I IFN response core (broad ISG reference).
    "Type-I IFN response (canonical core)": [
        "IFI27", "IFI44", "IFI44L", "IFI6", "IFI35", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL",
        "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3", "RSAD2", "USP18", "HERC5", "HERC6",
        "STAT1", "STAT2", "IRF7", "XAF1", "BST2", "CMPK2", "HELZ2", "SAMD9", "SAMD9L", "DDX58",
        "PARP9", "PARP12", "DTX3L", "TRIM22", "GBP1", "EPSTI1", "LY6E", "PLSCR1"],
    # CONTRAST: T-cell exhaustion — should NOT be enriched (module is IFN-specific, not generically inflammatory).
    "T-cell exhaustion (contrast)": [
        "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "ENTPD1", "CD160", "BTLA", "VSIR",
        "EOMES", "NR4A1", "NR4A2", "TNFRSF9", "LAYN", "CXCL13"],
}


def main():
    var = pd.read_csv(CACHE / "fingerprint_var.csv")
    universe = set(var["gene_name"].values)
    N = len(universe)
    mod = pd.read_csv(TAB / "convergent_module_genes.csv")
    module = set(mod["measured_gene"]) & universe
    n = len(module)
    msum = json.loads((TAB / "convergent_module_summary.json").read_text())

    rows = []
    for name, genes in SIGNATURES.items():
        sig = set(genes) & universe
        K = len(sig)
        ov = module & sig
        k = len(ov)
        p = float(hypergeom.sf(k - 1, N, K, n)) if K and k else 1.0
        fold = (k / n) / (K / N) if K else 0.0
        rows.append({"signature": name, "sig_genes_total": len(genes), "sig_in_universe": K,
                     "module_size": n, "overlap": k, "fold_enrichment": round(fold, 1),
                     "hypergeom_p": p, "is_contrast": name.endswith("(contrast)"),
                     "overlap_genes": ";".join(sorted(ov))})
        print(f"  {name:42s} K={K:3d} overlap={k:2d} fold={fold:5.1f}x p={p:.1e}")
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "module_disease_overlap.csv", index=False)

    # direction payload (from the module summary)
    print(f"\n  DIRECTION: module ISG edges all positive = {msum.get('ISG_edges_all_positive')}; "
          f"module positive-share = {msum.get('module_edges_positive_share')} "
          f"→ knockdown DE-REPRESSES these ISGs; the coactivators RESTRAIN the interferon axis.")
    summary = {
        "module_size": n, "background": N,
        "overlaps": {r["signature"]: {"overlap": r["overlap"], "fold": r["fold_enrichment"],
                                      "p": r["hypergeom_p"]} for r in rows},
        "direction": "de-repressive (knockdown raises ISGs)",
        "isg_edges_all_positive": msum.get("ISG_edges_all_positive"),
        "framing": ("overlap is confirmatory by construction; the payload is the de-repressive "
                    "direction (coactivators restrain the clinical IFN axis) + upstream nomination. "
                    "The enrichment magnitude is not SAGA-specific (see chromatin_stress_control) — "
                    "nomination, not proven disease control."),
    }
    json.dump(summary, open(TAB / "module_disease_overlap_summary.json", "w"), indent=2)
    _figure(df)
    print("  figure → docs/figures/30_module_disease_overlap.png")


def _figure(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9.5, "axes.spines.top": False, "axes.spines.right": False})
    d = df.iloc[::-1]
    y = np.arange(len(d))
    colors = ["#8b95a8" if c else "#c0392b" for c in d["is_contrast"]]
    fig, ax = plt.subplots(figsize=(8.4, 4.2))
    ax.barh(y, d["fold_enrichment"], color=colors, height=0.66)
    ax.set_yticks(y); ax.set_yticklabels(d["signature"], fontsize=9)
    for yi, f, k, K, p in zip(y, d["fold_enrichment"], d["overlap"], d["sig_in_universe"], d["hypergeom_p"]):
        sig = "***" if p < 1e-10 else ("**" if p < 1e-4 else ("*" if p < 0.05 else "ns"))
        ax.text(f + max(d["fold_enrichment"]) * 0.01, yi, f"  {k}/{K}  {sig}", va="center", fontsize=8, color="0.3")
    ax.axvline(1, color="#999", ls="--", lw=1)
    ax.set_xlabel("Fold-enrichment of the 163-gene convergent module in each signature\n"
                  "(background = 10,282 measured genes)")
    ax.set_title("Module vs clinical interferon signatures — confirmatory overlap; exhaustion contrast is not enriched",
                 fontsize=10, fontweight="bold", loc="left")
    ax.set_xlim(0, max(d["fold_enrichment"]) * 1.25)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "30_module_disease_overlap.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
