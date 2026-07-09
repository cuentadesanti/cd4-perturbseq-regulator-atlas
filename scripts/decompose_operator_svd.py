#!/usr/bin/env python3
"""Step 1 — gene programs = right singular vectors of the operator matrix (z-score).

Recovers V (the fingerprint PCA kept only U). Orients each program by an ISG anchor,
optionally varimax-rotates for interpretability, enriches both tails offline, and —
new — runs the SAME power gate as CP: spearman(left factor u_k, per-row n_cells). A
leading program that is really a power axis gets flagged.

    python scripts/decompose_operator_svd.py --k 10 --tail-pct 2 --rotate

Outputs: docs/tables/operator_svd_programs.csv, operator_svd_enrichment.csv,
         operator_svd_power.csv, docs/figures/32_operator_svd_scree.png
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import _opkernels as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"
GENESETS = ROOT / "data" / "genesets"    # optional, gitignored; merged if present
ANCHOR = ["ISG15", "MX1", "OAS1", "IFIT1", "STAT1", "IFI6", "IRF7"]

# Inline curated gene sets (offline default; matches the COMPLEXES dict convention
# in analyze_fingerprints.py). /data/ is gitignored so these are NOT a committed file.
FALLBACK_GENESETS = {
    "IFN_ISG": ["ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT2",
                "IFIT3", "IFI6", "IRF7", "STAT1", "STAT2", "IFI44", "IFI44L",
                "RSAD2", "USP18"],
    "TCR_PROXIMAL": ["ZAP70", "LCK", "LAT", "CD3D", "CD3E", "CD3G", "CD247", "FYN",
                     "ITK", "PLCG1", "LCP2", "VAV1", "PIK3CD", "PRKCQ", "CARD11",
                     "BCL10", "MALT1"],
    "CHROMATIN_SAGA": ["TADA1", "TADA2A", "TADA2B", "TADA3", "SUPT20H", "SUPT7L",
                       "TAF5L", "TAF6L", "SGF29", "ATXN7", "ATXN7L3", "USP22",
                       "ENY2", "KAT2A", "KAT2B", "SUPT3H"],
    "MEDIATOR": ["MED1", "MED12", "MED13", "MED14", "MED23", "MED24", "CDK8",
                 "CDK19", "CCNC"],
}


def load_genesets():
    """Inline curated sets, plus any optional .gmt files under data/genesets/."""
    sets = dict(FALLBACK_GENESETS)
    if GENESETS.exists():
        for gmt in sorted(GENESETS.glob("*.gmt")):
            for line in gmt.read_text().splitlines():
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    sets[p[0]] = p[2:]
    return sets


def build_matrix(d):
    """Fully-observed fingerprint rows -> (rows, G) matrix, row index, per-row n_cells."""
    tensor, mask, n_cells = d["tensor"], d["mask"], d["n_cells"]
    R, G, C = tensor.shape
    rows, ncell = [], []
    for i in range(R):
        for c in range(C):
            if mask[i, :, c].all():
                rows.append(tensor[i, :, c]); ncell.append(n_cells[i, c])
    M = np.asarray(rows, np.float64)
    return M - M.mean(axis=0, keepdims=True), np.asarray(ncell)


def orient(v, genes):
    idx = np.isin(genes, ANCHOR)
    s = v[idx].sum() if idx.any() else v[np.argmax(np.abs(v))]
    return v if s >= 0 else -v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--tail-pct", type=float, default=2.0)
    ap.add_argument("--rotate", action="store_true")
    args = ap.parse_args()

    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    genes = d["genes"].astype(str)
    M, ncell = build_matrix(d)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    V = np.column_stack([orient(Vt[j], genes) for j in range(args.k)])   # (G,k)

    # power gate on the LEFT factors (rows == fingerprints)
    pwr = [dict(pc=j + 1, power_rho=float(spearmanr(U[:, j], ncell).statistic),
                power_confounded=bool(abs(spearmanr(U[:, j], ncell).statistic) > 0.3))
           for j in range(args.k)]

    variants = {"raw": V}
    if args.rotate:
        variants["varimax"] = op.varimax(V)[0]

    genesets, bg = load_genesets(), list(genes)
    prog_rows, enr_rows = [], []
    for rotation, VV in variants.items():
        for j in range(args.k):
            v = VV[:, j]
            hi, lo = np.percentile(v, 100 - args.tail_pct), np.percentile(v, args.tail_pct)
            for tail, sel in [("top", v >= hi), ("bottom", v <= lo)]:
                for g, val in zip(genes[sel], v[sel]):
                    prog_rows.append(dict(rotation=rotation, pc=j + 1, gene=g,
                                          loading=float(val), tail=tail))
                enr = op.hypergeometric_enrichment(list(genes[sel]), genesets, bg)
                enr.insert(0, "rotation", rotation); enr.insert(1, "pc", j + 1)
                enr.insert(2, "tail", tail); enr_rows.append(enr)

    TAB.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(prog_rows).to_csv(TAB / "operator_svd_programs.csv", index=False)
    pd.concat(enr_rows, ignore_index=True).to_csv(TAB / "operator_svd_enrichment.csv", index=False)
    pd.DataFrame(pwr).to_csv(TAB / "operator_svd_power.csv", index=False)

    FIG.mkdir(exist_ok=True, parents=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    ev = (S ** 2 / (S ** 2).sum())[:30]
    ax.plot(np.arange(1, len(ev) + 1), ev, "o-")
    ax.set_xlabel("component"); ax.set_ylabel("variance explained")
    ax.set_title("Operator matrix SVD (z-score) — scree")
    fig.tight_layout(); fig.savefig(FIG / "32_operator_svd_scree.png", dpi=150)

    e = pd.concat(enr_rows, ignore_index=True)
    clean_pcs = sorted(set(e[(e.rotation == "raw") & (e.fdr < 0.05) & (e.pc <= 5)]["pc"])
                       & set(p["pc"] for p in pwr if not p["power_confounded"]))
    print(f"top-5 PCs with FDR<0.05 label AND clean power gate: {clean_pcs}")


if __name__ == "__main__":
    main()
