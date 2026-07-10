#!/usr/bin/env python3
"""Step 6, 2nd pass · Task A — empirical noise-edge + modularity null (the guardrail).

The MP edge lam+ (1.49) is a LOWER bound because the 6000 feature columns are not i.i.d.
(same gene x 3 conditions), so the "336 signal eigenvalues" over-counts. This fixes the edge
EMPIRICALLY: column-permute the real fingerprint (destroys reg-reg correlation, preserves each
column's marginal), recompute the correlation, take its top eigenvalue -> the empirical noise
edge; 95th percentile over >=100 nulls is the honest cutoff. Also: does the consensus
partition's modularity exceed a degree-preserving rewired-graph null (z-score)?

    python scripts/operator_community_null.py --nperm 100

Output: docs/tables/operator_community_null_3106.csv
"""
import argparse, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy.sparse.linalg import eigsh
import operator_communities as oc

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"


def top_eig(C):
    return float(eigsh(C, k=1, which="LA", return_eigenvectors=False)[0])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--nperm", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0); args = ap.parse_args()
    t0 = time.time(); rng = np.random.default_rng(args.seed)

    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    X, regs = oc.build_feature_matrix(d)
    N, Tt = X.shape; q = N / Tt
    C = np.corrcoef(X)
    C_clean, spec, sigma2, lam_plus, lam0 = oc.rmt_clean(C, q)
    ev = spec.eigenvalue.values
    n_sig_theory = int(((ev > lam_plus) & (~spec.is_global.values)).sum())
    print(f"[A] theoretical lam+={lam_plus:.3f}  n_signal(theory)={n_sig_theory}", flush=True)

    # --- empirical noise edge: top eigenvalue of column-permuted null correlations ---
    ncol = np.arange(Tt)
    edges = np.empty(args.nperm)
    for b in range(args.nperm):
        perm = np.argsort(rng.random((N, Tt)), axis=0)      # independent per-column row shuffle
        Xn = X[perm, ncol]
        edges[b] = top_eig(np.corrcoef(Xn))
        if (b + 1) % 20 == 0:
            print(f"  null {b+1}/{args.nperm}  edge~{edges[b]:.3f}  ({time.time()-t0:.0f}s)", flush=True)
    lam_emp = float(np.quantile(edges, 0.95))
    n_sig_emp = int(((ev > lam_emp) & (~spec.is_global.values)).sum())
    print(f"[A] empirical lam+(95pct)={lam_emp:.3f} [{edges.min():.3f}-{edges.max():.3f}]  "
          f"n_signal(empirical)={n_sig_emp}", flush=True)

    # --- modularity null: consensus partition vs degree-preserving rewired graphs ---
    import igraph as ig
    g = oc.build_knn_graph(C_clean, regs, k=15)
    lab = pd.read_csv(TAB / "operator_communities_3106.csv").set_index("regulator").reindex(regs)["community"].to_numpy()
    lab = np.asarray(lab, int)
    Q_real = g.modularity(lab, weights="weight")
    Q_null = np.empty(min(args.nperm, 200))                          # label-permutation null
    for b in range(len(Q_null)):                                     # (same community sizes, shuffled)
        Q_null[b] = g.modularity(rng.permutation(lab), weights="weight")
    z_mod = (Q_real - Q_null.mean()) / (Q_null.std(ddof=1) + 1e-12)
    print(f"[A] modularity: real Q={Q_real:.3f}  null={Q_null.mean():.3f}±{Q_null.std(ddof=1):.3f}  z={z_mod:.1f}", flush=True)

    pd.DataFrame([
        dict(metric="noise_edge_theoretical_MP", value=round(lam_plus, 4), n_signal=n_sig_theory),
        dict(metric="noise_edge_empirical_p95", value=round(lam_emp, 4), n_signal=n_sig_emp),
        dict(metric="modularity_real", value=round(Q_real, 4), n_signal=np.nan),
        dict(metric="modularity_null_mean", value=round(float(Q_null.mean()), 4), n_signal=np.nan),
        dict(metric="modularity_null_z", value=round(float(z_mod), 2), n_signal=np.nan),
    ]).to_csv(TAB / "operator_community_null_3106.csv", index=False)
    print(f"[A DONE] -> operator_community_null_3106.csv  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
