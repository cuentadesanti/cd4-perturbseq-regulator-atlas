#!/usr/bin/env python3
"""Step 6 (core) — RMT-clean regulator communities.

A third, independent read of the operator: build the regulator–regulator correlation from the
fingerprint z-score, denoise it with random-matrix theory (remove the global co-variation mode
and the Marchenko–Pastur noise bulk), then find communities of co-acting regulators with a
Leiden resolution sweep wrapped in multi-seed consensus, gated by co-assignment stability >= 0.8.

Pipeline: tensor -> 3106x6000 fingerprint -> reg x reg Pearson correlation -> RMT clean
(global-mode deflation + MP-edge signal reconstruction) -> kNN graph -> Leiden sweep x seeds
-> consensus co-assignment -> stable communities.

    pip install leidenalg python-igraph
    python scripts/operator_communities.py --k 15 --res 0.2:2.0:0.2 --seeds 50 --stab-threshold 0.8

Outputs: docs/tables/operator_communities_3106.csv, docs/tables/operator_community_spectrum_3106.csv,
         docs/figures/37_operator_communities_3106.png

Validation (enrichment, RMT/community nulls, CP/SVD convergence) is a deliberate 2nd pass.
"""
import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _figstyle as fs

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"
TAB = ROOT / "docs" / "tables"
FIG = ROOT / "docs" / "figures"


def build_feature_matrix(d):
    """(3106, 2000, 3) pooled z-score tensor -> standardized 3106 x 6000 feature matrix."""
    tensor = d["tensor"].astype(np.float64)
    R, G, C = tensor.shape
    X = tensor.reshape(R, G * C)                       # concatenate the 3 conditions
    assert np.isfinite(X).all(), "tensor has non-finite entries (mask not all-observed?)"
    mu, sd = X.mean(0), X.std(0)                        # standardize each (gene x condition) column
    sd[sd == 0] = 1.0
    return (X - mu) / sd, d["regulators"].astype(str)


def rmt_clean(C, q, n_iter=8):
    """Global-mode deflation + Marchenko-Pastur denoising.

    Returns C_clean (signal-only reconstruction), a per-eigenvalue spectrum DataFrame, and the
    fitted (sigma2, lam_plus, lam0). Global mode is removed BEFORE sigma2 is estimated; the noise
    variance is the mean of the sub-lam_plus bulk (MP: mean of the bulk density == sigma2),
    iterated to convergence. Signal = eigenvalues in (lam_plus, lam0).
    """
    N = C.shape[0]
    w, V = np.linalg.eigh(C)                            # ascending
    order = np.argsort(w)[::-1]                         # -> descending
    w, V = w[order], V[:, order]
    lam0 = float(w[0])                                  # global co-variation mode
    noise = np.ones(N, bool); noise[0] = False          # global never counts as noise
    sigma2 = float(w[noise].mean())
    for _ in range(n_iter):
        lam_plus = sigma2 * (1 + np.sqrt(q)) ** 2
        new = (w <= lam_plus); new[0] = False
        if np.array_equal(new, noise):
            noise = new
            break
        noise = new
        sigma2 = float(w[noise].mean())
    lam_plus = sigma2 * (1 + np.sqrt(q)) ** 2
    signal = (w > lam_plus); signal[0] = False           # signal = above the MP edge, not global
    C_clean = (V[:, signal] * w[signal]) @ V[:, signal].T
    spec = pd.DataFrame({
        "index": np.arange(N), "eigenvalue": np.round(w, 6),
        "is_global": np.arange(N) == 0, "is_signal": signal,
        "lam_plus": round(lam_plus, 6), "sigma2": round(sigma2, 6),
    })
    return C_clean, spec, sigma2, lam_plus, lam0


def build_knn_graph(C_clean, regs, k=15):
    """Cleaned correlation -> positive-weight union-kNN igraph (dedup by max weight)."""
    import igraph as ig
    S = C_clean.copy()
    np.fill_diagonal(S, 0.0)
    n = S.shape[0]
    edge_set = {}                                       # (a,b) sorted -> max positive weight
    for i in range(n):
        nbrs = np.argsort(S[i])[::-1][:k]               # top-k by (signed) cleaned correlation
        for j in nbrs:
            wij = S[i, j]
            if wij <= 0:
                continue                                # positive (co-regulation) edges only
            a, b = (i, j) if i < j else (j, i)
            e = (int(a), int(b))
            edge_set[e] = max(edge_set.get(e, 0.0), float(wij))
    edges = list(edge_set.keys())
    weights = [edge_set[e] for e in edges]
    g = ig.Graph(n=n, edges=edges, edge_attrs={"weight": weights})
    g.vs["name"] = list(regs)
    return g


def leiden_sweep(g, gammas, seeds):
    """RBConfiguration Leiden across (resolution x seed) -> list of label vectors (n,)."""
    import leidenalg as la
    parts = []
    for gamma in gammas:
        for s in seeds:
            p = la.find_partition(g, la.RBConfigurationVertexPartition, weights="weight",
                                  resolution_parameter=float(gamma), seed=int(s), n_iterations=-1)
            parts.append(np.asarray(p.membership))
    return parts


def consensus_matrix(partitions, n):
    """Co-assignment frequency P_ij = mean_r 1[label_r(i)==label_r(j)] via one-hot Gram sums."""
    P = np.zeros((n, n), np.float64)
    for lab in partitions:
        k = int(lab.max()) + 1
        Z = np.zeros((n, k), np.float64)
        Z[np.arange(n), lab] = 1.0
        P += Z @ Z.T
    return P / len(partitions)


def final_communities(P, regs, tau=0.5, gamma=1.0, seed=0):
    """Leiden on the consensus graph (edges P_ij >= tau, weight P_ij)."""
    import igraph as ig, leidenalg as la
    n = P.shape[0]
    A = P.copy(); np.fill_diagonal(A, 0.0)
    iu = np.triu_indices(n, 1)
    m = A[iu] >= tau
    edges = list(zip(iu[0][m].tolist(), iu[1][m].tolist()))
    weights = A[iu][m].tolist()
    g = ig.Graph(n=n, edges=edges, edge_attrs={"weight": weights})
    g.vs["name"] = list(regs)
    p = la.find_partition(g, la.RBConfigurationVertexPartition, weights="weight",
                          resolution_parameter=gamma, seed=seed, n_iterations=-1)
    return np.asarray(p.membership)


def community_stability(P, labels):
    """s_c = mean intra-community co-assignment; s_i = mean co-assignment to own community."""
    n = len(labels)
    s_i = np.full(n, np.nan)
    s_c = {}
    for c in np.unique(labels):
        idx = np.where(labels == c)[0]
        if len(idx) < 2:
            s_c[c] = np.nan
            continue
        sub = P[np.ix_(idx, idx)].copy()
        np.fill_diagonal(sub, np.nan)
        s_c[c] = float(np.nanmean(sub))
        s_i[idx] = np.nanmean(sub, axis=1)
    return s_c, s_i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--res", default="0.2:2.0:0.2", help="lo:hi:step for the Leiden resolution sweep")
    ap.add_argument("--seeds", type=int, default=50)
    ap.add_argument("--stab-threshold", type=float, default=0.8)
    ap.add_argument("--consensus-tau", type=float, default=0.5)
    args = ap.parse_args()
    try:
        import igraph, leidenalg  # noqa: F401
    except ImportError:
        sys.exit("[communities] needs Leiden: pip install leidenalg python-igraph")
    fs.apply_rc()

    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    X, regs = build_feature_matrix(d)
    N, T = X.shape
    q = N / T
    assert q < 1, f"MP needs q<1; got q={q:.3f} (use concatenated conditions, not nan-mean)"
    C = np.corrcoef(X)
    print(f"[communities] N={N} regulators, T={T} features, q={q:.3f}", flush=True)

    C_clean, spec, sigma2, lam_plus, lam0 = rmt_clean(C, q)
    n_signal = int(spec["is_signal"].sum())
    print(f"[RMT] sigma2={sigma2:.3f}  lam_plus={lam_plus:.3f}  lam0(global)={lam0:.1f}  "
          f"signal modes={n_signal}", flush=True)

    lo, hi, step = (float(x) for x in args.res.split(":"))
    gammas = np.round(np.arange(lo, hi + 1e-9, step), 3)
    seeds = list(range(args.seeds))
    g = build_knn_graph(C_clean, regs, k=args.k)
    print(f"[graph] kNN k={args.k}: {g.vcount()} nodes, {g.ecount()} edges | "
          f"Leiden sweep {len(gammas)} res x {len(seeds)} seeds = {len(gammas)*len(seeds)} runs", flush=True)
    parts = leiden_sweep(g, gammas, seeds)
    P = consensus_matrix(parts, N)
    labels = final_communities(P, regs, tau=args.consensus_tau, seed=0)
    s_c, s_i = community_stability(P, labels)

    # ---- per-regulator table (+ cheap cross-ref joins) ----
    hb = pd.read_csv(TAB / "hub_ranking_bayes.csv")
    dr = set(hb.loc[hb.get("donor_robust", False) == True, "target_contrast_gene_name"].astype(str))
    cisclean = set(np.load(CACHE / "operator_tensor_cisclean.npz", allow_pickle=True)["regulators"].astype(str))
    df = pd.DataFrame({
        "regulator": regs, "community": labels,
        "community_stability": np.round([s_c[c] for c in labels], 4),
        "regulator_stability": np.round(s_i, 4),
        "is_stable": np.array([s_c[c] for c in labels]) >= args.stab_threshold,
        "donor_robust": [r in dr for r in regs],
        "in_cisclean": [r in cisclean for r in regs],
    }).sort_values(["is_stable", "community_stability", "community"], ascending=[False, False, True])
    TAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(TAB / "operator_communities_3106.csv", index=False)
    spec.to_csv(TAB / "operator_community_spectrum_3106.csv", index=False)

    stable = sorted([c for c, v in s_c.items() if np.isfinite(v) and v >= args.stab_threshold],
                    key=lambda c: -s_c[c])
    sizes = {c: int((labels == c).sum()) for c in np.unique(labels)}
    make_figure(C_clean, spec, labels, s_c, args.stab_threshold, sigma2, lam_plus, lam0, q)

    print(f"[FLAGSHIP] communities={len(np.unique(labels))} | stable(s_c>={args.stab_threshold})="
          f"{len(stable)} sizes={[sizes[c] for c in stable]} | sigma2={sigma2:.3f} lam_plus={lam_plus:.3f} "
          f"lam0={lam0:.1f} signal_modes={n_signal}", flush=True)
    print("  wrote operator_communities_3106.csv, operator_community_spectrum_3106.csv, "
          "figures/37_operator_communities_3106.png", flush=True)


def make_figure(C_clean, spec, labels, s_c, thr, sigma2, lam_plus, lam0, q):
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    # (L) MP spectrum: bulk histogram + theoretical MP density + edge + global-mode arrow
    ev = spec.eigenvalue.values
    bulk = ev[(~spec.is_global.values) & (ev <= lam_plus * 1.6)]
    a.hist(bulk, bins=80, density=True, color=fs.GENERIC, edgecolor="white", linewidth=.3, label="eigenvalues")
    lam_minus = sigma2 * (1 - np.sqrt(q)) ** 2
    xs = np.linspace(max(lam_minus, 1e-6), lam_plus, 400)
    mp = np.sqrt(np.clip((lam_plus - xs) * (xs - lam_minus), 0, None)) / (2 * np.pi * sigma2 * q * xs)
    a.plot(xs, mp, color=fs.ISG, lw=1.8, label="Marchenko–Pastur")
    a.axvline(lam_plus, color=fs.MUTED, ls="--", lw=1.2)
    ymax = a.get_ylim()[1]
    a.text(lam_plus + 0.03, ymax * 0.30, f"λ+ = {lam_plus:.2f}\n(noise edge)",
           ha="left", va="center", fontsize=8.5, color=fs.MUTED, fontweight="bold")
    n_sig = int(spec.is_signal.sum())
    a.text(0.5, 0.965, f"global mode λ₀ = {lam0:.0f}  (off-scale)\n{n_sig} signal eigenvalues above λ+",
           transform=a.transAxes, ha="center", va="top", fontsize=8.5, color=fs.INK,
           bbox=dict(boxstyle="round,pad=0.4", fc="#f6f7f9", ec="#dfe3e8"))
    a.set_xlabel("eigenvalue"); a.set_ylabel("density")
    a.set_title("RMT spectrum — global mode + MP noise removed")
    a.legend(frameon=False, fontsize=8, loc="upper right")
    # (R) cleaned correlation heatmap, block-ordered by community (stable first)
    stable = {c for c, v in s_c.items() if np.isfinite(v) and v >= thr}
    ordkey = np.array([(0 if labels[i] in stable else 1, labels[i], -C_clean[i].sum()) for i in range(len(labels))],
                      dtype=[("s", int), ("c", int), ("m", float)])
    order = np.argsort(ordkey, order=["s", "c", "m"])
    M = C_clean[np.ix_(order, order)]
    vlim = np.percentile(np.abs(M), 99)
    im = b.imshow(M, cmap="RdBu_r", vmin=-vlim, vmax=vlim, aspect="auto", interpolation="nearest")
    b.set_xticks([]); b.set_yticks([])
    b.set_title(f"cleaned correlation, community-ordered ({len(stable)} stable)")
    fig.colorbar(im, ax=b, fraction=0.046, pad=0.02, label="signal correlation")
    fs.footnote(fig, "RMT-clean regulator communities · Leiden RBConfiguration sweep + multi-seed consensus, "
                     "co-assignment stability ≥ %.1f" % thr)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG / "37_operator_communities_3106.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
