#!/usr/bin/env python3
"""Step 3 — is the z-score operator recoverably low-rank? Held-out prediction.

3b (FLAGSHIP): hold out (regulator, Stim48hr) fibers, PREFERENTIALLY for out-of-panel
regulators, and predict them from that regulator's Rest+Stim8hr via the low-rank fit
on the other regulators. Baseline = persistence (Stim48hr := Stim8hr).
3a (SANITY): random entry-wise completion; report the effective-rank elbow and R²
vs a rank-1 baseline (beating the gene-mean baseline is near-trivial and not claimed).

    python scripts/operator_completion.py --max-rank 12 --holdout 0.2

Outputs: docs/tables/operator_completion_condition.csv,
         operator_completion_entrywise.csv, docs/figures/35_*
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _opkernels as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"; FIG = ROOT / "docs" / "figures"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def to_matrix(d):
    tensor, mask = d["tensor"], d["mask"]
    R, G, C = tensor.shape
    rows = [tensor[i, :, c] for i in range(R) for c in range(C) if mask[i, :, c].all()]
    return np.asarray(rows, float)


def entrywise(M, max_rank, holdout, seed):
    rng = np.random.default_rng(seed)
    hide = rng.random(M.shape) < holdout
    train = ~hide
    Mc, _ = op.train_test_standardize(M, train)
    yt = Mc[hide]
    base_gm = float(np.mean(yt ** 2))                        # gene-mean == 0 (centered)
    r1 = op.soft_impute(Mc, train, 1, n_iter=200)
    base_r1 = float(np.mean((yt - r1[hide]) ** 2))
    rows = []
    for r in range(1, max_rank + 1):
        print(f"[3a] entry-wise rank {r}/{max_rank}", flush=True)
        yp = op.soft_impute(Mc, train, r, n_iter=200)[hide]
        mse = float(np.mean((yt - yp) ** 2))
        rows.append(dict(rank=r, r2_vs_geneMean=1 - mse / (base_gm + 1e-12),
                         r2_vs_rank1=1 - mse / (base_r1 + 1e-12),
                         beats_rank1=bool(mse < base_r1)))
    return pd.DataFrame(rows)


def condition_extrap(d, max_rank, seed):
    tensor, mask, in_orig = d["tensor"], d["mask"], d["in_original_panel"]
    R, G, C = tensor.shape
    late, early = COND_ORDER.index("Stim48hr"), COND_ORDER.index("Stim8hr")
    full = np.array([i for i in range(R) if mask[i, :, :].all(axis=0).all()])
    novel_full = full[~in_orig[full]]
    rng = np.random.default_rng(seed)
    # hold out preferentially from out-of-panel regulators
    n_test = max(20, len(full) // 5)
    test = rng.permutation(novel_full)[:n_test]
    if len(test) < n_test:                                   # top up from any full reg
        extra = rng.permutation(np.setdiff1d(full, test))[: n_test - len(test)]
        test = np.concatenate([test, extra])
    X = np.concatenate([tensor[full, :, 0], tensor[full, :, early], tensor[full, :, late]], axis=1)
    obs = np.ones_like(X, bool)
    test_pos = np.isin(full, test)
    novel_pos = test_pos & (~in_orig[full])
    obs[test_pos, 2 * G:3 * G] = False
    Xc, _ = op.train_test_standardize(X, obs)
    rows = []
    for r in range(1, max_rank + 1):
        print(f"[3b] condition rank {r}/{max_rank}", flush=True)
        Xhat = op.soft_impute(Xc, obs, r, n_iter=200)
        def block(posmask):
            yt = Xc[np.ix_(posmask, np.arange(2 * G, 3 * G))]
            yp = Xhat[np.ix_(posmask, np.arange(2 * G, 3 * G))]
            pers = Xc[np.ix_(posmask, np.arange(G, 2 * G))]
            ss = np.sum((yt - yt.mean()) ** 2) + 1e-12
            return (float(1 - np.sum((yt - yp) ** 2) / ss),
                    float(1 - np.sum((yt - pers) ** 2) / ss))
        rm, rp = block(test_pos)
        rmn, rpn = block(novel_pos)
        rows.append(dict(rank=r, r2_model=rm, r2_persistence=rp, beats_persistence=bool(rm > rp),
                         r2_model_novel=rmn, r2_persistence_novel=rpn,
                         beats_persistence_novel=bool(rmn > rpn),
                         n_test=int(test_pos.sum()), n_test_novel=int(novel_pos.sum())))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-rank", type=int, default=12)
    ap.add_argument("--holdout", type=float, default=0.2)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    print("[note] pure low-rank cannot predict a regulator with ZERO observed entries; "
          "3b holds out a condition FIBER, never a whole regulator. No genome-scale claim.")
    ew = entrywise(to_matrix(d), args.max_rank, args.holdout, args.seed)
    ce = condition_extrap(d, args.max_rank, args.seed)
    TAB.mkdir(exist_ok=True, parents=True)
    ew.to_csv(TAB / "operator_completion_entrywise.csv", index=False)
    ce.to_csv(TAB / "operator_completion_condition.csv", index=False)

    FIG.mkdir(exist_ok=True, parents=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ce["rank"], ce["r2_model_novel"], "s-", label="3b model R² (out-of-panel)")
    ax.plot(ce["rank"], ce["r2_persistence_novel"], "--", c="crimson", label="3b persistence (out-of-panel)")
    ax.plot(ew["rank"], ew["r2_vs_rank1"], "o-", alpha=0.5, label="3a R² vs rank-1 (sanity)")
    ax.axhline(0, ls=":", c="grey"); ax.set_xlabel("rank"); ax.set_ylabel("held-out R²")
    ax.legend(fontsize=8); ax.set_title("Operator low-rank recoverability (flagship = 3b out-of-panel)")
    fig.tight_layout(); fig.savefig(FIG / "35_operator_completion_curve.png", dpi=150)

    print(f"[3b FLAGSHIP] out-of-panel n={int(ce.n_test_novel.iloc[0])}; "
          f"beats persistence at ranks={ce.loc[ce.beats_persistence_novel,'rank'].tolist()}; "
          f"max model R²={ce.r2_model_novel.max():.3f} vs persistence={ce.r2_persistence_novel.max():.3f}")
    print(f"[3a sanity] beats rank-1 at ranks={ew.loc[ew.beats_rank1,'rank'].tolist()}; "
          f"elbow ~rank {int(ew.r2_vs_rank1.idxmax()) + 1}")


if __name__ == "__main__":
    main()
