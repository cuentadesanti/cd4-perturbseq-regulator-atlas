#!/usr/bin/env python3
"""Block A · Task 1a — cross-validated PREDICTIVE rank of the operator (Owen–Perry spirit).

The completion flagship (Step 3b) reported the predictive rank from a single seed-0 holdout.
This replaces it with a DISTRIBUTION: bi-cross-validation of the condition-extrapolation task —
fold over regulators (rows) x hold out one condition block (columns), predict it from the other
conditions via the low-rank fit, and record the rank that MAXIMISES held-out R^2. Over the
seed x (held-out condition) x (holdout fraction) grid this gives the predictive rank with a CI.

  Primary config: hold out Stim48hr for 20% of regulators (the 3b setup), 20 seeds -> rank CI.
  Robustness: also hold out Rest and Stim8hr, and sweep holdout 10/20/40% — if the optimal rank
  moves with which fiber/fraction you hold out, it is config-specific and we say so.

The point (contrast with guardrail A): the predictive rank (~7) is FAR below the 92 signal modes
over the empirical noise edge (2.95). Only ~7 of ~92 signal directions GENERALISE across the
condition split; the rest is regulator-idiosyncratic — which honestly explains the R^2~=0.07
ceiling as mostly irreducible, not a model failure.

    python scripts/operator_bcv_rank.py --max-rank 15 --seeds 20

Output: docs/tables/operator_bcv_rank_3106.csv
"""
import argparse, time
from pathlib import Path
import numpy as np, pandas as pd
import _opkernels as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"
COND = ["Rest", "Stim8hr", "Stim48hr"]


def peak_rank(T, io, cond_idx, hold, ranks, seed):
    """Hold out condition `cond_idx` for a random `hold` fraction of the OUT-OF-PANEL regulators
    (matching the 3b flagship's generalisation-to-unseen), predict from the other conditions via
    soft_impute; return (rank of peak held-out R2, peak R2, n_test)."""
    R, G, C = T.shape
    X = np.concatenate([T[:, :, k] for k in range(C)], axis=1)   # condition-BLOCKED [Rest|Stim8|Stim48]
    cols = np.arange(cond_idx * G, (cond_idx + 1) * G)           # so this is condition cond_idx's block
    rng = np.random.default_rng(seed)
    oop = np.where(~io)[0]                                        # out-of-panel pool (3023)
    test = np.zeros(R, bool); test[rng.choice(oop, int(hold * len(oop)), replace=False)] = True
    obs = np.ones_like(X, bool); obs[np.ix_(test, cols)] = False
    Xs, _ = op.train_test_standardize(X, obs)
    yt = Xs[np.ix_(test, cols)]; ss = np.sum((yt - yt.mean()) ** 2) + 1e-12
    r2s = []
    for r in ranks:
        yp = np.nan_to_num(op.soft_impute(Xs, obs, r, n_iter=100)[np.ix_(test, cols)],
                           nan=0.0, posinf=0.0, neginf=0.0)
        r2 = float(1 - np.sum((yt - yp) ** 2) / ss)
        r2s.append(r2 if np.isfinite(r2) else -1e9)
    r2s = np.array(r2s)
    return int(ranks[int(np.argmax(r2s))]), round(float(r2s.max()), 4), int(test.sum())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--max-rank", type=int, default=15)
    ap.add_argument("--seeds", type=int, default=20); args = ap.parse_args()
    t0 = time.time()
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    T = d["tensor"].astype(np.float64); io = d["in_original_panel"].astype(bool)
    ranks = list(range(1, args.max_rank + 1))
    rows = []

    # ---- primary: hold out Stim48hr @ 20%, many seeds -> predictive-rank distribution ----
    late = COND.index("Stim48hr")
    for s in range(args.seeds):
        r_opt, r2, ntest = peak_rank(T, io, late, 0.20, ranks, s)
        rows.append(dict(config="primary", held_out=COND[late], holdout_frac=0.20, seed=s,
                         optimal_rank=r_opt, peak_r2=r2, n_test=ntest))
    prim = pd.DataFrame([r for r in rows if r["config"] == "primary"]).optimal_rank
    print(f"[BCV primary] Stim48hr@20%, {args.seeds} seeds: predictive rank median={int(prim.median())} "
          f"mean={prim.mean():.1f} 95%CI=[{np.percentile(prim,2.5):.0f},{np.percentile(prim,97.5):.0f}]  "
          f"({time.time()-t0:.0f}s)", flush=True)

    # ---- robustness: each held-out condition @20%, and Stim48hr @10/40% (3 seeds each) ----
    grid = [(0, 0.2), (1, 0.2), (2, 0.1), (2, 0.4)]   # Rest@20, Stim8@20, Stim48@10, Stim48@40
    for c, hold in grid:
        opt = []
        for s in range(3):
            r_opt, r2, ntest = peak_rank(T, io, c, hold, ranks, 100 + s)
            opt.append(r_opt)
            rows.append(dict(config="robustness", held_out=COND[c], holdout_frac=hold, seed=100 + s,
                             optimal_rank=r_opt, peak_r2=r2, n_test=ntest))
        print(f"[robust] hold {COND[c]:8s} @ {int(hold*100)}%: rank median={int(np.median(opt))} "
              f"[{min(opt)}-{max(opt)}]  ({time.time()-t0:.0f}s)", flush=True)

    TAB.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TAB / "operator_bcv_rank_3106.csv", index=False)
    print(f"\n[FLAGSHIP] cross-validated predictive rank = {int(prim.median())} "
          f"[95%CI {np.percentile(prim,2.5):.0f}-{np.percentile(prim,97.5):.0f}] vs 92 signal modes "
          f"(guardrail A): only ~{int(prim.median())} of ~92 signal directions generalise across "
          f"conditions. -> operator_bcv_rank_3106.csv ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
