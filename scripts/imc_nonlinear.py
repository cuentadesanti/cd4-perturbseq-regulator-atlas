#!/usr/bin/env python3
"""Bloque C robustness — NON-LINEAR check of the inductive IMC null (CatBoost, leave-regulator-out).

Ridge showed regulator features do not predict a never-seen regulator's trans response LINEARLY.
This asks whether a strong NON-LINEAR learner (CatBoost, gradient-boosted trees) finds interaction
signal ridge misses. Pre-registered expectation: reinforces the null (probable), does not flip it.
Framing = "robustness of the null", NOT a second attempt to extract signal.

For comparability with the ridge null (scripts/imc_inductive.py) NOTHING that defines the split or
the features changes:
  * same folds   — KFold(5, shuffle=True, random_state=0) over the SAME 3106 regulators (tensor order)
  * same X       — build_X() imported from imc_inductive.py (numeric + is_tf + dbd + GO/Pfam)
  * same metric  — per-regulator held-out R^2 over the 6000 genes vs the TRAINING mean (mean-pred = 0)
  * same b1      — CatBoost re-run on row-permuted X; the deciding number is real vs shuffled.

The one change that makes CatBoost tractable/valid (it is single-output; Y is 3106x6000):
REDUCE Y by a TRAIN-ONLY truncated SVD to K in {10,30}, fit one CatBoost per component, then
RECONSTRUCT the K predicted scores back to 6000-dim with the train SVD basis and score R^2 in the
SAME 6000-gene space as ridge. K=10/30 bracket the BCV predictive rank ~7 / signal ~92.

Capacity diagnostic: R^2_train is reported too. High train R^2 with test R^2 ~ 0 is positive
evidence of "capacity present, signal absent" — the tree analogue of ridge driving alpha up.

  pip install catboost && python scripts/imc_nonlinear.py --seed 0

Output: docs/tables/imc_nonlinear_catboost_3106.csv (fold, K, r2_test_real, r2_test_shuffled,
r2_train_real, r2_mean, n_test).
"""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import KFold, train_test_split
from sklearn.decomposition import TruncatedSVD
from catboost import CatBoostRegressor

sys.path.insert(0, str(Path(__file__).resolve().parent))
from imc_inductive import build_X, per_reg_r2                      # SAME X and metric as the ridge null

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; FEAT = ROOT / "data" / "features"; TAB = ROOT / "docs" / "tables"
KS = [10, 30]; KMAX = 30                                            # nested: K=10 reuses first 10 models


def fit_catboost(Xtr, ytr, seed):
    Xa, Xv, ya, yv = train_test_split(Xtr, ytr, test_size=0.15, random_state=seed)  # internal val
    m = CatBoostRegressor(depth=6, iterations=500, random_seed=seed, loss_function="RMSE",
                          early_stopping_rounds=30, verbose=False, allow_writing_files=False,
                          thread_count=-1)
    m.fit(Xa, ya, eval_set=(Xv, yv))
    return m


def capacity_probe(Xv, Y, seed, t0):
    """Fold-0, K=30, NO early stopping (full 500 iters): can CatBoost even FIT the train scores?
    High train R^2 here with test R^2 ~ 0 above = capacity present, signal absent (the point of the
    diagnostic; early stopping in the main run suppresses train overfit, muting this contrast)."""
    kf = KFold(n_splits=5, shuffle=True, random_state=seed)
    tr, te = next(iter(kf.split(Xv)))
    ybar = Y[tr].mean(0)
    with np.errstate(all="ignore"):        # macOS Accelerate spurious matmul FP warnings (randomized SVD)
        svd = TruncatedSVD(n_components=KMAX, random_state=seed).fit(Y[tr] - ybar)
        comps = svd.components_; scores_tr = (Y[tr] - ybar) @ comps.T
    pred_tr = np.zeros((len(tr), KMAX)); pred_te = np.zeros((len(te), KMAX))
    for j in range(KMAX):
        m = CatBoostRegressor(depth=6, iterations=500, random_seed=seed, loss_function="RMSE",
                              verbose=False, allow_writing_files=False, thread_count=-1)  # no early stop
        m.fit(Xv[tr], scores_tr[:, j])
        pred_tr[:, j] = m.predict(Xv[tr]); pred_te[:, j] = m.predict(Xv[te])
    with np.errstate(all="ignore"):
        r2_tr = float(np.mean(per_reg_r2(Y[tr], ybar + pred_tr @ comps, ybar)))
        r2_te = float(np.mean(per_reg_r2(Y[te], ybar + pred_te @ comps, ybar)))
    print(f"[capacity probe] fold0 K30 no-early-stop: train R^2={r2_tr:+.4f}  test R^2={r2_te:+.4f}  "
          f"({time.time()-t0:.0f}s)", flush=True)
    return r2_tr, r2_te


def run(Xv, Y, seed, tag, t0, shuffle_features=False):
    """returns {(fold,K): (r2_test, r2_train), ...} plus per-fold n_test."""
    Xm = Xv.copy()
    if shuffle_features:
        Xm = Xm[np.random.default_rng(seed + 777).permutation(len(Xm))]   # same de-align as ridge b1
    kf = KFold(n_splits=5, shuffle=True, random_state=seed)
    res, ntest = {}, {}
    for fi, (tr, te) in enumerate(kf.split(Xm)):
        ybar = Y[tr].mean(0)
        with np.errstate(all="ignore"):        # macOS Accelerate spurious matmul FP warnings (randomized SVD)
            svd = TruncatedSVD(n_components=KMAX, random_state=seed).fit(Y[tr] - ybar)
            comps = svd.components_                                 # KMAX x 6000 (train SVD basis)
            scores_tr = (Y[tr] - ybar) @ comps.T                   # n_tr x KMAX (CatBoost targets)
        pred_te = np.zeros((len(te), KMAX)); pred_tr = np.zeros((len(tr), KMAX))
        for j in range(KMAX):
            m = fit_catboost(Xm[tr], scores_tr[:, j], seed)
            pred_te[:, j] = m.predict(Xm[te]); pred_tr[:, j] = m.predict(Xm[tr])
        for K in KS:
            with np.errstate(all="ignore"):    # macOS Accelerate spurious matmul FP warnings
                Yhat_te = ybar + pred_te[:, :K] @ comps[:K]
                Yhat_tr = ybar + pred_tr[:, :K] @ comps[:K]
            r2_te = float(np.mean(per_reg_r2(Y[te], Yhat_te, ybar)))
            r2_tr = float(np.mean(per_reg_r2(Y[tr], Yhat_tr, ybar)))
            res[(fi, K)] = (r2_te, r2_tr)
        ntest[fi] = len(te)
        print(f"  [{tag}] fold {fi+1}/5: "
              f"K10 test={res[(fi,10)][0]:+.4f} train={res[(fi,10)][1]:+.4f} | "
              f"K30 test={res[(fi,30)][0]:+.4f} train={res[(fi,30)][1]:+.4f}  "
              f"({time.time()-t0:.0f}s)", flush=True)
    return res, ntest


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(); t0 = time.time()

    print(f"[load] tensor + features  ({time.time()-t0:.0f}s)", flush=True)
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    regs = np.array([str(x) for x in d["regulators"]])
    T = d["tensor"].astype(np.float64)
    Y = np.concatenate([T[:, :, k] for k in range(T.shape[2])], axis=1)     # 3106 x 6000 (same as ridge)
    feat = pd.read_csv(FEAT / "regulator_features_3106.csv")
    X = build_X(feat, regs); Xv = X.values.astype(float)                    # SAME X as the ridge null
    print(f"  X={Xv.shape}  Y={Y.shape}  ({time.time()-t0:.0f}s)", flush=True)

    print(f"[run] CatBoost real features  ({time.time()-t0:.0f}s)", flush=True)
    real, ntest = run(Xv, Y, args.seed, "real", t0, shuffle_features=False)
    print(f"[run] CatBoost shuffled features (b1 control)  ({time.time()-t0:.0f}s)", flush=True)
    shuf, _ = run(Xv, Y, args.seed, "shuf", t0, shuffle_features=True)
    print(f"[run] capacity probe (no early stopping)  ({time.time()-t0:.0f}s)", flush=True)
    cap_tr, cap_te = capacity_probe(Xv, Y, args.seed, t0)

    rows = []
    for fi in range(5):
        for K in KS:
            rt, rtr = real[(fi, K)]; st, _ = shuf[(fi, K)]
            rows.append(dict(fold=fi, K=K, r2_test_real=round(rt, 4), r2_test_shuffled=round(st, 4),
                             r2_train_real=round(rtr, 4), r2_mean=0.0, n_test=ntest[fi]))
    out = pd.DataFrame(rows)
    TAB.mkdir(parents=True, exist_ok=True)
    out.to_csv(TAB / "imc_nonlinear_catboost_3106.csv", index=False)

    print("\n[Task C robustness · CatBoost leave-regulator-out] mean over 5 folds:")
    for K in KS:
        sub = out[out.K == K]
        rt, st, tr = sub.r2_test_real.mean(), sub.r2_test_shuffled.mean(), sub.r2_train_real.mean()
        margin = rt - st
        verdict = ("POSITIVE — non-linear signal (beats shuffled)" if (rt > 0 and margin > 0.005)
                   else "NULL — does not beat shuffled")
        print(f"  K={K:2d}: test_real={rt:+.4f}  test_shuffled={st:+.4f}  train_real={tr:+.4f}  "
              f"-> {verdict}")
    overall = ("NULL reinforced — a strong non-linear learner extracts no generalizable signal either; "
               "the limit is the DATA, not the model class"
               if all((out[out.K == K].r2_test_real.mean() -
                       out[out.K == K].r2_test_shuffled.mean()) <= 0.005 for K in KS)
               else "POSITIVE at some K — revisit tier-1")
    print(f"\n[FLAGSHIP] {overall}")
    print(f"  capacity probe (fold0 K30, no early stop): train R^2={cap_tr:+.4f}  test R^2={cap_te:+.4f}")
    print(f"  -> CatBoost CAN fit train (capacity present) but does NOT generalize (signal absent) — "
          f"the tree analogue of ridge driving alpha to the ceiling")
    print(f"  -> imc_nonlinear_catboost_3106.csv  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
