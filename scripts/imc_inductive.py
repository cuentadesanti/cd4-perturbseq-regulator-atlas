#!/usr/bin/env python3
"""Bloque C · MVP — Inductive Matrix Completion: can we predict a NEVER-SEEN regulator?

Today's completion (Step 3b) is TRANSDUCTIVE: it only knows regulators present in the tensor and
holds out a condition *fiber* of an otherwise-observed regulator. This asks the harder, inductive
question that closes gap #3: predict the trans response of a regulator whose EVERY condition is
absent from training, using only its features (sequence/network/annotation). Y ~= X W Z^T; for a
held-out regulator i the prediction X_i W Z^T is computable WITHOUT ever seeing a row of i.

The design is SYMMETRIC on purpose — this is the one analysis that can come out WORSE than what we
have, and a null must be clean and reportable, not hidden:
  b0 mean-predictor   -> R^2 = 0 by construction (predict the training-mean response). IMC must beat this.
  b1 features-shuffled -> permute X rows (features de-aligned from Y) and re-run. If IMC does not beat
                          this, the features carry no signal and the result is NULL — reported as such.
  b2 (reference only)  -> the transductive R^2 (~0.06, BCV rank-7). Inductive will almost certainly be
                          LOWER (harder task); matching it is not the point.

MVP = ridge (a full-rank IMC and a legitimate inductive baseline). Low-rank X W Z^T is added ONLY if
ridge beats shuffled; if ridge is already null the MVP ends there — a clean, honest negative.

  python scripts/imc_inductive.py --folds 5 --seed 0

Inputs (already prepared, not re-downloaded): data/features/regulator_features_3106.csv,
data/cache/operator_tensor.npz.
Outputs: docs/tables/imc_feature_matrix_3106.csv, docs/tables/imc_leave_regulator_out_3106.csv.
"""
import argparse, time
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; FEAT = ROOT / "data" / "features"; TAB = ROOT / "docs" / "tables"
ALPHAS = np.logspace(2, 6, 13)           # ridge strength grid (selected by efficient LOO/GCV on train)
GO_MIN = 10; GO_CAP = 500; PF_MIN = 10; PF_CAP = 200


def _multihot(series, min_freq, cap, sep=";"):
    """binary columns for the most frequent tokens (freq >= min_freq, capped at `cap`)."""
    toks = series.fillna("").apply(lambda s: [t.strip() for t in str(s).split(sep) if t.strip()])
    freq = pd.Series([t for row in toks for t in row]).value_counts()
    keep = list(freq[freq >= min_freq].index[:cap])
    cols = {t: np.array([t in row for row in toks], float) for t in keep}
    return pd.DataFrame(cols), keep


NUMERIC = ["length", "string_degree", "string_wdegree", "n_go", "n_pfam"]


def build_X(feat, regs):
    """R^{3106 x d}: raw numerics (cols 0..4) + is_tf + dbd_family one-hot + sparse GO/Pfam multi-hot.

    Numerics are kept RAW here and standardized WITHIN each training fold (in run_cv) so there is no
    train/test leakage; the binary columns (is_tf/dbd/GO/Pfam) are left as 0/1 and never rescaled —
    rescaling rare 0/1 columns by their std amplifies them and blows ridge up numerically."""
    f = feat.set_index("Unnamed: 0").reindex(regs)
    num = f[NUMERIC].astype(float)
    num = num.fillna(num.median()).reset_index(drop=True)
    is_tf = f["is_tf"].astype(float).to_frame("is_tf").reset_index(drop=True)
    dbd = pd.get_dummies(f["dbd_family"].fillna("None"), prefix="dbd").astype(float).reset_index(drop=True)
    go, go_keep = _multihot(f["go_ids"].reset_index(drop=True), GO_MIN, GO_CAP)
    pf, pf_keep = _multihot(f["pfam_ids"].reset_index(drop=True), PF_MIN, PF_CAP)
    X = pd.concat([num, is_tf, dbd, go, pf], axis=1)      # numeric cols first (len(NUMERIC))
    X.index = regs
    print(f"  X built: {X.shape[1]} features = {len(NUMERIC)} numeric + 1 is_tf + {dbd.shape[1]} dbd + "
          f"{len(go_keep)} GO + {len(pf_keep)} Pfam", flush=True)
    return X


def per_reg_r2(Ytest, Yhat, ybar):
    """held-out R^2 per test regulator against the TRAINING mean (=> mean-predictor R^2 = 0)."""
    num = ((Ytest - Yhat) ** 2).sum(1)
    den = ((Ytest - ybar) ** 2).sum(1) + 1e-12
    return 1.0 - num / den


def run_cv(X, Y, folds, seed, tag, t0, shuffle_features=False):
    """5-fold leave-REGULATOR-out ridge; returns (per-fold [mean R^2, n_test, alpha], all per-reg R^2)."""
    Xv = X.values.astype(float)
    if shuffle_features:
        Xv = Xv[np.random.default_rng(seed + 777).permutation(len(Xv))]   # de-align features from Y
    nnum = len(NUMERIC)                  # standardize ONLY the numeric cols; leave binaries as 0/1
    kf = KFold(n_splits=folds, shuffle=True, random_state=seed)
    fold_mean, allr2 = [], []
    for i, (tr, te) in enumerate(kf.split(Xv)):
        Xtr, Xte = Xv[tr].copy(), Xv[te].copy()
        mu, sd = Xtr[:, :nnum].mean(0), Xtr[:, :nnum].std(0); sd[sd == 0] = 1.0
        Xtr[:, :nnum] = (Xtr[:, :nnum] - mu) / sd; Xte[:, :nnum] = (Xte[:, :nnum] - mu) / sd
        ybar = Y[tr].mean(0)
        with np.errstate(all="ignore"):        # macOS Accelerate emits spurious matmul FP warnings
            m = RidgeCV(alphas=ALPHAS).fit(Xtr, Y[tr] - ybar)  # alpha via efficient LOO on TRAIN only
            Yhat = m.predict(Xte) + ybar
        r2 = per_reg_r2(Y[te], Yhat, ybar)
        allr2.append(r2); fold_mean.append((float(np.mean(r2)), len(te), float(m.alpha_)))
        print(f"  [{tag}] fold {i+1}/{folds}: mean R^2={np.mean(r2):+.4f}  n_test={len(te)}  "
              f"alpha={m.alpha_:.1e}  ({time.time()-t0:.0f}s)", flush=True)
    return fold_mean, np.concatenate(allr2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=int, default=5); ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(); t0 = time.time()

    print(f"[load] tensor + features  ({time.time()-t0:.0f}s)", flush=True)
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    regs = np.array([str(x) for x in d["regulators"]])
    T = d["tensor"].astype(np.float64)
    Y = np.concatenate([T[:, :, k] for k in range(T.shape[2])], axis=1)     # 3106 x 6000 (cond-blocked)
    feat = pd.read_csv(FEAT / "regulator_features_3106.csv")
    assert set(feat["Unnamed: 0"]) == set(regs), "feature/tensor regulator set mismatch"

    print(f"[build X]  ({time.time()-t0:.0f}s)", flush=True)
    X = build_X(feat, regs)
    TAB.mkdir(parents=True, exist_ok=True)
    X.to_csv(TAB / "imc_feature_matrix_3106.csv")
    print(f"  Y target: {Y.shape} (concat of {T.shape[2]} conditions)  ({time.time()-t0:.0f}s)", flush=True)

    print(f"[run] real features, {args.folds}-fold leave-regulator-out  ({time.time()-t0:.0f}s)", flush=True)
    real_fold, real_all = run_cv(X, Y, args.folds, args.seed, "IMC", t0, shuffle_features=False)
    print(f"[run] b1 shuffled features (control)  ({time.time()-t0:.0f}s)", flush=True)
    shuf_fold, shuf_all = run_cv(X, Y, args.folds, args.seed, "shuf", t0, shuffle_features=True)

    rows = []
    for i, ((rm, ntest, alpha), (sm, _, _)) in enumerate(zip(real_fold, shuf_fold)):
        rows.append(dict(fold=i, rank="full_ridge", r2_imc=round(rm, 4), r2_ridge=round(rm, 4),
                         r2_shuffled=round(sm, 4), r2_mean=0.0, n_test=ntest, ridge_alpha=alpha))
    out = pd.DataFrame(rows)
    out.to_csv(TAB / "imc_leave_regulator_out_3106.csv", index=False)

    r_mean, r_med = float(real_all.mean()), float(np.median(real_all))
    s_mean = float(shuf_all.mean())
    pcts = np.percentile(real_all, [10, 25, 50, 75, 90])
    print("\n[Task C · leave-regulator-out] per-regulator held-out R^2 (n=%d):" % len(real_all))
    print(f"  IMC ridge   : mean={r_mean:.4f}  median={r_med:.4f}  "
          f"pct[10/25/50/75/90]={np.round(pcts,4).tolist()}")
    print(f"  b1 shuffled : mean={s_mean:.4f}   (features de-aligned from Y)")
    print(f"  b0 mean     : 0.0000  (by construction)")
    margin = r_mean - s_mean
    verdict = ("POSITIVE — features carry inductive signal (beats shuffled by %.4f)" % margin
               if (r_mean > 0 and margin > 0.005) else
               "NULL — features do NOT beat shuffled; unseen-regulator prediction unresolved on this data")
    print(f"\n[FLAGSHIP] inductive R^2={r_mean:.4f} vs shuffled={s_mean:.4f} -> {verdict}")
    print(f"  (reference: transductive BCV rank-7 R^2 ~ 0.06 — inductive is expected lower, that is correct)")
    print(f"  -> imc_leave_regulator_out_3106.csv, imc_feature_matrix_3106.csv  ({time.time()-t0:.0f}s)")
    print("[NEXT] add low-rank X W Z^T ONLY if ridge beats shuffled; else MVP ends on a clean negative.")


if __name__ == "__main__":
    main()
