#!/usr/bin/env python3
"""Step 3b, multi-seed — error bars on the out-of-panel completion margin.

The single-seed completion (operator_completion_stratified.py) reports the margin from
ONE random 20% holdout (seed 0). A reviewer flagged that the flagship +0.379 rests on a
single seed with no error bar. This repeats the exact same held-out prediction over K
independent random holdouts and reports the margin as MEAN ± SD per stratum, at the
rank-7 elbow, so "+0.379" becomes "+0.379 ± SD".

Same setup as operator_completion.condition_extrap (held-out (regulator, Stim48hr)
fibers for out-of-panel regulators, predicted from Rest+Stim8hr via the low-rank fit,
baseline = persistence). Strata are recomputed within each seed's held-out set:
STRONG/WEAK = median split by mean |ontarget_effect_size|; PANEL-MATCHED = strength >=
the in_original_panel median.

    python scripts/operator_completion_multiseed.py --seeds 20 --max-rank 12

Output: docs/tables/operator_completion_multiseed_3106.csv
        (stratum, rank, n_seeds, mean_margin, sd_margin, mean_r2_model, mean_r2_persistence)
"""
import argparse, time
from pathlib import Path
import numpy as np, pandas as pd
import _opkernels as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]
RANK_HL = 7   # elbow established in the single-seed run; headline reported here


def strength_map(regs):
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv")
    obs = obs[obs["ontarget_significant"].astype(str).str.strip().eq("True")]
    s = obs.assign(a=obs["ontarget_effect_size"].abs()).groupby("target_contrast_gene_name")["a"].mean()
    return np.array([s.get(r, np.nan) for r in regs])


def run_seed(tensor, mask, in_orig, strength_full, full, panel_med, seed, max_rank):
    R, G, C = tensor.shape
    late, early = COND_ORDER.index("Stim48hr"), COND_ORDER.index("Stim8hr")
    novel_full = full[~in_orig[full]]
    rng = np.random.default_rng(seed)
    n_test = max(20, len(full) // 5)
    test = rng.permutation(novel_full)[:n_test]
    if len(test) < n_test:
        test = np.concatenate([test, rng.permutation(np.setdiff1d(full, test))[: n_test - len(test)]])
    X = np.concatenate([tensor[full, :, 0], tensor[full, :, early], tensor[full, :, late]], axis=1)
    obs = np.ones_like(X, bool)
    test_pos = np.isin(full, test)
    novel_pos = test_pos & (~in_orig[full])
    obs[test_pos, 2 * G:3 * G] = False
    Xc, _ = op.train_test_standardize(X, obs)
    med = float(np.nanmedian(strength_full[novel_pos]))
    strata = {"all_novel": novel_pos,
              "strong": novel_pos & (strength_full >= med),
              "weak": novel_pos & (strength_full < med),
              "strong_panel_matched": novel_pos & (strength_full >= panel_med)}

    def block(Xhat, m):
        yt = Xc[np.ix_(m, np.arange(2 * G, 3 * G))]; yp = Xhat[np.ix_(m, np.arange(2 * G, 3 * G))]
        pers = Xc[np.ix_(m, np.arange(G, 2 * G))]; ss = np.sum((yt - yt.mean()) ** 2) + 1e-12
        return float(1 - np.sum((yt - yp) ** 2) / ss), float(1 - np.sum((yt - pers) ** 2) / ss)

    out = {}
    for r in range(1, max_rank + 1):
        Xhat = op.soft_impute(Xc, obs, r, n_iter=200)
        for name, m in strata.items():
            rm, rp = block(Xhat, m)
            out[(name, r)] = (rm, rp, rm - rp)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seeds", type=int, default=20)
    ap.add_argument("--max-rank", type=int, default=12); args = ap.parse_args()
    t0 = time.time()
    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    tensor, mask, in_orig = d["tensor"], d["mask"], d["in_original_panel"]
    regs = d["regulators"].astype(str)
    R = tensor.shape[0]
    full = np.array([i for i in range(R) if mask[i, :, :].all(axis=0).all()])
    strength = strength_map(regs); strength_full = strength[full]
    panel_med = float(np.nanmedian(strength[in_orig.astype(bool)]))
    print(f"full={len(full)} panel_med_strength={panel_med:.3f} seeds={args.seeds}", flush=True)

    acc = {}   # (stratum, rank) -> list of margins ; and model/pers
    for s in range(args.seeds):
        res = run_seed(tensor, mask, in_orig, strength_full, full, panel_med, s, args.max_rank)
        for k, (rm, rp, mg) in res.items():
            acc.setdefault(k, []).append((rm, rp, mg))
        print(f"[seed {s+1}/{args.seeds}] done ({time.time()-t0:.0f}s)", flush=True)

    rows = []
    for (name, r), vals in sorted(acc.items()):
        a = np.array(vals)  # cols: model, pers, margin
        rows.append(dict(stratum=name, rank=r, n_seeds=len(a),
                         mean_margin=round(a[:, 2].mean(), 4), sd_margin=round(a[:, 2].std(ddof=1), 4),
                         mean_r2_model=round(a[:, 0].mean(), 4), mean_r2_persistence=round(a[:, 1].mean(), 4)))
    df = pd.DataFrame(rows)
    TAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(TAB / "operator_completion_multiseed_3106.csv", index=False)
    print(f"\nHEADLINE (rank {RANK_HL}, mean ± SD over {args.seeds} seeds):", flush=True)
    for name in ["all_novel", "weak", "strong", "strong_panel_matched"]:
        row = df[(df.stratum == name) & (df["rank"] == RANK_HL)].iloc[0]
        print(f"  {name:22s} margin {row.mean_margin:+.3f} ± {row.sd_margin:.3f}  "
              f"(model {row.mean_r2_model:+.3f}, persistence {row.mean_r2_persistence:+.3f})", flush=True)
    print(f"\n-> {TAB/'operator_completion_multiseed_3106.csv'}  [{time.time()-t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
