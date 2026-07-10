#!/usr/bin/env python3
"""Step 3b, stratified by KD strength — the persisted source for the strong/weak margin.

The aggregate out-of-panel result lives in `operator_completion_condition_3106.csv`
(621 novel regulators, margin +0.379 @ rank 7, beats persistence at every rank). The
strong/weak *stratified* margin was previously computed inline and never saved — this
script regenerates it as a committed artifact so any margin cited for a stratum has a
reproducible source on disk.

It replicates `operator_completion.condition_extrap` EXACTLY (same seed, same held-out
621, same low-rank fit) so the `all_novel` stratum reproduces the aggregate CSV, then
splits those same 621 regulators into strata and reports per-stratum R²/margin.

STRATUM DEFINITION (hardcoded and explicit — this is what was missing before):
    strength(regulator) = mean |ontarget_effect_size| over its significant rows.
    STRONG = out-of-panel test regulators with strength >= the MEDIAN strength of the
             held-out set;  WEAK = below the median.  (~half/half split.)
    Rationale: the median is a threshold with no free parameter and no dependence on the
    pilot's selection; we report whatever this reproducible split yields — the numbers
    are NOT forced to any previously-quoted value.

    python scripts/operator_completion_stratified.py --max-rank 12

Output: docs/tables/operator_completion_stratified_3106.csv
        (columns: stratum, n_regulators, rank, r2_model, r2_persistence, margin)
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import _opkernels as op

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "data" / "cache"; TAB = ROOT / "docs" / "tables"
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def regulator_strength(regs):
    """mean |ontarget_effect_size| per regulator (index-aligned to `regs`)."""
    obs = pd.read_csv(CACHE / "fingerprint_obs.csv")
    obs = obs[obs["ontarget_significant"].astype(str).str.strip().eq("True")]
    s = obs.assign(a=obs["ontarget_effect_size"].abs()).groupby("target_contrast_gene_name")["a"].mean()
    return np.array([s.get(r, np.nan) for r in regs])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--max-rank", type=int, default=12)
    ap.add_argument("--seed", type=int, default=0); args = ap.parse_args()

    d = np.load(CACHE / "operator_tensor.npz", allow_pickle=True)
    tensor, mask, in_orig = d["tensor"], d["mask"], d["in_original_panel"]
    regs = d["regulators"].astype(str)
    R, G, C = tensor.shape
    late, early = COND_ORDER.index("Stim48hr"), COND_ORDER.index("Stim8hr")

    # --- replicate condition_extrap setup EXACTLY (same seed => same held-out 621) ---
    full = np.array([i for i in range(R) if mask[i, :, :].all(axis=0).all()])
    novel_full = full[~in_orig[full]]
    rng = np.random.default_rng(args.seed)
    n_test = max(20, len(full) // 5)
    test = rng.permutation(novel_full)[:n_test]
    if len(test) < n_test:
        extra = rng.permutation(np.setdiff1d(full, test))[: n_test - len(test)]
        test = np.concatenate([test, extra])
    X = np.concatenate([tensor[full, :, 0], tensor[full, :, early], tensor[full, :, late]], axis=1)
    obs = np.ones_like(X, bool)
    test_pos = np.isin(full, test)
    novel_pos = test_pos & (~in_orig[full])
    obs[test_pos, 2 * G:3 * G] = False
    Xc, _ = op.train_test_standardize(X, obs)

    # --- strata among the held-out novel regulators ---
    # (1) MEDIAN split by strength (primary; no free parameter):
    strength = regulator_strength(regs)
    strength_full = strength[full]
    med = float(np.nanmedian(strength_full[novel_pos]))
    strong = novel_pos & (strength_full >= med)          # STRONG: strength >= held-out median
    weak = novel_pos & (strength_full < med)             # WEAK: strength < held-out median
    # (2) ORIGINAL-PANEL-MATCHED strong stratum (the apples-to-apples-to-pilot framing):
    #     novel regulators whose strength is at least the MEDIAN strength of the
    #     in_original_panel regulators — i.e. "as strong as the original panel".
    panel_med = float(np.nanmedian(strength[in_orig.astype(bool)]))
    panel_matched = novel_pos & (strength_full >= panel_med)
    strata = {"all_novel": novel_pos, "strong": strong, "weak": weak,
              "strong_panel_matched": panel_matched}
    print(f"held-out novel={int(novel_pos.sum())} | median strength(held-out)={med:.3f} "
          f"| strong={int(strong.sum())} weak={int(weak.sum())}")
    print(f"original-panel median strength={panel_med:.3f} | panel-matched strong={int(panel_matched.sum())}")

    def block(Xhat, posmask):
        yt = Xc[np.ix_(posmask, np.arange(2 * G, 3 * G))]
        yp = Xhat[np.ix_(posmask, np.arange(2 * G, 3 * G))]
        pers = Xc[np.ix_(posmask, np.arange(G, 2 * G))]
        ss = np.sum((yt - yt.mean()) ** 2) + 1e-12
        return float(1 - np.sum((yt - yp) ** 2) / ss), float(1 - np.sum((yt - pers) ** 2) / ss)

    rows = []
    for r in range(1, args.max_rank + 1):
        print(f"[3b-strat] rank {r}/{args.max_rank}", flush=True)
        Xhat = op.soft_impute(Xc, obs, r, n_iter=200)
        for name, m in strata.items():
            rm, rp = block(Xhat, m)
            rows.append(dict(stratum=name, n_regulators=int(m.sum()), rank=r,
                             r2_model=round(rm, 4), r2_persistence=round(rp, 4), margin=round(rm - rp, 4)))
    df = pd.DataFrame(rows)
    TAB.mkdir(parents=True, exist_ok=True)
    df.to_csv(TAB / "operator_completion_stratified_3106.csv", index=False)

    print("\npeak margin per stratum (max over ranks):")
    for name in strata:
        s = df[df.stratum == name]; b = s.loc[s.margin.idxmax()]
        print(f"  {name:9s} n={int(b.n_regulators):3d}  peak margin {b.margin:+.3f} @ rank {int(b['rank'])} "
              f"(model {b.r2_model:+.3f}, persistence {b.r2_persistence:+.3f})")
    print(f"\n-> {TAB / 'operator_completion_stratified_3106.csv'}")


if __name__ == "__main__":
    main()
