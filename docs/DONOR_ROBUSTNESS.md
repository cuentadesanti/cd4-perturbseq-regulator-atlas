# Donor-level robustness — the real replication unit (critique #4)

Donors are the unit of replication in this screen (4 donors), so cross-donor reproducibility is the
robustness question that matters most. The original audit reweighted the EB score with donor metadata
present on only **19%** of contrasts — neutral-weighted (inert) for the other 81%. This replaces that
gesture with the dedicated per-donor DE object.

## Data

`GWCD4i.DE_stats.by_donors.h5mu` (streamed, obs only — no 15.7 GB download). Its `obs` carries
`donor_correlation_hits_{mean,min}` on **100%** of contrasts (3,993 KD-gated) — the mean / worst of
the **6 pairwise donor correlations** (4 donors → C(4,2)=6) of each perturbation's effect vector over
its hit genes. `donor_corr_hits_min` (worst pair) is the strict analog of "concordant in ≥3/4 donors".
Reproduce: `scripts/donor_concordance.py` → `docs/tables/donor_concordance.csv`.

## Findings

1. **Full coverage, real spread.** 3,993 KD-gated regulators, 100% with a donor-reproducibility
   value (vs 19% before). Distribution: median 0.74, 10th pct 0.27, range −0.72…0.99 — discriminating,
   not uniformly high.
2. **An independent axis.** Donor reproducibility is ~uncorrelated with effect size
   (ρ=0.045 vs `n_downstream`) and cell number (ρ=0.033) — it is genuinely new information, not a
   proxy for rank or power.
3. **The top regulators are donor-robust as a class.** Top-50 median 0.79 vs 0.73 for the rest
   (Mann-Whitney p=1.3e-4); **29/30** top regulators pass the hard flag (worst-pair ≥ 0.5). The core
   chromatin/TCR hits (TADA2B, CD3E, ZAP70, SGF29, TAF6L…) replicate across donors — the main claim
   survives the replication test.
4. **The payoff — large-effect hubs that do NOT replicate across donors** (demote candidates):

   | gene | condition | n_downstream | donor_corr_hits_mean | worst pair |
   |---|---|--:|--:|--:|
   | UBXN1 | Rest | 4439 | 0.28 | 0.21 |
   | DOP1B | Stim48hr | 3240 | 0.36 | 0.21 |
   | UBE2L3 | Stim8hr | 2873 | 0.43 | 0.33 |
   | ATF7IP2 | Stim8hr | 2730 | 0.36 | 0.27 |
   | SMG1 | Stim8hr | 2683 | 0.45 | 0.37 |
   | EIF1AX | Stim8hr | 2576 | 0.09 | −0.02 |

   These are exactly the false positives a count-based ranking misses: thousands of DE genes but no
   cross-donor replication (EIF1AX is essentially *anti*-correlated between donors). They should be
   demoted regardless of their rank.

## How to use it

`donor_robust` (worst-pair donor correlation ≥ 0.5) is a **hard survival flag**, not a continuous
reweight — a binary call on real per-donor data beats a weight that was inert for 81% of genes. Read
the ranking as: a large effect that is *also* donor-robust. The six hubs above fail that second test.
