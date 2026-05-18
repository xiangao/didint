# Monte Carlo validation findings (Phase 4)

Reproduce: `Rscript inst/sims/mc_validation.R [reps]` from the package root.
Defaults to 200 reps per scenario; takes ~10 minutes.

## What works

**Point estimation (2×2 base case).** Bias → 0 as N grows; with `trim = 0.01`,
the estimator is approximately unbiased at N = 1500 (bias ≈ 0.001) and the
SE is mildly conservative. Coverage 0.92–0.94 across N = 500, 1500, 3000.

**Estimand correctness.** On a DGP where the direct effect depends on z
(true full-population value 2.75; true within-G=g value 3.00), the
estimator recovers 2.745 — the paper's full-population estimand, not the
within-stratum mean. Locked in by `test-2x2-zdep-estimand.R`.

**Dynamic event study.** Per-period DR estimates are unbiased; the
aggregated mean across post-periods has SE ≈ empirical SD and coverage
0.90 (slightly under nominal, likely from the assumption that per-period
IFs across periods are independent — they aren't, because the same units
appear in every period).

## Surprise finding: DR absorbs spatial structure

On a DGP with strong spatial correlation in `Y_post`
(Gaussian random field, range = 4, σ = 2 on a 10×10 grid), Conley spatial
HAC SEs are *identical to iid SEs* to three decimals (0.225 vs 0.225).
Coverage with iid is 0.955.

Direct diagnostic (`mean(IF_i × IF_j)` by distance bin, units 0–1 apart):
- raw `Y_post` residual: 3.14 (vs Var = 5.10) — strong spatial correlation
- doubly-robust IF: -0.16 (vs Var = 17.5) — essentially zero

The DR machinery — inverse weighting by three propensities + regression
adjustment — whitens the residual spatial structure when nuisances are
correctly specified. **Conley HAC may matter more under misspecification
or with weaker regression adjustment** (e.g., when the IPW terms don't
adequately balance the spatial covariates).

This is good news for inference robustness in this setting; it also means
that if a reviewer demands Conley SEs, you can show them they don't move
the estimate or interval.

## Staggered aggregated SE — fixed via joint-IF stacking

An earlier version of the staggered aggregation treated `(c, t)` cells
as independent. This gave SEs that undercovered:

| Estimator | empSD | meanSE | Coverage |
|---|---|---|---|
| staggered, 3 cohorts, N=2500 (independent cells)  | 0.051 | 0.037 | 0.840 |
| staggered, 3 cohorts, N=2500 (joint-IF stacking)  | 0.051 | 0.051 | 0.960 |

The fix: each cell's influence function is now stacked into the full
universe of unit IDs (zeros for units outside `S_M`), weighted by the
aggregation weight `w_k = n_k / Σ n_l`, summed across cells per unit,
and `Var(theta) ≈ Σ h_i²` where `h_i` is unit `i`'s contribution. This
correctly captures the positive cell-cell covariance from shared units
(the never-treated comparison reappears in every cell).

## What's not in this validation yet

- Misspecified nuisances (e.g., omitting an interaction in the outcome
  model). DR is supposed to be robust if at least one of (propensity,
  outcome) models is correctly specified; verifying this empirically is
  Phase 4 follow-up.
- Conley HAC under propensity-model misspecification — should make the
  spatial structure visible in the IF.
- Coverage of per-cohort and per-event-time aggregates in
  `did_int_staggered()`. Same independent-cells assumption applies; same
  bias direction expected.
- Time-varying covariates. Both estimators assume time-invariant `z`.
