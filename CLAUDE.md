# didint — project notes for Claude

R package implementing the doubly robust DiD with spatial interference
estimators of Xu (2023, arXiv:2306.12003) and Xu (2026, *AEA P&P* 116:
58–63).

Companion Julia port: [`DidInterference.jl`](https://github.com/xiangao/DidInterference.jl).
Both packages are tested against the same DGPs and should produce
matching estimates up to MC noise — if you change the DR core in one,
mirror the change in the other.

## What's where

- `R/dr_atte.R` — internal `.dr_atte()` helper: the doubly-robust core
  (3 propensity scores + 2 outcome regressions + plug-in formula +
  influence-function SE). Returns `keep_idx` so callers can align IF
  values back to their input rows after PS trimming.
- `R/did_int_2x2.R` — Xu (2023) 2×2 case. Thin wrapper around `.dr_atte()`.
- `R/did_int_dynamic.R` — Xu (2026) §I event study under common adoption.
- `R/did_int_staggered.R` — Xu (2026) §II staggered adoption.
  **Crucial subtlety**: the cross-cell aggregation uses joint-IF
  stacking (per-cell IFs aligned to unit IDs and summed across cells),
  not independent-cells SEs. Cells share the never-treated comparison,
  so independent-cells SEs underestimate by ~30% (verified in
  `inst/sims/`).
- `inst/sims/mc_validation.R` — Monte Carlo bias/coverage harness.
- `inst/sims/findings.md` — write-up of MC results plus the surprise
  that Conley HAC ≡ iid SE on this DR estimator (DR machinery whitens
  the IF's spatial structure when nuisances are correct).
- `vignettes/brazil_amazon.Rmd` — Brazil Amazon Priority List
  replication (Xu 2026 Section III) using the public Assunção et al.
  ReStud Zenodo archive. `eval = FALSE` because the data isn't bundled.

## Estimand pin

`did_int_2x2()` averages the DR signal over the **full sample** (S_M),
not the `{G = g}` stratum. This is the paper's definition. There's a
regression test in `tests/testthat/test-2x2-zdep-estimand.R` that
pins this — don't drop it. (Original bug had the average over
`{G = g}`; only caught when treatment effect varied with z.)

## Workflows

- **Tests**: small-N smoke tests with `compute_se = FALSE` for speed
  (~30s total). Real validation is in `inst/sims/mc_validation.R`.
- **Docs**: pkgdown auto-builds and deploys to `gh-pages` on push to
  master. Live at <https://xiangao.github.io/didint/>. Roxygen man
  pages are committed (regenerate with `roxygen2::roxygenise()`).
- **Pages config**: GH Pages source = `gh-pages` branch root (one-time
  manual setup; already done).

## When changing the DR core

1. Edit `R/dr_atte.R`.
2. Re-run tests (`testthat::test_dir("tests/testthat")`).
3. Mirror the same change in `DidInterference.jl/src/dr_atte.jl`.
4. Update `inst/sims/findings.md` if numerical results materially shift.
