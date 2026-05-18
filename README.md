# didint

Doubly robust difference-in-differences with spatial interference, following Xu (2023, 2026).

## What this package does

Standard DiD assumes one unit's outcome doesn't depend on another's treatment. When that's wrong — a treated municipality affecting its neighbours, a vaccine reducing transmission across the social network — the canonical DiD estimand loses its causal interpretation.

`didint` implements the doubly robust estimators of Ruonan Xu:

- **`did_int_2x2()`** — two-period, common-adoption-timing case from Xu (2023). Estimates the direct ATT at a chosen exposure level `g`.
- **`did_int_dynamic()`** — event study with common adoption timing (Xu 2026, Section I). Per-period direct ATTs plus a simple cross-period average.
- **`did_int_staggered()`** — staggered adoption with not-yet-treated comparison groups (Xu 2026, Section II). Per (cohort, time) cells plus three aggregations (simple, event-time, by-cohort) with joint-IF stacking across cells that share units.

Standard errors come from the empirical influence function. With `coords` and a `cutoff`, they are Conley spatial-HAC. An optional `trim` argument drops units with extreme propensities (Xu 2026 uses 0.01 in the Brazil application).

## Installation

```r
# install.packages("remotes")
remotes::install_github("xiangao/didint")
```

## Minimal example

```r
library(didint)
res <- did_int_2x2(
  data       = my_panel,
  yname      = "Y_post",
  yname_pre  = "Y_pre",
  treat      = "W",
  exposure   = "G",
  g          = 1,
  covariates = c("z1", "z2"),
  trim       = 0.01
)
print(res)
```

## Vignettes & examples

| Resource | Description |
|---|---|
| [`vignettes/brazil_amazon.Rmd`](https://github.com/xiangao/didint/blob/master/vignettes/brazil_amazon.Rmd) | End-to-end real-data replication of Xu (2026) Section III: Brazil Amazon *Lista de Municípios Prioritários* with `did_int_staggered()`, using the public Assunção-McMillan-Murphy-Souza-Rodrigues replication archive on Zenodo. Includes a [pre-rendered event-study figure](https://github.com/xiangao/didint/blob/master/vignettes/figures/brazil_event_study.png). |
| [`inst/sims/mc_validation.R`](https://github.com/xiangao/didint/blob/master/inst/sims/mc_validation.R) | Reproducible Monte Carlo validation script: bias, SE accuracy, and 95% CI coverage across N = 500–3000, with constant- and z-dependent-effect DGPs and a spatially-correlated-errors variant. Run with `Rscript inst/sims/mc_validation.R [reps]`. |
| [`inst/sims/findings.md`](https://github.com/xiangao/didint/blob/master/inst/sims/findings.md) | Write-up of MC findings, including the surprise that the DR estimator's influence function is essentially spatially uncorrelated even under strong spatial errors (Conley HAC ≡ iid here). |
| [`tests/testthat/`](https://github.com/xiangao/didint/tree/master/tests/testthat) | testthat suite (22 tests). Each `test-*.R` is a small worked example. |

A Julia port with identical estimators is at [DidInterference.jl](https://github.com/xiangao/DidInterference.jl).

## References

- Xu, Ruonan (2023). "Difference-in-Differences with Interference." [arXiv:2306.12003](https://arxiv.org/abs/2306.12003).
- Xu, Ruonan (2026). "Dynamic Difference-in-Differences with Interference." *AEA Papers and Proceedings* 116: 58–63.

## License

MIT
