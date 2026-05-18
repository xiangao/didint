# didint

Doubly robust difference-in-differences with spatial interference, following Xu (2023, 2026).

## What this package does

Standard DiD assumes one unit's outcome doesn't depend on another's treatment. When that's wrong — a treated municipality affecting its neighbours, a vaccine reducing transmission across the social network — the canonical DiD estimand loses its causal interpretation.

`didint` implements the doubly robust estimators of Ruonan Xu:

- **`did_int_2x2()`** — two-period, common-adoption-timing case from Xu (2023). Estimates the direct ATT at a chosen exposure level `g`.
- *Coming in subsequent versions:*
  - `did_int_dynamic()` — event study with common adoption timing (Xu 2026, Section I).
  - `did_int_staggered()` — staggered adoption with not-yet-treated comparison groups (Xu 2026, Section II).
  - Spillover-effect contrasts across exposure levels.

Standard errors come from the empirical influence function. With `coords` and a `cutoff`, they are Conley spatial-HAC.

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
  coords     = my_panel[, c("lon", "lat")],
  cutoff     = 100         # km, with spherical distance
)
print(res)
```

## References

- Xu, Ruonan (2023). "Difference-in-Differences with Interference." arXiv:2306.12003.
- Xu, Ruonan (2026). "Dynamic Difference-in-Differences with Interference." *AEA Papers and Proceedings* 116: 58–63.

## License

MIT
