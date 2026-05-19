# didint

[![pkgdown](https://img.shields.io/badge/pkgdown-site-blue.svg)](https://xiangao.github.io/didint/)

Doubly robust difference-in-differences with spatial interference,
following Xu (2023, 2026).

## What this package does

Standard DiD assumes one unit’s outcome doesn’t depend on another’s
treatment. When that’s wrong — a treated municipality affecting its
neighbours, a vaccine reducing transmission across the social network —
the canonical DiD estimand loses its causal interpretation.

`didint` implements the doubly robust estimators of Ruonan Xu:

- **[`did_int_2x2()`](https://xiangao.github.io/didint/reference/did_int_2x2.md)**
  — two-period, common-adoption-timing case from Xu (2023). Estimates
  the direct ATT at a chosen exposure level `g`.
- **[`did_int_dynamic()`](https://xiangao.github.io/didint/reference/did_int_dynamic.md)**
  — event study with common adoption timing (Xu 2026, Section I).
- **[`did_int_staggered()`](https://xiangao.github.io/didint/reference/did_int_staggered.md)**
  — staggered adoption with not-yet-treated comparison groups (Xu 2026,
  Section II), with joint-IF aggregation across cohort-time cells.

Standard errors come from the empirical influence function. With
`coords` and a `cutoff`, they are Conley spatial-HAC. An optional `trim`
argument drops units with extreme propensities (Xu 2026 uses 0.01 in the
Brazil application).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("xiangao/didint")
```

## Documentation & vignettes

Full documentation: **<https://xiangao.github.io/didint/>**

| Page | Description |
|----|----|
| [Brazil Amazon Priority List](https://xiangao.github.io/didint/articles/brazil_amazon.html) | End-to-end replication of Xu (2026) Section III on the public Assunção-McMillan-Murphy-Souza-Rodrigues archive |
| [`did_int_2x2()`](https://xiangao.github.io/didint/reference/did_int_2x2.html) | 2x2 base case (Xu 2023) |
| [`did_int_dynamic()`](https://xiangao.github.io/didint/reference/did_int_dynamic.html) | Dynamic event study (Xu 2026 §I) |
| [`did_int_staggered()`](https://xiangao.github.io/didint/reference/did_int_staggered.html) | Staggered adoption with joint-IF aggregation (Xu 2026 §II) |
| [Reference index](https://xiangao.github.io/didint/reference/index.html) | All functions on one page |

A Julia port with identical estimators is at
[DidInterference.jl](https://github.com/xiangao/DidInterference.jl)
(docs: <https://xiangao.github.io/DidInterference.jl/>).

## References

- Xu, Ruonan (2023). “Difference-in-Differences with Interference.”
  [arXiv:2306.12003](https://arxiv.org/abs/2306.12003).
- Xu, Ruonan (2026). “Dynamic Difference-in-Differences with
  Interference.” *AEA Papers and Proceedings* 116: 58–63.

## License

MIT
