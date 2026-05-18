# Dynamic DR DATT with interference (event study, common adoption)

Implements Section I of Xu (2026). With common adoption timing `c`,
computes the doubly robust DATT at exposure level `g` for each
post-treatment period, taking the long difference `Y_t - Y_{c-1}` as the
outcome of interest.

## Usage

``` r
did_int_dynamic(
  data,
  yname_pre,
  ynames,
  treat,
  exposure,
  g,
  covariates,
  event_time = NULL,
  coords = NULL,
  cutoff = NULL,
  dist_fn = c("spherical", "euclidean"),
  trim = NULL,
  alpha = 0.05,
  aggregate = TRUE
)
```

## Arguments

- data:

  A data frame in wide format. Must contain a single pre-period outcome
  column and one column per post-period outcome.

- yname_pre:

  Character. Column name of the pre-period outcome (period `c - 1`).

- ynames:

  Character vector of column names for the post-period outcomes, ordered
  chronologically (period `c, c+1, ..., T`).

- treat, exposure, g, covariates:

  See
  [`did_int_2x2()`](https://xiangao.github.io/didint/reference/did_int_2x2.md).

- event_time:

  Optional integer vector matching `ynames`, giving event time relative
  to treatment (e.g. `0:3` for c, c+1, c+2, c+3). Used to label the
  per-period results. Defaults to `seq_along(ynames) - 1`.

- coords, cutoff, dist_fn, trim, alpha:

  See
  [`did_int_2x2()`](https://xiangao.github.io/didint/reference/did_int_2x2.md).

- aggregate:

  Logical. If `TRUE` (default), also returns a simple average across
  post-periods.

## Value

A list of class `"didint_dynamic"` with:

- per_period:

  Data frame: one row per post-period with `event_time`, `estimate`,
  `se`, `ci_lo`, `ci_hi`.

- agg:

  If `aggregate = TRUE`, a list with `simple_avg` (mean of per-period
  estimates) and `se` (SE of the average, accounting for shared
  influence functions across periods).

- models:

  Per-period model objects (returned by `did_int_2x2`).

## Details

Because adoption timing is common, exposure does not drift across
post-periods: the 2x2 DR estimator is applied period-by-period with the
same treatment indicator and same exposure variable. The same
parallel-trends assumption must hold for *each* post-period separately
(Assumption 1 of Xu 2026 applied to every `t >= c`).

## See also

[`did_int_2x2()`](https://xiangao.github.io/didint/reference/did_int_2x2.md)
for the 2x2 building block.
