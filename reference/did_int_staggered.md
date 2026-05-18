# Staggered-adoption DR DATT with interference (Xu 2026, Section II)

Computes the doubly robust direct ATT at exposure level `g` for each
(cohort, period) cell with `t >= c`, using the not-yet-directly- treated
comparison group `{C > t}`. Returns per-cell estimates and simple
aggregations (simple average, by event-time, by cohort).

## Usage

``` r
did_int_staggered(
  data,
  yname,
  time,
  id,
  cohort,
  exposure,
  g,
  covariates,
  pre_period = NULL,
  cohorts = NULL,
  times = NULL,
  coords_cols = NULL,
  cutoff = NULL,
  dist_fn = c("spherical", "euclidean"),
  trim = NULL,
  alpha = 0.05
)
```

## Arguments

- data:

  Long-format panel: one row per `(id, time)`.

- yname:

  Outcome column.

- time:

  Time-period column.

- id:

  Unit identifier column.

- cohort:

  Cohort column; numeric, with `Inf` or `NA` for never-treated units.
  Treated units must have `cohort = c` for all their rows (i.e., cohort
  is time-invariant).

- exposure:

  Time-varying exposure column (one value per `(id, time)`).

- g:

  Target exposure level.

- covariates:

  Character vector of time-invariant attribute columns. Values at the
  post-period `t` are used (which equal the pre-period values when the
  column is truly time-invariant).

- pre_period:

  Baseline period. Defaults to `min(finite cohorts) - 1`.

- cohorts:

  Optional vector restricting which cohorts to estimate. Default: all
  finite cohorts.

- times:

  Optional vector restricting which post-periods to estimate. Default:
  all periods `>= min(cohorts)`.

- coords_cols:

  Optional length-2 character vector `c(lon, lat)` for spatial-HAC SEs.

- cutoff, dist_fn, trim, alpha:

  See
  [`did_int_2x2()`](https://xiangao.github.io/didint/reference/did_int_2x2.md).

## Value

A list of class `"didint_staggered"` with:

- per_cell:

  Data frame with one row per estimated `(c, t)` cell: `cohort`, `time`,
  `event_time = t - c`, `estimate`, `se`, `ci_lo`, `ci_hi`, `n_total`,
  `n_at_g`, `n_dropped`.

- agg:

  List of aggregated estimates with stacked-IF SEs: `simple` (average
  over all cells), `event_time` (data frame over `event_time`), `cohort`
  (data frame over `cohort`).

- influence:

  List of per-cell influence functions, indexed by the cell's row in
  `per_cell`. Each IF is aligned to the cell's own `S_M` subset, so they
  cannot be stacked unit-wise across cells; the aggregated SEs are
  computed by averaging within-cell contributions, weighted by cell
  size.

## Details

For each cell `(c, t)`:

1.  Restrict to `S_M = { i : C_i = c OR C_i > t }`.

2.  Compute `dY = Y_t - Y_{c_underbar - 1}` using `pre_period` (defaults
    to `min(finite cohorts) - 1`).

3.  Run the DR estimator (Xu 2026, eq. 5) with `W = 1{C_i = c}` and
    `Ig = 1{G_it = g}`.

Exposure is allowed to vary across periods (the column passed in
`exposure` should hold the time-varying `G_it`).

## See also

[`did_int_2x2()`](https://xiangao.github.io/didint/reference/did_int_2x2.md),
[`did_int_dynamic()`](https://xiangao.github.io/didint/reference/did_int_dynamic.md).
