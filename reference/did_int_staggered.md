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

## Examples

``` r
# 3 cohorts (t = 2, 3, 4) plus a never-treated group.
set.seed(7)
N <- 600; T <- 5
lon <- runif(N, 0, 10); lat <- runif(N, 0, 10)
z   <- 0.3 * lon + 0.2 * lat + rnorm(N)
p_t <- plogis(-0.5 + 0.5 * z)
is_t <- rbinom(N, 1, p_t) == 1
cohort <- rep(Inf, N)
cohort[is_t] <- sample(2:4, sum(is_t), replace = TRUE,
                       prob = c(0.4, 0.4, 0.2))
dij <- as.matrix(dist(cbind(lon, lat)))
A   <- (dij < 1.5) & (dij > 0)
deg <- pmax(rowSums(A), 1)
rows <- vector("list", N * T)
k <- 1L
for (i in seq_len(N)) for (t in seq_len(T)) {
  W_t <- as.integer(cohort[i] <= t)
  share_t <- sum(A[i, ] * (cohort <= t)) / deg[i]
  G_t <- as.integer(share_t > 0.3)
  Y <- 0.8 * z[i] + 0.1 * t * z[i] + 1.5 * W_t + 0.5 * G_t * W_t + rnorm(1)
  rows[[k]] <- data.frame(id = i, time = t, cohort = cohort[i],
                          z = z[i], Y = Y, G = G_t)
  k <- k + 1L
}
d <- do.call(rbind, rows)

# DR DATT at high exposure (g = 1) across cohort-time cells
res <- did_int_staggered(
  d, yname = "Y", time = "time", id = "id",
  cohort = "cohort", exposure = "G", g = 1, covariates = "z")
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: algorithm did not converge
head(res$per_cell)
#>   cohort time event_time estimate        se    ci_lo    ci_hi n_total n_at_g
#> 1      2    2          0 2.137885 0.1785172 1.787997 2.487772     600    336
#> 2      2    3          1 2.103313 0.1361382 1.836487 2.370139     452    452
#> 3      2    4          2 2.426680 0.1534515 2.125920 2.727439     377    377
#> 4      2    5          3 2.159710 0.1483221 1.869004 2.450416     377    377
#> 5      3    3          0 2.020326 0.1469968 1.732218 2.308435     411    411
#> 6      3    4          1 2.070073 0.1781373 1.720931 2.419216     336    336
#>   n_dropped
#> 1         0
#> 2         0
#> 3         0
#> 4         0
#> 5         0
#> 6         0
res$agg$simple   # joint-IF aggregate; truth is 2.0
#> $label
#> [1] "all cells"
#> 
#> $estimate
#> [1] 2.140156
#> 
#> $se
#> [1] 0.09898323
#> 
#> $ci
#> [1] 1.946152 2.334159
#> 
#> $n_cells
#> [1] 9
#> 
```
