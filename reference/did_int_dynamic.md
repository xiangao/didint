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

## Examples

``` r
# Common adoption, 4 post-periods. Same DGP as did_int_2x2() example
# but with multiple post outcomes.
set.seed(2)
N <- 800; T_post <- 4
lon <- runif(N, 0, 10); lat <- runif(N, 0, 10)
z   <- 0.3 * lon + 0.2 * lat + rnorm(N)
W <- rbinom(N, 1, plogis(-0.5 + 0.6 * z))
dij <- as.matrix(dist(cbind(lon, lat)))
A <- (dij < 1.5) & (dij > 0)
share <- (A %*% W) / pmax(rowSums(A), 1)
G <- as.integer(share > median(share))
Y_pre <- 0.8 * z + rnorm(N)
df <- data.frame(W = W, G = G, z = z, Y_pre = Y_pre,
                 lon = lon, lat = lat)
for (k in 1:T_post)
  df[[paste0("Y_post_", k)]] <- Y_pre + 0.2 * k * z + 1.5 * W +
                                 0.5 * G * W + rnorm(N)

res <- did_int_dynamic(df, yname_pre = "Y_pre",
                       ynames = paste0("Y_post_", 1:T_post),
                       treat = "W", exposure = "G", g = 1,
                       covariates = "z", trim = 0.01)
res$per_period      # per-period direct ATTs
#>   event_time estimate        se    ci_lo    ci_hi
#> 1          0 2.224815 0.1695760 1.892452 2.557177
#> 2          1 2.251935 0.1389546 1.979589 2.524281
#> 3          2 2.197539 0.1500105 1.903524 2.491554
#> 4          3 1.949359 0.1665906 1.622847 2.275871
res$agg$simple_avg  # cross-period average; truth is 2.0
#> [1] 2.155912
```
