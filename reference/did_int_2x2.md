# Doubly robust DiD with interference (2x2 base case)

Implements the doubly robust estimator of Xu (2023) for the direct
average treatment effect on the treated (DATT) at a given exposure level
`g`, in the two-period, common-adoption-timing setting.

## Usage

``` r
did_int_2x2(
  data,
  yname,
  yname_pre,
  treat,
  exposure,
  g,
  covariates,
  coords = NULL,
  cutoff = NULL,
  dist_fn = c("spherical", "euclidean"),
  trim = NULL,
  alpha = 0.05
)
```

## Arguments

- data:

  A data frame.

- yname:

  Character. Column name for the post-period outcome.

- yname_pre:

  Character. Column name for the pre-period outcome.

- treat:

  Character. Column name of the binary treatment indicator `W_i` (post
  period; 0/1).

- exposure:

  Character. Column name of the exposure variable `G_i` (an integer or
  factor). Effects are computed at exposure level `g`.

- g:

  The exposure level at which to estimate the DATT.

- covariates:

  Character vector of column names for the attributes `z_i` used in all
  five working models.

- coords:

  Optional 2-column matrix or data frame of unit coordinates (e.g.
  longitude, latitude). When supplied together with `cutoff`, the
  standard error is the Conley spatial-HAC of the influence function.

- cutoff:

  Distance cutoff for the Conley kernel, in the same units as `coords`
  (km if coords are lon/lat with `dist_fn = "spherical"`).

- dist_fn:

  Either `"spherical"` (great-circle, expects lon/lat) or `"euclidean"`.
  Default `"spherical"`.

- trim:

  Optional propensity-score trimming threshold. If supplied, units with
  `p_hat <= trim` or `p_hat >= 1 - trim` are dropped before computing
  the DR estimate. Matches the trim at 0.01 used by Xu (2026) in the
  Brazil application. Default `NULL` (no trimming).

- alpha:

  Significance level for the CI; default 0.05.

## Value

A list of class `"didint_2x2"` with:

- estimate:

  The DR estimate of DATT at exposure `g`.

- se:

  Standard error (iid by default; Conley if `coords` supplied).

- ci:

  Two-element numeric vector with lower and upper CI bounds.

- n_treated:

  Number of treated units used.

- n_control:

  Number of control units used.

- n_total:

  Total units satisfying the inclusion mask.

- call:

  The matched call.

## Details

Under correctly specified exposure mapping, conditional parallel trends
at the chosen exposure level, overlap, and no anticipation,
`did_int_2x2()` returns a consistent estimate of \$\$\tau_g = E\[
y\_{i,1}(1, g) - y\_{i,1}(0, g) \| z_i, W_i = 1, G_i = g \].\$\$

Three propensity models and two outcome-change models are fit:

- `p(z) = P(W=1 | z)` — cohort propensity

- `pi_1g(z) = P(G=g | z, W=1)` — exposure prop. among treated

- `pi_0g(z) = P(G=g | z, W=0)` — exposure prop. among controls

- `m_1g(z) = E[dY | z, W=1, G=g]` — outcome change for treated

- `m_0g(z) = E[dY | z, W=0, G=g]` — outcome change for controls

The DR estimator is doubly robust: consistent if either all three
propensities OR both outcome models are correctly specified.

Standard errors come from the empirical influence function. For spatial
inference, pass `coords` and `cutoff`; SEs are then computed via the
Conley spatial-HAC variance on the influence-function vector (requires
the `conleyreg` package).

## References

Xu, Ruonan (2023). "Difference-in-Differences with Interference."
arXiv:2306.12003.

Xu, Ruonan (2026). "Dynamic Difference-in-Differences with
Interference." AEA Papers and Proceedings 116: 58–63.
