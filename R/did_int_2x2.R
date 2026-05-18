#' Doubly robust DiD with interference (2x2 base case)
#'
#' Implements the doubly robust estimator of Xu (2023) for the direct
#' average treatment effect on the treated (DATT) at a given exposure
#' level `g`, in the two-period, common-adoption-timing setting.
#'
#' Under correctly specified exposure mapping, conditional parallel
#' trends at the chosen exposure level, overlap, and no anticipation,
#' [`did_int_2x2()`] returns a consistent estimate of
#' \deqn{\tau_g = E[ y_{i,1}(1, g) - y_{i,1}(0, g) | z_i, W_i = 1, G_i = g ].}
#'
#' Three propensity models and two outcome-change models are fit:
#' \itemize{
#'   \item `p(z)  = P(W=1 | z)` — cohort propensity
#'   \item `pi_1g(z) = P(G=g | z, W=1)` — exposure prop. among treated
#'   \item `pi_0g(z) = P(G=g | z, W=0)` — exposure prop. among controls
#'   \item `m_1g(z) = E[dY | z, W=1, G=g]` — outcome change for treated
#'   \item `m_0g(z) = E[dY | z, W=0, G=g]` — outcome change for controls
#' }
#' The DR estimator is doubly robust: consistent if either all three
#' propensities OR both outcome models are correctly specified.
#'
#' Standard errors come from the empirical influence function. For
#' spatial inference, pass `coords` and `cutoff`; SEs are then computed
#' via the Conley spatial-HAC variance on the influence-function vector
#' (requires the `conleyreg` package).
#'
#' @param data A data frame.
#' @param yname Character. Column name for the post-period outcome.
#' @param yname_pre Character. Column name for the pre-period outcome.
#' @param treat Character. Column name of the binary treatment indicator
#'   `W_i` (post period; 0/1).
#' @param exposure Character. Column name of the exposure variable `G_i`
#'   (an integer or factor). Effects are computed at exposure level `g`.
#' @param g The exposure level at which to estimate the DATT.
#' @param covariates Character vector of column names for the attributes
#'   `z_i` used in all five working models.
#' @param coords Optional 2-column matrix or data frame of unit
#'   coordinates (e.g. longitude, latitude). When supplied together with
#'   `cutoff`, the standard error is the Conley spatial-HAC of the
#'   influence function.
#' @param cutoff Distance cutoff for the Conley kernel, in the same
#'   units as `coords` (km if coords are lon/lat with `dist_fn = "spherical"`).
#' @param dist_fn Either `"spherical"` (great-circle, expects lon/lat) or
#'   `"euclidean"`. Default `"spherical"`.
#' @param trim Optional propensity-score trimming threshold. If supplied,
#'   units with `p_hat <= trim` or `p_hat >= 1 - trim` are dropped before
#'   computing the DR estimate. Matches the trim at 0.01 used by Xu (2026)
#'   in the Brazil application. Default `NULL` (no trimming).
#' @param alpha Significance level for the CI; default 0.05.
#'
#' @return A list of class `"didint_2x2"` with:
#' \describe{
#'   \item{estimate}{The DR estimate of DATT at exposure `g`.}
#'   \item{se}{Standard error (iid by default; Conley if `coords` supplied).}
#'   \item{ci}{Two-element numeric vector with lower and upper CI bounds.}
#'   \item{n_treated}{Number of treated units used.}
#'   \item{n_control}{Number of control units used.}
#'   \item{n_total}{Total units satisfying the inclusion mask.}
#'   \item{call}{The matched call.}
#' }
#'
#' @references
#' Xu, Ruonan (2023). "Difference-in-Differences with Interference."
#' arXiv:2306.12003.
#'
#' Xu, Ruonan (2026). "Dynamic Difference-in-Differences with
#' Interference." AEA Papers and Proceedings 116: 58–63.
#'
#' @export
did_int_2x2 <- function(data, yname, yname_pre, treat, exposure, g,
                        covariates, coords = NULL, cutoff = NULL,
                        dist_fn = c("spherical", "euclidean"),
                        trim = NULL,
                        alpha = 0.05) {

  dist_fn <- match.arg(dist_fn)
  call <- match.call()

  # --- assemble vectors ----------------------------------------------------
  W  <- as.integer(data[[treat]])
  Gv <- data[[exposure]]
  dY <- data[[yname]] - data[[yname_pre]]
  Z  <- as.data.frame(data[, covariates, drop = FALSE])
  Ig <- as.integer(Gv == g)

  if (anyNA(W) || anyNA(Gv) || anyNA(dY) || anyNA(Z))
    stop("did_int_2x2: missing values in W, G, outcomes, or covariates")
  if (any(!(W %in% c(0, 1))))
    stop("did_int_2x2: treat must be 0/1")

  # Delegate the DR machinery to .dr_atte (shared with did_int_staggered).
  out <- .dr_atte(W = W, Ig = Ig, Z = Z, dY = dY,
                  coords = coords, cutoff = cutoff,
                  dist_fn = dist_fn, trim = trim, alpha = alpha)

  structure(
    c(out, list(exposure_g = g, call = call)),
    class = "didint_2x2"
  )
}

#' @export
print.didint_2x2 <- function(x, digits = 4, ...) {
  cat("Doubly robust DATT (Xu 2023, 2x2 case)\n")
  cat(sprintf("  Exposure level g = %s\n", format(x$exposure_g)))
  cat(sprintf("  N total = %d (treated %d, control %d), of which %d at exposure g\n",
              x$n_total, x$n_treated, x$n_control, x$n_at_g))
  if (!is.null(x$n_dropped) && x$n_dropped > 0)
    cat(sprintf("  Dropped by PS trim: %d\n", x$n_dropped))
  cat(sprintf("  DATT     = %.*f\n", digits, x$estimate))
  cat(sprintf("  SE       = %.*f\n", digits, x$se))
  cat(sprintf("  95%% CI  = [%.*f, %.*f]\n",
              digits, x$ci[1], digits, x$ci[2]))
  invisible(x)
}
