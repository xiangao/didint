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

  # Early mask check: the DATT averages over units with G = g; if none exist,
  # nothing is identified. Catch before fitting outcome regressions on 0 obs.
  if (!any(Ig == 1))
    stop("did_int_2x2: no units with exposure level g; nothing to estimate")
  if (!any(W == 1 & Ig == 1))
    stop("did_int_2x2: no treated units at exposure level g")
  if (!any(W == 0 & Ig == 1))
    stop("did_int_2x2: no control units at exposure level g")

  # --- fit working models --------------------------------------------------
  rhs <- paste(c("1", covariates), collapse = " + ")

  # cohort propensity p(z) = P(W=1 | z)
  fml_p <- as.formula(paste("W ~", rhs))
  fit_p <- glm(fml_p, data = cbind(W = W, Z), family = binomial())
  p_hat <- predict(fit_p, type = "response")

  # exposure propensity among treated, among controls
  treated_idx <- which(W == 1)
  control_idx <- which(W == 0)

  fit_pi1 <- glm(as.formula(paste("Ig ~", rhs)),
                 data = cbind(Ig = Ig[treated_idx], Z[treated_idx, , drop = FALSE]),
                 family = binomial())
  fit_pi0 <- glm(as.formula(paste("Ig ~", rhs)),
                 data = cbind(Ig = Ig[control_idx], Z[control_idx, , drop = FALSE]),
                 family = binomial())
  pi1g_hat <- predict(fit_pi1, newdata = Z, type = "response")
  pi0g_hat <- predict(fit_pi0, newdata = Z, type = "response")

  # Optional propensity-score trim (Xu 2026 uses 0.01 in the Brazil example).
  # Drop units whose cohort *or* exposure propensities are in the tails;
  # extreme exposure propensities are the more common source of heavy-tailed
  # DR estimates in our simulation experience.
  n_dropped <- 0L
  if (!is.null(trim)) {
    keep <- p_hat    > trim & p_hat    < (1 - trim) &
            pi1g_hat > trim & pi1g_hat < (1 - trim) &
            pi0g_hat > trim & pi0g_hat < (1 - trim)
    n_dropped <- sum(!keep)
    if (n_dropped > 0) {
      W <- W[keep]; Gv <- Gv[keep]; dY <- dY[keep]
      Z <- Z[keep, , drop = FALSE]; Ig <- Ig[keep]
      p_hat <- p_hat[keep]; pi1g_hat <- pi1g_hat[keep]; pi0g_hat <- pi0g_hat[keep]
      if (!is.null(coords)) coords <- coords[keep, , drop = FALSE]
    }
  }

  # outcome change regressions on (W=1, G=g) and (W=0, G=g) subsets
  fit_m1 <- lm(as.formula(paste("dY ~", rhs)),
               data = cbind(dY = dY[W == 1 & Ig == 1],
                            Z[W == 1 & Ig == 1, , drop = FALSE]))
  fit_m0 <- lm(as.formula(paste("dY ~", rhs)),
               data = cbind(dY = dY[W == 0 & Ig == 1],
                            Z[W == 0 & Ig == 1, , drop = FALSE]))
  m1_hat <- predict(fit_m1, newdata = Z)
  m0_hat <- predict(fit_m0, newdata = Z)

  # --- DR influence-function contributions ---------------------------------
  # Mask: we form the DATT averaging over units with G_i = g (estimand restricts
  # the conditioning event to the {G=g} stratum). The early-exit check above
  # has already verified the mask is non-empty.
  mask <- Ig == 1

  # IF components (one per unit; zero outside mask), aligned to full N
  if_treated <- W * Ig / (p_hat * pi1g_hat) * (dY - m1_hat)
  if_control <- (1 - W) * Ig / ((1 - p_hat) * pi0g_hat) * (dY - m0_hat)
  if_reg     <- Ig * (m1_hat - m0_hat)

  psi <- if_treated - if_control + if_reg              # length n
  N   <- sum(mask)                                     # |{G = g}|
  est <- sum(psi) / N                                  # DR estimate

  # Empirical influence function: psi_i / Pr(G=g)  minus  est
  # In practice we use psi/N - est * (mask/N) so that mean(if_emp) = 0.
  if_emp <- (psi - est * Ig) / (N / length(W))         # scale to per-unit
  # The variance formula: Var(est) = (1/n^2) * sum(if^2) for iid;
  # Conley HAC otherwise.
  n <- length(W)

  # --- standard error: iid or spatial Conley HAC ---------------------------
  if (!is.null(coords) && !is.null(cutoff)) {
    if (!requireNamespace("conleyreg", quietly = TRUE))
      stop("did_int_2x2: install package 'conleyreg' for spatial-HAC SEs")
    # Regress IF on a constant; the (1,1) entry of the Conley vcov is the SE^2
    df_hac <- data.frame(if_emp = if_emp,
                         x_lon  = coords[, 1],
                         x_lat  = coords[, 2])
    fit_hac <- conleyreg::conleyreg(
      formula  = if_emp ~ 1,
      data     = df_hac,
      dist_cutoff = cutoff,
      lat = "x_lat", lon = "x_lon",
      kernel = "bartlett",
      dist_comp = if (dist_fn == "spherical") "spherical" else "planar"
    )
    se <- sqrt(diag(fit_hac$vcov)[1] / n)
  } else {
    se <- sqrt(sum(if_emp^2) / n^2)
  }

  z <- qnorm(1 - alpha / 2)
  ci <- est + c(-1, 1) * z * se

  structure(
    list(
      estimate    = est,
      se          = se,
      ci          = ci,
      n_treated   = sum(W == 1),
      n_control   = sum(W == 0),
      n_total     = n,
      n_dropped   = n_dropped,
      n_at_g      = N,
      exposure_g  = g,
      influence   = if_emp,
      models      = list(p = fit_p, pi1 = fit_pi1, pi0 = fit_pi0,
                         m1 = fit_m1, m0 = fit_m0),
      call        = call
    ),
    class = "didint_2x2"
  )
}

#' @export
print.didint_2x2 <- function(x, digits = 4, ...) {
  cat("Doubly robust DATT (Xu 2023, 2x2 case)\n")
  cat(sprintf("  Exposure level g = %s\n", format(x$exposure_g)))
  cat(sprintf("  N treated = %d, N control = %d, N at exposure g = %d\n",
              x$n_treated, x$n_control, x$n_at_g))
  cat(sprintf("  DATT     = %.*f\n", digits, x$estimate))
  cat(sprintf("  SE       = %.*f\n", digits, x$se))
  cat(sprintf("  95%% CI  = [%.*f, %.*f]\n",
              digits, x$ci[1], digits, x$ci[2]))
  invisible(x)
}
