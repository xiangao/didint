# Internal: core DR estimator for direct ATT at exposure level g.
#
# Takes already-prepared per-unit vectors and returns the estimate, SE,
# influence function, and basic counts. Both did_int_2x2() and
# did_int_staggered() prepare inputs in their own way (one is just W
# vs. (1 - W), the other is 1{C = c} vs. 1{C > t} within S_M), then
# call this helper.
#
# Inputs are assumed already validated; the helper performs minimal
# checks so callers don't double-validate.
#
# @param W   integer 0/1 indicator of treated group ("treated cohort").
# @param Ig  integer 0/1 indicator of exposure G == g.
# @param Z   data.frame of covariates (time-invariant or evaluated at t).
# @param dY  numeric outcome change Y_t - Y_pre.
# @param coords,cutoff,dist_fn,trim,alpha  as in did_int_2x2.
#
# @return a list with estimate, se, ci, influence (length n), n_total,
#   n_treated, n_control, n_at_g, n_dropped, working-model fits.
.dr_atte <- function(W, Ig, Z, dY,
                     coords = NULL, cutoff = NULL,
                     dist_fn = "spherical",
                     trim = NULL, alpha = 0.05) {

  if (length(W) != length(Ig) || length(W) != length(dY) ||
      length(W) != nrow(Z))
    stop(".dr_atte: input lengths do not match")
  if (!any(Ig == 1))
    stop(".dr_atte: no units with G = g")
  if (!any(W == 1 & Ig == 1))
    stop(".dr_atte: no treated units at G = g (cannot fit m_1g)")
  if (!any(W == 0 & Ig == 1))
    stop(".dr_atte: no control units at G = g (cannot fit m_0g)")

  covs <- names(Z)
  rhs  <- paste(c("1", covs), collapse = " + ")

  # Cohort propensity p(z) = P(W = 1 | z) within the sample
  fit_p <- glm(as.formula(paste("W ~", rhs)),
               data = cbind(W = W, Z), family = binomial())
  p_hat <- predict(fit_p, type = "response")

  # Exposure propensities, fit separately on W==1 and W==0 subsamples,
  # then predicted for all units.
  treated_idx <- which(W == 1)
  control_idx <- which(W == 0)
  fit_pi1 <- glm(as.formula(paste("Ig ~", rhs)),
                 data = cbind(Ig = Ig[treated_idx],
                              Z[treated_idx, , drop = FALSE]),
                 family = binomial())
  fit_pi0 <- glm(as.formula(paste("Ig ~", rhs)),
                 data = cbind(Ig = Ig[control_idx],
                              Z[control_idx, , drop = FALSE]),
                 family = binomial())
  pi1g_hat <- predict(fit_pi1, newdata = Z, type = "response")
  pi0g_hat <- predict(fit_pi0, newdata = Z, type = "response")

  # Optional PS trim on all three propensities (Xu 2026 uses 0.01).
  # Track which input rows survive so callers can align the returned
  # influence function back to their original sample.
  keep_idx <- seq_along(W)
  n_dropped <- 0L
  if (!is.null(trim)) {
    keep <- p_hat    > trim & p_hat    < (1 - trim) &
            pi1g_hat > trim & pi1g_hat < (1 - trim) &
            pi0g_hat > trim & pi0g_hat < (1 - trim)
    n_dropped <- sum(!keep)
    if (n_dropped > 0) {
      keep_idx <- keep_idx[keep]
      W <- W[keep]; Ig <- Ig[keep]; dY <- dY[keep]
      Z <- Z[keep, , drop = FALSE]
      p_hat <- p_hat[keep]
      pi1g_hat <- pi1g_hat[keep]; pi0g_hat <- pi0g_hat[keep]
      if (!is.null(coords)) coords <- coords[keep, , drop = FALSE]
    }
  }

  # Post-trim checks: trim can empty the (W=1, G=g) or (W=0, G=g) outcome-
  # regression subsets even when they were non-empty pre-trim.
  if (!is.null(trim)) {
    if (!any(W == 1 & Ig == 1))
      stop(".dr_atte: trim emptied the (W=1, G=g) subset; lower trim or skip cell")
    if (!any(W == 0 & Ig == 1))
      stop(".dr_atte: trim emptied the (W=0, G=g) subset; lower trim or skip cell")
  }

  # Outcome-change regressions on the (W=1, G=g) and (W=0, G=g) subsets,
  # predicted at every unit's z.
  fit_m1 <- lm(as.formula(paste("dY ~", rhs)),
               data = cbind(dY = dY[W == 1 & Ig == 1],
                            Z[W == 1 & Ig == 1, , drop = FALSE]))
  fit_m0 <- lm(as.formula(paste("dY ~", rhs)),
               data = cbind(dY = dY[W == 0 & Ig == 1],
                            Z[W == 0 & Ig == 1, , drop = FALSE]))
  m1_hat <- predict(fit_m1, newdata = Z)
  m0_hat <- predict(fit_m0, newdata = Z)

  # DR signal per unit (Xu 2026 eq. 5). Regression term contributes to all
  # units in S_M, not just G = g.
  if_treated <- W * Ig / (p_hat * pi1g_hat) * (dY - m1_hat)
  if_control <- (1 - W) * Ig / ((1 - p_hat) * pi0g_hat) * (dY - m0_hat)
  if_reg     <- m1_hat - m0_hat
  psi <- if_treated - if_control + if_reg

  n   <- length(W)
  est <- sum(psi) / n
  if_emp <- psi - est

  if (!is.null(coords) && !is.null(cutoff)) {
    if (!requireNamespace("conleyreg", quietly = TRUE))
      stop(".dr_atte: install package 'conleyreg' for spatial-HAC SEs")
    # conleyreg requires at least one RHS variable; we regress the
    # influence-function values on a constant column (with intercept = FALSE)
    # to recover the spatial-HAC variance of the mean.
    df_hac <- data.frame(if_emp = if_emp,
                         const  = 1,
                         x_lon  = coords[, 1],
                         x_lat  = coords[, 2])
    vcov_mat <- conleyreg::conleyreg(
      formula  = if_emp ~ const, data = df_hac,
      dist_cutoff = cutoff,
      lat = "x_lat", lon = "x_lon", kernel = "bartlett",
      intercept = FALSE,
      dist_comp = if (dist_fn == "spherical") "spherical" else "planar",
      vcov = TRUE, verbose = FALSE
    )
    # conleyreg with vcov = TRUE returns the variance-covariance matrix
    # directly. The coefficient on 'const' equals the mean of if_emp;
    # its sampling variance estimated by conleyreg already corresponds
    # to the spatial-HAC variance of that mean.
    se <- sqrt(as.numeric(vcov_mat[1, 1]))
  } else {
    se <- sqrt(sum(if_emp^2) / n^2)
  }

  z <- qnorm(1 - alpha / 2)
  list(
    estimate  = est,
    se        = se,
    ci        = est + c(-1, 1) * z * se,
    influence = if_emp,
    keep_idx  = keep_idx,           # input rows surviving the trim
    n_total   = n,
    n_treated = sum(W == 1),
    n_control = sum(W == 0),
    n_at_g    = sum(Ig == 1),
    n_dropped = n_dropped,
    models    = list(p = fit_p, pi1 = fit_pi1, pi0 = fit_pi0,
                     m1 = fit_m1, m0 = fit_m0)
  )
}
