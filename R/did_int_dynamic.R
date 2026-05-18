#' Dynamic DR DATT with interference (event study, common adoption)
#'
#' Implements Section I of Xu (2026). With common adoption timing `c`,
#' computes the doubly robust DATT at exposure level `g` for each
#' post-treatment period, taking the long difference `Y_t - Y_{c-1}`
#' as the outcome of interest.
#'
#' Because adoption timing is common, exposure does not drift across
#' post-periods: the 2x2 DR estimator is applied period-by-period with
#' the same treatment indicator and same exposure variable. The same
#' parallel-trends assumption must hold for *each* post-period
#' separately (Assumption 1 of Xu 2026 applied to every `t >= c`).
#'
#' @param data A data frame in wide format. Must contain a single
#'   pre-period outcome column and one column per post-period outcome.
#' @param yname_pre Character. Column name of the pre-period outcome
#'   (period `c - 1`).
#' @param ynames Character vector of column names for the post-period
#'   outcomes, ordered chronologically (period `c, c+1, ..., T`).
#' @param treat,exposure,g,covariates See [`did_int_2x2()`].
#' @param event_time Optional integer vector matching `ynames`, giving
#'   event time relative to treatment (e.g. `0:3` for c, c+1, c+2, c+3).
#'   Used to label the per-period results. Defaults to `seq_along(ynames) - 1`.
#' @param coords,cutoff,dist_fn,trim,alpha See [`did_int_2x2()`].
#' @param aggregate Logical. If `TRUE` (default), also returns a simple
#'   average across post-periods.
#'
#' @return A list of class `"didint_dynamic"` with:
#' \describe{
#'   \item{per_period}{Data frame: one row per post-period with
#'     `event_time`, `estimate`, `se`, `ci_lo`, `ci_hi`.}
#'   \item{agg}{If `aggregate = TRUE`, a list with `simple_avg`
#'     (mean of per-period estimates) and `se` (SE of the average,
#'     accounting for shared influence functions across periods).}
#'   \item{models}{Per-period model objects (returned by `did_int_2x2`).}
#' }
#'
#' @seealso [`did_int_2x2()`] for the 2x2 building block.
#'
#' @export
did_int_dynamic <- function(data, yname_pre, ynames, treat, exposure, g,
                            covariates, event_time = NULL,
                            coords = NULL, cutoff = NULL,
                            dist_fn = c("spherical", "euclidean"),
                            trim = NULL, alpha = 0.05,
                            aggregate = TRUE) {

  dist_fn <- match.arg(dist_fn)
  if (is.null(event_time)) event_time <- seq_along(ynames) - 1L
  stopifnot(length(event_time) == length(ynames))

  per <- vector("list", length(ynames))
  rows <- vector("list", length(ynames))

  for (k in seq_along(ynames)) {
    fit <- did_int_2x2(
      data        = data,
      yname       = ynames[k],
      yname_pre   = yname_pre,
      treat       = treat,
      exposure    = exposure,
      g           = g,
      covariates  = covariates,
      coords      = coords,
      cutoff      = cutoff,
      dist_fn     = dist_fn,
      trim        = trim,
      alpha       = alpha
    )
    per[[k]] <- fit
    rows[[k]] <- data.frame(
      event_time = event_time[k],
      estimate   = fit$estimate,
      se         = fit$se,
      ci_lo      = fit$ci[1],
      ci_hi      = fit$ci[2]
    )
  }

  per_period <- do.call(rbind, rows)
  rownames(per_period) <- NULL

  agg <- NULL
  if (aggregate) {
    # Simple average estimate
    est_avg <- mean(per_period$estimate)
    # SE for the average: stack the *per-period* influence functions and
    # average them, then take the iid (or Conley) variance of the stacked
    # average. With K periods sharing the same N units, the IF for the
    # average is (1/K) sum_k IF_k(i), and the variance is the iid sum.
    K <- length(per)
    n_units <- length(per[[1]]$influence)
    if_avg <- rowMeans(sapply(per, function(x) x$influence))
    se_avg <- sqrt(sum(if_avg^2) / n_units^2)
    z <- qnorm(1 - alpha / 2)
    agg <- list(
      simple_avg = est_avg,
      se         = se_avg,
      ci         = est_avg + c(-1, 1) * z * se_avg
    )
  }

  structure(
    list(per_period = per_period, agg = agg, models = per,
         exposure_g = g, alpha = alpha),
    class = "didint_dynamic"
  )
}

#' @export
print.didint_dynamic <- function(x, digits = 4, ...) {
  cat("Dynamic DR DATT (Xu 2026, common adoption)\n")
  cat(sprintf("  Exposure level g = %s\n", format(x$exposure_g)))
  cat("\nPer post-period:\n")
  pp <- x$per_period
  pp[, c("estimate", "se", "ci_lo", "ci_hi")] <-
    round(pp[, c("estimate", "se", "ci_lo", "ci_hi")], digits)
  print(pp, row.names = FALSE)
  if (!is.null(x$agg)) {
    cat(sprintf("\nSimple average across post-periods: %.*f (SE %.*f, %.0f%% CI [%.*f, %.*f])\n",
                digits, x$agg$simple_avg, digits, x$agg$se,
                100 * (1 - x$alpha),
                digits, x$agg$ci[1], digits, x$agg$ci[2]))
  }
  invisible(x)
}
