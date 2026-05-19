#' Staggered-adoption DR DATT with interference (Xu 2026, Section II)
#'
#' Computes the doubly robust direct ATT at exposure level `g` for each
#' (cohort, period) cell with `t >= c`, using the not-yet-directly-
#' treated comparison group `{C > t}`. Returns per-cell estimates and
#' simple aggregations (simple average, by event-time, by cohort).
#'
#' For each cell `(c, t)`:
#' \enumerate{
#'   \item Restrict to `S_M = { i : C_i = c  OR  C_i > t }`.
#'   \item Compute `dY = Y_t - Y_{c_underbar - 1}` using `pre_period`
#'         (defaults to `min(finite cohorts) - 1`).
#'   \item Run the DR estimator (Xu 2026, eq. 5) with
#'         `W = 1{C_i = c}` and `Ig = 1{G_it = g}`.
#' }
#'
#' Exposure is allowed to vary across periods (the column passed in
#' `exposure` should hold the time-varying `G_it`).
#'
#' @param data Long-format panel: one row per `(id, time)`.
#' @param yname Outcome column.
#' @param time Time-period column.
#' @param id Unit identifier column.
#' @param cohort Cohort column; numeric, with `Inf` or `NA` for
#'   never-treated units. Treated units must have `cohort = c` for all
#'   their rows (i.e., cohort is time-invariant).
#' @param exposure Time-varying exposure column (one value per
#'   `(id, time)`).
#' @param g Target exposure level.
#' @param covariates Character vector of time-invariant attribute
#'   columns. Values at the post-period `t` are used (which equal the
#'   pre-period values when the column is truly time-invariant).
#' @param pre_period Baseline period. Defaults to
#'   `min(finite cohorts) - 1`.
#' @param cohorts Optional vector restricting which cohorts to
#'   estimate. Default: all finite cohorts.
#' @param times Optional vector restricting which post-periods to
#'   estimate. Default: all periods `>= min(cohorts)`.
#' @param coords_cols Optional length-2 character vector
#'   `c(lon, lat)` for spatial-HAC SEs.
#' @param cutoff,dist_fn,trim,alpha See [`did_int_2x2()`].
#'
#' @return A list of class `"didint_staggered"` with:
#' \describe{
#'   \item{per_cell}{Data frame with one row per estimated `(c, t)`
#'     cell: `cohort`, `time`, `event_time = t - c`, `estimate`,
#'     `se`, `ci_lo`, `ci_hi`, `n_total`, `n_at_g`, `n_dropped`.}
#'   \item{agg}{List of aggregated estimates with stacked-IF SEs:
#'     `simple` (average over all cells), `event_time` (data frame
#'     over `event_time`), `cohort` (data frame over `cohort`).}
#'   \item{influence}{List of per-cell influence functions, indexed by
#'     the cell's row in `per_cell`. Each IF is aligned to the cell's
#'     own `S_M` subset, so they cannot be stacked unit-wise across
#'     cells; the aggregated SEs are computed by averaging within-cell
#'     contributions, weighted by cell size.}
#' }
#'
#' @seealso [`did_int_2x2()`], [`did_int_dynamic()`].
#'
#' @examples
#' # 3 cohorts (t = 2, 3, 4) plus a never-treated group.
#' set.seed(7)
#' N <- 600; T <- 5
#' lon <- runif(N, 0, 10); lat <- runif(N, 0, 10)
#' z   <- 0.3 * lon + 0.2 * lat + rnorm(N)
#' p_t <- plogis(-0.5 + 0.5 * z)
#' is_t <- rbinom(N, 1, p_t) == 1
#' cohort <- rep(Inf, N)
#' cohort[is_t] <- sample(2:4, sum(is_t), replace = TRUE,
#'                        prob = c(0.4, 0.4, 0.2))
#' dij <- as.matrix(dist(cbind(lon, lat)))
#' A   <- (dij < 1.5) & (dij > 0)
#' deg <- pmax(rowSums(A), 1)
#' rows <- vector("list", N * T)
#' k <- 1L
#' for (i in seq_len(N)) for (t in seq_len(T)) {
#'   W_t <- as.integer(cohort[i] <= t)
#'   share_t <- sum(A[i, ] * (cohort <= t)) / deg[i]
#'   G_t <- as.integer(share_t > 0.3)
#'   Y <- 0.8 * z[i] + 0.1 * t * z[i] + 1.5 * W_t + 0.5 * G_t * W_t + rnorm(1)
#'   rows[[k]] <- data.frame(id = i, time = t, cohort = cohort[i],
#'                           z = z[i], Y = Y, G = G_t)
#'   k <- k + 1L
#' }
#' d <- do.call(rbind, rows)
#'
#' # DR DATT at high exposure (g = 1) across cohort-time cells
#' res <- did_int_staggered(
#'   d, yname = "Y", time = "time", id = "id",
#'   cohort = "cohort", exposure = "G", g = 1, covariates = "z")
#' head(res$per_cell)
#' res$agg$simple   # joint-IF aggregate; truth is 2.0
#'
#' @export
did_int_staggered <- function(data, yname, time, id, cohort, exposure, g,
                              covariates,
                              pre_period = NULL,
                              cohorts = NULL, times = NULL,
                              coords_cols = NULL, cutoff = NULL,
                              dist_fn = c("spherical", "euclidean"),
                              trim = NULL, alpha = 0.05) {

  dist_fn <- match.arg(dist_fn)

  needed <- c(yname, time, id, cohort, exposure, covariates,
              if (!is.null(coords_cols)) coords_cols)
  miss <- setdiff(needed, names(data))
  if (length(miss))
    stop("did_int_staggered: missing columns: ",
         paste(miss, collapse = ", "))

  # Identify finite cohorts and time range.
  C_vals <- data[[cohort]]
  finite_cohorts <- sort(unique(C_vals[is.finite(C_vals)]))
  if (length(finite_cohorts) == 0)
    stop("did_int_staggered: no finite cohorts; nothing to estimate")
  c_underbar <- min(finite_cohorts)
  if (is.null(pre_period)) pre_period <- c_underbar - 1L
  if (is.null(cohorts))    cohorts    <- finite_cohorts
  T_max <- max(data[[time]])
  if (is.null(times))      times      <- seq.int(c_underbar, T_max)

  # Pre-period snapshot: one row per id with Y_pre.
  pre <- data[data[[time]] == pre_period, c(id, yname), drop = FALSE]
  names(pre)[2] <- ".Y_pre"

  # Container for results
  rows     <- list()
  ifs      <- list()    # per-cell influence-function vector
  cell_ids <- list()    # per-cell unit IDs aligned to ifs[[k]]
  k        <- 0L

  # All unit IDs (for joint-IF aggregation: stack to this universe)
  all_ids <- unique(data[[id]])

  for (c_val in cohorts) {
    for (t_val in times) {
      if (t_val < c_val) next

      # Snapshot at time t
      dt <- data[data[[time]] == t_val, ]
      m  <- merge(dt, pre, by = id, all.x = FALSE)
      if (nrow(m) == 0) next

      C_i  <- m[[cohort]]
      G_it <- m[[exposure]]
      Y_t  <- m[[yname]]
      Y_pre <- m[[".Y_pre"]]
      Z    <- as.data.frame(m[, covariates, drop = FALSE])
      ids_t <- m[[id]]
      coords <- if (!is.null(coords_cols))
        as.matrix(m[, coords_cols, drop = FALSE]) else NULL

      # S_M = { C = c OR C > t }
      in_sm <- (C_i == c_val) | (C_i > t_val)
      if (sum(in_sm) == 0) next

      W  <- as.integer(C_i[in_sm] == c_val)
      Ig <- as.integer(G_it[in_sm] == g)
      dY <- Y_t[in_sm] - Y_pre[in_sm]
      Z_sm <- Z[in_sm, , drop = FALSE]
      ids_sm <- ids_t[in_sm]
      coords_sm <- if (!is.null(coords)) coords[in_sm, , drop = FALSE] else NULL

      # Skip cells with empty {W=1, G=g} or {W=0, G=g} subsets — these
      # are not identifiable for this cell. Record a warning.
      if (!any(Ig == 1) || !any(W == 1 & Ig == 1) ||
          !any(W == 0 & Ig == 1)) {
        warning(sprintf(
          "did_int_staggered: cell (c=%s, t=%s) lacks units in {W=%s, G=%s}; skipped",
          format(c_val), format(t_val), "{0,1}", format(g)
        ))
        next
      }

      out <- tryCatch(
        .dr_atte(W = W, Ig = Ig, Z = Z_sm, dY = dY,
                 coords = coords_sm, cutoff = cutoff,
                 dist_fn = dist_fn, trim = trim, alpha = alpha),
        error = function(e) {
          warning(sprintf(
            "did_int_staggered: cell (c=%s, t=%s) failed: %s",
            format(c_val), format(t_val), conditionMessage(e)
          ))
          NULL
        }
      )
      if (is.null(out)) next

      k <- k + 1L
      rows[[k]] <- data.frame(
        cohort     = c_val,
        time       = t_val,
        event_time = t_val - c_val,
        estimate   = out$estimate,
        se         = out$se,
        ci_lo      = out$ci[1],
        ci_hi      = out$ci[2],
        n_total    = out$n_total,
        n_at_g     = out$n_at_g,
        n_dropped  = out$n_dropped
      )
      ifs[[k]]      <- out$influence
      cell_ids[[k]] <- ids_sm[out$keep_idx]   # IDs aligned to IF post-trim
    }
  }

  if (k == 0L)
    stop("did_int_staggered: no cell could be estimated")

  per_cell <- do.call(rbind, rows)
  rownames(per_cell) <- NULL

  # --- aggregation ---------------------------------------------------------
  # Joint-IF aggregation across cells that share units.
  #
  # Each cell estimate is tau_k = (1/n_k) sum_{i in S_M^k} psi_{k,i}.
  # A weighted average theta = sum_k w_k * tau_k has the per-unit
  # contribution h_i = sum_k w_k * (psi_{k,i} / n_k) * 1{i in S_M^k}.
  # Treating units as independent (correct for iid sampling):
  #   Var(theta) = sum_i Var(h_i) ≈ sum_i h_i^2
  # No need to multiply by 1/n_universe because h_i is already the
  # contribution to theta (an average), not a sum.
  agg_one <- function(idx, label) {
    if (length(idx) == 0) return(NULL)
    ests <- per_cell$estimate[idx]
    ns   <- per_cell$n_total[idx]
    w    <- ns / sum(ns)
    est_avg <- sum(w * ests)

    # Build per-unit contribution h_i across cells in this aggregate.
    h <- numeric(length(all_ids))
    names(h) <- as.character(all_ids)
    for (j in seq_along(idx)) {
      kk    <- idx[j]
      psi_k <- ifs[[kk]] + per_cell$estimate[kk]   # un-centered DR signal
      n_k   <- per_cell$n_total[kk]
      ids_k <- as.character(cell_ids[[kk]])
      # Each cell-unit pair's centered contribution to its cell estimate
      # is (psi_{k,i} - tau_k) / n_k; aggregated contribution to theta is
      # w[j] times this.
      contrib <- w[j] * (psi_k - per_cell$estimate[kk]) / n_k
      h[ids_k] <- h[ids_k] + contrib
    }
    se_avg <- sqrt(sum(h^2))
    z      <- qnorm(1 - alpha / 2)
    list(label = label, estimate = est_avg, se = se_avg,
         ci = est_avg + c(-1, 1) * z * se_avg,
         n_cells = length(idx))
  }

  agg_simple <- agg_one(seq_len(nrow(per_cell)), "all cells")

  agg_event <- do.call(rbind, lapply(
    sort(unique(per_cell$event_time)),
    function(et) {
      a <- agg_one(which(per_cell$event_time == et),
                   sprintf("event_time=%s", et))
      data.frame(event_time = et, estimate = a$estimate, se = a$se,
                 ci_lo = a$ci[1], ci_hi = a$ci[2], n_cells = a$n_cells)
    }
  ))
  rownames(agg_event) <- NULL

  agg_cohort <- do.call(rbind, lapply(
    sort(unique(per_cell$cohort)),
    function(c_val) {
      a <- agg_one(which(per_cell$cohort == c_val),
                   sprintf("cohort=%s", c_val))
      data.frame(cohort = c_val, estimate = a$estimate, se = a$se,
                 ci_lo = a$ci[1], ci_hi = a$ci[2], n_cells = a$n_cells)
    }
  ))
  rownames(agg_cohort) <- NULL

  structure(
    list(per_cell = per_cell,
         agg = list(simple = agg_simple,
                    event_time = agg_event,
                    cohort = agg_cohort),
         influence = ifs,
         exposure_g = g,
         pre_period = pre_period,
         alpha = alpha),
    class = "didint_staggered"
  )
}

#' @export
print.didint_staggered <- function(x, digits = 4, ...) {
  cat("Staggered DR DATT with interference (Xu 2026)\n")
  cat(sprintf("  Exposure level g = %s; pre-period = %s; %d (c,t) cells estimated\n",
              format(x$exposure_g), format(x$pre_period), nrow(x$per_cell)))
  cat("\nSimple average across cells:\n")
  with(x$agg$simple, cat(sprintf(
    "  estimate = %.*f   se = %.*f   %.0f%% CI = [%.*f, %.*f]\n",
    digits, estimate, digits, se, 100 * (1 - x$alpha),
    digits, ci[1], digits, ci[2]
  )))
  cat("\nBy event time:\n")
  pp <- x$agg$event_time
  pp[, c("estimate", "se", "ci_lo", "ci_hi")] <-
    round(pp[, c("estimate", "se", "ci_lo", "ci_hi")], digits)
  print(pp, row.names = FALSE)
  invisible(x)
}
