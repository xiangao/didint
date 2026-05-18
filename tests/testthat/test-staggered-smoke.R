# Smoke test for did_int_staggered:
# - 3 cohorts (treated at t=2, t=3, t=4) plus never-treated
# - 5 time periods (t=1..5); pre-period = 1
# - Time-varying exposure based on share of treated neighbours at each t
# - Constant direct effect = 1.5, constant spillover = 0.5
#   => DATT at g = 1 is 2.0 for every (c, t) cell.

simulate_staggered <- function(N = 800, T = 5, c_first = 2,
                               direct = 1.5, spill = 0.5,
                               seed = 7) {
  set.seed(seed)
  lon <- runif(N, 0, 10); lat <- runif(N, 0, 10)
  z   <- 0.3 * lon + 0.2 * lat + rnorm(N, sd = 1)
  # Assign cohort: probability of being treated in each cohort group is
  # decreasing in z so all groups have meaningful overlap.
  # We'll split treated units across c=2, 3, 4 with thirds, and never-treated.
  p_treated <- plogis(-0.5 + 0.5 * z)
  is_treated <- rbinom(N, 1, p_treated) == 1
  cohort <- rep(Inf, N)
  cohort_choice <- sample(c(c_first, c_first + 1, c_first + 2),
                          sum(is_treated), replace = TRUE,
                          prob = c(0.4, 0.4, 0.2))
  cohort[is_treated] <- cohort_choice

  # Neighbour structure (fixed across t)
  dij <- as.matrix(dist(cbind(lon, lat)))
  A   <- (dij < 1.5) & (dij > 0)
  deg <- pmax(rowSums(A), 1)

  rows <- vector("list", N * T)
  k <- 1L
  for (i in seq_len(N)) {
    for (t in seq_len(T)) {
      # Treatment indicator at time t
      W_t <- as.integer(cohort[i] <= t)
      # Exposure: share of neighbours treated at time t (binary on median later)
      share_t_i <- sum(A[i, ] * (cohort <= t)) / deg[i]
      # Outcome
      Y <- 0.8 * z[i] + 0.1 * t * z[i] +
           direct * W_t + spill * (share_t_i > 0.3) * W_t + rnorm(1)
      rows[[k]] <- data.frame(id = i, time = t, cohort = cohort[i],
                              z = z[i], lon = lon[i], lat = lat[i],
                              Y = Y, share = share_t_i)
      k <- k + 1L
    }
  }
  d <- do.call(rbind, rows)
  d$G <- as.integer(d$share > 0.3)
  d
}

test_that("did_int_staggered runs and recovers truth across all (c,t) cells", {
  d <- simulate_staggered(N = 2500, T = 5, c_first = 2, seed = 7)
  res <- did_int_staggered(
    data = d, yname = "Y", time = "time", id = "id",
    cohort = "cohort", exposure = "G", g = 1, covariates = "z"
    # trim = NULL by default; in staggered settings trim easily empties
    # the small (W, G=g) outcome subsets within each cell.
  )
  expect_s3_class(res, "didint_staggered")
  # With 3 cohorts (c=2,3,4) and times 2..5, we expect 9 cells:
  #   c=2 @ t=2,3,4,5 (4 cells); c=3 @ t=3,4,5 (3); c=4 @ t=4,5 (2)
  expect_equal(nrow(res$per_cell), 9)
  expect_true(all(res$per_cell$cohort %in% c(2, 3, 4)))
  expect_true(all(res$per_cell$time >= res$per_cell$cohort))

  # Aggregations exist
  expect_true(!is.null(res$agg$simple))
  expect_true(nrow(res$agg$event_time) > 0)
  expect_true(nrow(res$agg$cohort) > 0)

  # Truth (direct + spillover at g=1) is 2.0. Bias of simple average:
  bias <- res$agg$simple$estimate - 2.0
  expect_lt(abs(bias), 0.20)
})

test_that("did_int_staggered errors when no finite cohorts", {
  d <- simulate_staggered(N = 200, T = 3, seed = 1)
  d$cohort <- Inf
  expect_error(
    did_int_staggered(d, "Y", "time", "id", "cohort", "G",
                       g = 1, covariates = "z"),
    "no finite cohorts"
  )
})
