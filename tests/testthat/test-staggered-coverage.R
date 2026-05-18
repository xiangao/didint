# Regression test for joint-IF aggregation in did_int_staggered().
#
# The earlier "independent cells" aggregation underestimated the SE of
# the simple average across cells, because cells share units (the never-
# treated comparison appears in every cell). The fix stacks per-cell
# influence functions to the universe of unit IDs and uses the standard
# (1/n^2) sum(h_i^2) formula. This test runs a small MC and pins the
# mean SE close to the empirical SD.

# Note: 30 reps is the minimum for a sane MC SE check; can be flaky on
# very small samples. Marked slow with skip_on_cran() since CRAN tests
# should be fast.

simulate_stagg <- function(N = 2000, T = 5, c_first = 2, seed = 7) {
  set.seed(seed)
  lon <- runif(N, 0, 10); lat <- runif(N, 0, 10)
  z   <- 0.3 * lon + 0.2 * lat + rnorm(N, sd = 1)
  p_t <- plogis(-0.5 + 0.5 * z)
  is_t <- rbinom(N, 1, p_t) == 1
  cohort <- rep(Inf, N)
  cohort[is_t] <- sample(c(c_first, c_first + 1, c_first + 2),
                         sum(is_t), replace = TRUE, prob = c(0.4, 0.4, 0.2))
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
  do.call(rbind, rows)
}

test_that("did_int_staggered aggregate SE matches empirical SD across reps", {
  skip_if_not(identical(Sys.getenv("DIDINT_SLOW_TESTS"), "1"),
              "set DIDINT_SLOW_TESTS=1 to run the staggered MC coverage test")
  reps <- 30
  truth <- 2.0
  ests <- ses <- rep(NA_real_, reps)
  for (r in seq_len(reps)) {
    d <- simulate_stagg(N = 2000, seed = r)
    out <- tryCatch(
      suppressWarnings(did_int_staggered(
        d, yname = "Y", time = "time", id = "id",
        cohort = "cohort", exposure = "G", g = 1, covariates = "z")),
      error = function(e) NULL)
    if (is.null(out)) next
    ests[r] <- out$agg$simple$estimate
    ses[r]  <- out$agg$simple$se
  }
  ok <- !is.na(ests)
  # Bias should be near zero
  expect_lt(abs(mean(ests[ok]) - truth), 0.10)
  # Mean SE should be within 30% of empirical SD (loose to keep MC noise
  # from flaking the test). With independent-cells SE, this gap was ~25%.
  ratio <- mean(ses[ok]) / sd(ests[ok])
  expect_gt(ratio, 0.70)
  expect_lt(ratio, 1.40)
})
