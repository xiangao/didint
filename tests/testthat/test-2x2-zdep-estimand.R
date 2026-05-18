# Regression test for the estimand definition.
#
# The original 2x2 implementation accidentally averaged the DR signal over
# the {G = g} stratum instead of over the full population, and multiplied the
# regression term by 1{G = g}. With a treatment effect that depends on z and
# a {G = g} stratum whose z distribution differs from the full population,
# these two recipes give materially different numbers.
#
# Per Xu (2026, eq. 5; common-adoption case S_M = D_M), the correct estimand
# averages the DR signal over the full population. This test pins that down
# by simulating a DGP where the two recipes are ~0.25 apart and asserts the
# estimator targets the full-population value.

simulate_zdep <- function(N = 2000, seed = 1) {
  set.seed(seed)
  lon <- runif(N, 0, 10); lat <- runif(N, 0, 10)
  z   <- 0.3 * lon + 0.2 * lat + rnorm(N, sd = 1)
  W <- rbinom(N, 1, plogis(-0.5 + 0.6 * z))
  dij <- as.matrix(dist(cbind(lon, lat)))
  A <- (dij < 1.5) & (dij > 0)
  share <- (A %*% W) / pmax(rowSums(A), 1)
  G <- as.integer(share > median(share))
  Y_pre  <- 0.8 * z + rnorm(N)
  # Direct effect on treated at exposure g=1: 1.5 + 0.3*z + 0.5 = 2.0 + 0.3*z.
  # Averaged over full z distribution this is ~2.75 in this DGP; averaged
  # within G=1 (which has higher mean z) this is ~3.00. The estimator must
  # target the former.
  Y_post <- Y_pre + 0.2 * z + (1.5 + 0.3 * z) * W + 0.5 * G * W + rnorm(N)
  data.frame(W = W, G = G, z = z, Y_pre = Y_pre, Y_post = Y_post)
}

test_that("did_int_2x2 targets the full-population estimand, not the G=g stratum mean", {
  reps <- 100
  ests <- numeric(reps)
  for (r in seq_len(reps)) {
    df <- simulate_zdep(N = 2000, seed = r)
    ests[r] <- did_int_2x2(df, "Y_post", "Y_pre", "W", "G", g = 1,
                            covariates = "z", trim = 0.01)$estimate
  }
  truth_full   <- 2.75   # avg over full population
  truth_at_g   <- 3.00   # avg within G=g stratum (would be a sign of the bug)
  # MC SE for 100 reps with empirical SD around 0.2 is ~0.02; allow 0.10 slack.
  expect_lt(abs(mean(ests) - truth_full), 0.10)
  expect_gt(abs(mean(ests) - truth_at_g), 0.15)
})
