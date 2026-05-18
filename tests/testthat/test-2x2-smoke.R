# Smoke test: simulate a small 2x2 DGP with spatial interference and check
# the DR estimator recovers the truth within a generous tolerance. Full
# Monte-Carlo coverage tests come in Phase 4.

simulate_2x2 <- function(N = 800, true_dte = 1.5, true_spill = 0.5,
                         seed = 42) {
  set.seed(seed)
  # Coordinates on a unit square
  lon <- runif(N, 0, 10)
  lat <- runif(N, 0, 10)
  # One attribute z (smooth in space, not too informative)
  z   <- 0.3 * lon + 0.2 * lat + rnorm(N, sd = 1)
  # Treatment: logistic in z
  W <- rbinom(N, 1, plogis(-0.5 + 0.6 * z))
  # Exposure G: share of treated neighbours within distance 1.5; binary by median
  dij <- as.matrix(dist(cbind(lon, lat)))
  A   <- (dij < 1.5) & (dij > 0)
  share_trt_nbrs <- (A %*% W) / pmax(rowSums(A), 1)
  G <- as.integer(share_trt_nbrs > median(share_trt_nbrs))
  # Outcomes:
  # Y_pre depends only on z (no treatment yet, no spillover)
  Y_pre  <- 0.8 * z + rnorm(N)
  # Y_post = Y_pre + parallel trend (= 0.2 * z) + direct + spillover
  Y_post <- Y_pre + 0.2 * z + true_dte * W + true_spill * G * W + rnorm(N)
  data.frame(W = W, G = G, z = z, Y_pre = Y_pre, Y_post = Y_post,
             lon = lon, lat = lat)
}

test_that("did_int_2x2 runs and returns plausible DATT at g = 1", {
  df  <- simulate_2x2(N = 1500, true_dte = 1.5, true_spill = 0.5, seed = 1)
  res <- did_int_2x2(
    data     = df,
    yname    = "Y_post",
    yname_pre = "Y_pre",
    treat    = "W",
    exposure = "G",
    g        = 1,
    covariates = "z"
  )
  expect_s3_class(res, "didint_2x2")
  expect_true(is.finite(res$estimate))
  expect_true(is.finite(res$se) && res$se > 0)
  # True DATT at g=1 is direct + spillover = 2.0; allow ~2x SE tolerance.
  truth <- 1.5 + 0.5
  expect_lt(abs(res$estimate - truth), 4 * res$se)
})

test_that("did_int_2x2 errors on missing exposure level", {
  df <- simulate_2x2(N = 300, seed = 2)
  df$G[df$G == 1] <- 0  # eliminate g = 1
  expect_error(
    did_int_2x2(df, "Y_post", "Y_pre", "W", "G", g = 1, covariates = "z"),
    "no units with exposure"
  )
})
