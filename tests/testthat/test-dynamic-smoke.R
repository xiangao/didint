# Smoke test for did_int_dynamic: simulate a multi-post-period DGP with
# common adoption timing, check per-period estimates recover truth and
# that the aggregated average has tighter CIs than any single period.

simulate_dynamic <- function(N = 1200, T_post = 4,
                             direct = 1.5, spill = 0.5,
                             seed = 100) {
  set.seed(seed)
  lon <- runif(N, 0, 10); lat <- runif(N, 0, 10)
  z   <- 0.3 * lon + 0.2 * lat + rnorm(N, sd = 1)
  W <- rbinom(N, 1, plogis(-0.5 + 0.6 * z))
  dij <- as.matrix(dist(cbind(lon, lat)))
  A <- (dij < 1.5) & (dij > 0)
  share <- (A %*% W) / pmax(rowSums(A), 1)
  G <- as.integer(share > median(share))

  Y_pre  <- 0.8 * z + rnorm(N)
  # Post outcomes: same true direct + spillover at each post-period,
  # plus a small linear pre-trend in z so PT is satisfied conditional on z.
  Y_post <- matrix(0, N, T_post)
  for (k in seq_len(T_post)) {
    Y_post[, k] <- Y_pre + 0.2 * k * z + direct * W + spill * G * W + rnorm(N)
  }
  colnames(Y_post) <- paste0("Y_post_", seq_len(T_post))

  data.frame(W = W, G = G, z = z, Y_pre = Y_pre, Y_post, lon = lon, lat = lat)
}

test_that("did_int_dynamic returns one row per post period with plausible DATT", {
  df  <- simulate_dynamic(N = 1500, T_post = 4, seed = 99)
  res <- did_int_dynamic(
    data      = df,
    yname_pre = "Y_pre",
    ynames    = paste0("Y_post_", 1:4),
    treat     = "W",
    exposure  = "G",
    g         = 1,
    covariates = "z",
    trim      = 0.01
  )
  expect_s3_class(res, "didint_dynamic")
  expect_equal(nrow(res$per_period), 4)
  expect_equal(res$per_period$event_time, 0:3)

  truth <- 1.5 + 0.5
  # Each period within ~3 SE of truth
  expect_true(all(abs(res$per_period$estimate - truth) <
                  3 * res$per_period$se))
  # Aggregated SE should be smaller than the worst single-period SE
  expect_lt(res$agg$se, max(res$per_period$se))
})

test_that("event_time can be supplied explicitly", {
  df <- simulate_dynamic(N = 800, T_post = 3, seed = 101)
  res <- did_int_dynamic(
    data = df, yname_pre = "Y_pre",
    ynames = paste0("Y_post_", 1:3),
    treat = "W", exposure = "G", g = 1, covariates = "z",
    event_time = c(0L, 2L, 5L),  # skipped years
    trim = 0.01
  )
  expect_equal(res$per_period$event_time, c(0L, 2L, 5L))
})
