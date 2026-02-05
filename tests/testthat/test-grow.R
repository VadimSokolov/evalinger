test_that("grow_lambda returns correct GROW fraction", {
  # For p_T = 0.45, p_C = 0.30:
  # a = 0.45 * 0.70 = 0.315, b = 0.55 * 0.30 = 0.165
  # lambda* = (0.315 - 0.165) / (0.315 + 0.165) = 0.15 / 0.48 = 0.3125
  lam <- grow_lambda(0.45, 0.30)
  expect_equal(lam, 0.3125, tolerance = 1e-10)
})

test_that("grow_lambda rejects invalid inputs", {
  expect_error(grow_lambda(0.3, 0.3))   # p_T must be > p_C
  expect_error(grow_lambda(0.2, 0.3))   # p_T < p_C
  expect_error(grow_lambda(1.1, 0.3))   # p_T > 1
  expect_error(grow_lambda(0.5, -0.1))  # p_C < 0
})

test_that("expected_growth_rate is positive under alternative", {
  lam <- grow_lambda(0.45, 0.30)
  g <- expected_growth_rate(lam, 0.45, 0.30)
  expect_gt(g, 0)
})

test_that("expected_growth_rate is negative under null", {
  # Under null p_T = p_C = 0.30, any lambda > 0 yields negative growth
  g <- expected_growth_rate(0.3, 0.30, 0.30)
  expect_lt(g, 0)
})

test_that("GROW lambda maximizes growth rate", {
  lam_opt <- grow_lambda(0.45, 0.30)
  g_opt <- expected_growth_rate(lam_opt, 0.45, 0.30)
  # Slightly perturbed lambdas should give lower growth
  g_lo <- expected_growth_rate(lam_opt - 0.05, 0.45, 0.30)
  g_hi <- expected_growth_rate(lam_opt + 0.05, 0.45, 0.30)
  expect_gt(g_opt, g_lo)
  expect_gt(g_opt, g_hi)
})

test_that("expected_stopping_time is finite and reasonable", {
  lam <- grow_lambda(0.45, 0.30)
  tau <- expected_stopping_time(lam, 0.45, 0.30, alpha = 0.025)
  expect_true(is.finite(tau))
  expect_gt(tau, 10)
  expect_lt(tau, 1000)
})

test_that("grow_lambda_grid returns a data frame", {
  df <- grow_lambda_grid(p_C = 0.30, delta_grid = c(0.10, 0.15))
  expect_s3_class(df, "data.frame")
  expect_true(all(c("delta", "lambda", "growth_rate", "expected_n") %in% names(df)))
  expect_gt(nrow(df), 0)
})
