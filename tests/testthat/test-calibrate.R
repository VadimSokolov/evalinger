test_that("edesign_binary returns correct structure", {
  des <- edesign_binary(p_C = 0.30, delta = 0.15, nrep = 1000, seed = 1)
  expect_true(is.list(des))
  expect_equal(des$p_C, 0.30)
  expect_equal(des$p_T, 0.45)
  expect_equal(des$delta, 0.15)
  expect_true(des$growth_rate > 0)
  expect_true(des$approx_stopping_time > 0)
  expect_true(des$simulated_power > 0)
  expect_true(des$simulated_type1 <= des$alpha + 0.02)
})

test_that("edesign_binary controls Type I error", {
  des <- edesign_binary(p_C = 0.30, delta = 0.15, nrep = 2000, seed = 42)
  expect_lte(des$simulated_type1, des$alpha + 0.01)
})
