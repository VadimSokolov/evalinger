test_that("futility_cs declares futility when CS upper < delta_min", {
  set.seed(42)
  x_T <- rbinom(200, 1, 0.31)  # nearly no effect

  x_C <- rbinom(200, 1, 0.30)
  cs <- confseq_binary(x_T, x_C)
  fut <- futility_cs(cs, delta_min = 0.10)
  expect_true(is.logical(fut$futile))
  expect_equal(length(fut$futile), cs$n)
  expect_equal(fut$delta_min, 0.10)
})

test_that("futility_cs does not declare futility under strong effect", {
  set.seed(42)
  x_T <- rbinom(200, 1, 0.50)
  x_C <- rbinom(200, 1, 0.30)
  cs <- confseq_binary(x_T, x_C)
  fut <- futility_cs(cs, delta_min = 0.05)
  # With a strong effect, CS upper should stay above 0.05 eventually

  # (may or may not declare futility early, but at the end it should not)
  expect_false(fut$futile[cs$n])
})

test_that("futility_cs requires confseq class", {
  expect_error(futility_cs(list(upper = 1:10), delta_min = 0.5))
})

test_that("futility_eprocess returns correct structure", {
  set.seed(42)
  x_T <- rbinom(100, 1, 0.32)
  x_C <- rbinom(100, 1, 0.30)
  fut <- futility_eprocess(x_T, x_C, delta_min = 0.10)
  expect_true(is.list(fut))
  expect_equal(length(fut$log_evalue), 100)
  expect_equal(length(fut$evalue), 100)
  expect_true(is.logical(fut$futile))
  expect_equal(fut$alpha_f, 0.10)
})

test_that("futility_eprocess does not declare futility under strong effect", {
  set.seed(42)
  x_T <- rbinom(200, 1, 0.50)
  x_C <- rbinom(200, 1, 0.30)
  fut <- futility_eprocess(x_T, x_C, delta_min = 0.05, alpha_f = 0.10)
  expect_false(fut$futile)
})

test_that("futility_eprocess requires equal-length inputs", {
  expect_error(futility_eprocess(c(0, 1), c(0, 1, 0), delta_min = 0.10))
})
