test_that("eprocess_binary returns correct class", {
  set.seed(1)
  x_T <- rbinom(100, 1, 0.45)
  x_C <- rbinom(100, 1, 0.30)
  ep <- eprocess_binary(x_T, x_C, lambda = 0.31)
  expect_s3_class(ep, "eprocess")
  expect_equal(ep$n, 100)
  expect_equal(ep$endpoint, "binary")
  expect_equal(length(ep$log_evalue), 100)
})

test_that("eprocess_binary e-values are nonnegative", {
  set.seed(2)
  x_T <- rbinom(200, 1, 0.45)
  x_C <- rbinom(200, 1, 0.30)
  ep <- eprocess_binary(x_T, x_C, lambda = 0.31)
  expect_true(all(ep$evalue >= 0))
})

test_that("eprocess_binary starts at 1", {
  # The e-process at time 0 is 1, so the first evalue is (1 + lambda*D_1)
  # which is > 0. Cumulative product at n=0 is 1 by convention.
  set.seed(3)
  x_T <- rbinom(50, 1, 0.30)
  x_C <- rbinom(50, 1, 0.30)
  ep <- eprocess_binary(x_T, x_C, lambda = 0.2)
  # The product starts from 1 and multiplies: first value = 1 * (1 + 0.2*D_1)
  expect_equal(ep$evalue[1], 1 + 0.2 * (x_T[1] - x_C[1]))
})

test_that("eprocess_binary controls Type I error (Monte Carlo)", {
  # Under the null, the rejection rate should be <= alpha
  set.seed(42)
  nrep <- 5000
  alpha <- 0.025
  n <- 100
  rejected <- 0
  for (i in seq_len(nrep)) {
    x_T <- rbinom(n, 1, 0.30)
    x_C <- rbinom(n, 1, 0.30)
    ep <- eprocess_binary(x_T, x_C, lambda = 0.31, alpha = alpha)
    if (ep$rejected) rejected <- rejected + 1
  }
  rej_rate <- rejected / nrep
  # Should be well below alpha (Ville's inequality)
  expect_lte(rej_rate, alpha + 0.01)  # allow small MC slack
})

test_that("eprocess_binary rejects under strong alternative", {
  set.seed(10)
  x_T <- rbinom(200, 1, 0.60)
  x_C <- rbinom(200, 1, 0.30)
  ep <- eprocess_binary(x_T, x_C, lambda = 0.31)
  expect_true(ep$rejected)
})

test_that("eprocess_binary uses GROW lambda when lambda is NULL", {
  set.seed(5)
  x_T <- rbinom(100, 1, 0.45)
  x_C <- rbinom(100, 1, 0.30)
  ep <- eprocess_binary(x_T, x_C, lambda = NULL,
                        p_T_design = 0.45, p_C_design = 0.30)
  expect_equal(ep$lambda, grow_lambda(0.45, 0.30))
})

test_that("eprocess_binary rejects bad inputs", {
  expect_error(eprocess_binary(c(0, 1, 2), c(0, 1, 0), lambda = 0.3))
  expect_error(eprocess_binary(c(0, 1), c(0, 1, 0), lambda = 0.3))
})
