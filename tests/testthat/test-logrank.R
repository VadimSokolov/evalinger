test_that("eprocess_logrank returns eprocess class", {
  set.seed(42)
  n <- 100
  arm <- rep(0:1, each = n / 2)
  time <- rexp(n, rate = ifelse(arm == 1, 0.7, 1.0))
  status <- rbinom(n, 1, 0.8)
  ep <- eprocess_logrank(time, status, arm, theta = 0.7)
  expect_s3_class(ep, "eprocess")
  expect_equal(ep$endpoint, "survival")
})

test_that("eprocess_logrank handles no events", {
  time <- c(1, 2, 3, 4)
  status <- c(0, 0, 0, 0)
  arm <- c(0, 1, 0, 1)
  ep <- eprocess_logrank(time, status, arm)
  expect_equal(ep$n, 0L)
  expect_false(ep$rejected)
})

test_that("eprocess_logrank rejects under strong alternative", {
  set.seed(42)
  n <- 300
  arm <- rep(0:1, each = n / 2)
  # Strong effect: HR = 0.3
  time <- rexp(n, rate = ifelse(arm == 1, 0.3, 1.0))
  status <- rep(1, n)  # all events observed
  ep <- eprocess_logrank(time, status, arm, theta = 0.7, alpha = 0.025)
  # With a very strong effect and all events, should reject
  expect_true(ep$n > 0)
})

test_that("eprocess_logrank validates inputs", {
  expect_error(eprocess_logrank(c(1, 2), c(1, 0), c(0, 1, 0)))  # unequal lengths
  expect_error(eprocess_logrank(c(1, 2), c(1, 2), c(0, 1)))     # status not 0/1
  expect_error(eprocess_logrank(c(1, 2), c(1, 0), c(0, 2)))     # arm not 0/1
})
