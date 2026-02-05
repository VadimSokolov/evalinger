test_that("simulate_comparison returns correct structure", {
  # Small simulation for speed
  cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 100,
                             n_looks = 2, nrep = 100, seed = 1)
  expect_s3_class(cmp, "ecomparison")
  expect_true(is.data.frame(cmp$results))
  expect_true(all(c("method", "null_rej", "alt_rej",
                     "avg_n_null", "avg_n_alt") %in% names(cmp$results)))
  expect_equal(nrow(cmp$results), 5)  # default 5 methods
})

test_that("e-value method controls Type I error", {
  cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 200,
                             n_looks = 4, nrep = 2000, seed = 42,
                             methods = "evalue")
  expect_lte(cmp$results$null_rej[1], 0.025 + 0.015)
})

test_that("naive_p inflates Type I error", {
  cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 200,
                             n_looks = 4, nrep = 2000, seed = 42,
                             methods = "naive_p")
  # Naive repeated p should inflate
  expect_gt(cmp$results$null_rej[1], 0.025)
})

test_that("simulate_comparison with subset of methods", {
  cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 100,
                             n_looks = 2, nrep = 50, seed = 1,
                             methods = c("evalue", "gs_obf"))
  expect_equal(nrow(cmp$results), 2)
})
