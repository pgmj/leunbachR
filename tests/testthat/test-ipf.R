test_that("leunbach_ipf returns a leunbach_ipf object with expected fields", {
  d <- make_paired_scores(n = 300, max1 = 6, max2 = 5)
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)

  expect_s3_class(fit, "leunbach_ipf")
  expect_true(fit$converged)
  expect_equal(length(fit$gamma), 7)
  expect_equal(length(fit$delta), 6)
  expect_equal(length(fit$sigma), 12)
  expect_true(is.numeric(fit$g_sq))
  expect_true(is.numeric(fit$df))
  expect_true(fit$df >= 0)
})

test_that("leunbach_ipf rejects bad input", {
  expect_error(leunbach_ipf(list(1, 2)), "data.frame or matrix")
  expect_error(leunbach_ipf(matrix(1:9, ncol = 3)), "exactly two columns")
})

test_that("identification constraints hold: gamma[xmin] = delta[ymin] = 1", {
  d <- make_paired_scores(n = 300, max1 = 6, max2 = 5)
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)
  expect_equal(unname(fit$gamma[as.character(fit$xmin)]), 1)
  expect_equal(unname(fit$delta[as.character(fit$ymin)]), 1)
})

test_that("print and summary methods run without error", {
  d <- make_paired_scores()
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)
  expect_output(print(fit), "Leunbach Score Parameter Estimation")
  expect_output(summary(fit), "Likelihood Ratio Test")
})
