test_that("leunbach_equate returns expected structure for both directions", {
  d <- make_paired_scores()
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)

  eq12 <- leunbach_equate(fit, direction = "1to2")
  expect_s3_class(eq12, "leunbach_equating")
  expect_true(is.data.frame(eq12$equating_table))
  expect_equal(ncol(eq12$equating_table), 4)
  expect_equal(eq12$source_name, "Test1")
  expect_equal(eq12$target_name, "Test2")

  eq21 <- leunbach_equate(fit, direction = "2to1")
  expect_equal(eq21$source_name, "Test2")
  expect_equal(eq21$target_name, "Test1")
})

test_that("equated scores are monotone within the observed range", {
  d <- make_paired_scores()
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)
  eq <- leunbach_equate(fit, direction = "1to2")
  tab <- eq$equating_table
  valid <- tab[, 1] >= fit$xmin & tab[, 1] <= fit$xmax
  expected <- tab[valid, 3]
  expect_false(any(is.na(expected)))
  expect_true(all(diff(expected) >= -1e-8))
})

test_that("optimize and newton methods give similar results", {
  d <- make_paired_scores()
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)
  eq_opt <- leunbach_equate(fit, direction = "1to2", method = "optimize")
  eq_new <- leunbach_equate(fit, direction = "1to2", method = "newton")
  diffs <- abs(eq_opt$equating_table[, 3] - eq_new$equating_table[, 3])
  expect_true(max(diffs, na.rm = TRUE) < 0.05)
})

test_that("leunbach_equate rejects non-leunbach_ipf input", {
  expect_error(leunbach_equate(list(a = 1)), "leunbach_ipf object")
})
