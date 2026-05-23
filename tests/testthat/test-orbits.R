test_that("analyze_orbits returns expected structure", {
  d <- make_paired_scores()
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)
  orb <- analyze_orbits(fit)

  expect_s3_class(orb, "leunbach_orbits")
  expect_true(is.matrix(orb$orbits))
  expect_equal(dim(orb$orbits), c(7, 6))
  expect_true(orb$n_significant >= 0)
  expect_true(orb$expected_critical >= 0)
})

test_that("get_orbit returns rows for a valid total score", {
  d <- make_paired_scores()
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)
  orb <- analyze_orbits(fit)
  res <- get_orbit(orb, total_score = 5)
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0)
  expect_true(all(res$test1 + res$test2 == 5))
})

test_that("analyze_orbits and get_orbit reject bad input", {
  expect_error(analyze_orbits(list()), "leunbach_ipf object")
  expect_error(get_orbit(list(), 5), "leunbach_orbits object")
})
