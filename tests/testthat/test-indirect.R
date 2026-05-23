make_indirect_fits <- function(seed = 1) {
  set.seed(seed)
  n <- 250
  theta1 <- rnorm(n)
  a <- pmin(pmax(round(3 + 1.5 * theta1 + rnorm(n, sd = 0.8)), 0), 6)
  b1 <- pmin(pmax(round(2.5 + 1.3 * theta1 + rnorm(n, sd = 0.7)), 0), 5)

  theta2 <- rnorm(n)
  b2 <- pmin(pmax(round(2.5 + 1.3 * theta2 + rnorm(n, sd = 0.7)), 0), 5)
  cc <- pmin(pmax(round(3 + 1.4 * theta2 + rnorm(n, sd = 0.8)), 0), 6)

  list(
    fit_ab = leunbach_ipf(data.frame(a, b1), max_score1 = 6, max_score2 = 5),
    fit_bc = leunbach_ipf(data.frame(b2, cc), max_score1 = 5, max_score2 = 6)
  )
}

test_that("leunbach_indirect_equate produces a valid table", {
  fits <- make_indirect_fits()
  ind <- leunbach_indirect_equate(fits$fit_ab, fits$fit_bc,
                                  direction_ab = "1to2",
                                  direction_bc = "1to2")
  expect_s3_class(ind, "leunbach_indirect")
  expect_true(is.data.frame(ind$equating_table))
  expect_true(all(c("source", "theta", "expected", "rounded") %in%
                    names(ind$equating_table)))
})

test_that("indirect equating respects target score bounds", {
  fits <- make_indirect_fits()
  ind <- leunbach_indirect_equate(fits$fit_ab, fits$fit_bc,
                                  direction_ab = "1to2",
                                  direction_bc = "1to2")
  expected <- ind$equating_table$expected
  expected <- expected[!is.na(expected)]
  expect_true(all(expected >= ind$target_min))
  expect_true(all(expected <= ind$target_max))
})

test_that("leunbach_indirect_equate rejects bad input", {
  expect_error(leunbach_indirect_equate(list(), list()),
               "fit_ab must be a leunbach_ipf object")
})
