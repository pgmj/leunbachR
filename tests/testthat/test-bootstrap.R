test_that("leunbach_bootstrap runs sequentially when parallel = FALSE", {
  skip_on_cran()
  d <- make_paired_scores(n = 200)
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)
  boot <- leunbach_bootstrap(fit, nsim = 10, parallel = FALSE, seed = 1)

  expect_s3_class(boot, "leunbach_bootstrap")
  expect_equal(boot$nsim, 10)
  expect_false(boot$parallel)
  expect_true(is.numeric(boot$p_lr))
  expect_true(boot$p_lr >= 0 && boot$p_lr <= 1)
})

test_that("get_equating_table returns clean data frame", {
  skip_on_cran()
  d <- make_paired_scores(n = 200)
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)
  boot <- leunbach_bootstrap(fit, nsim = 10, parallel = FALSE, seed = 1)
  tab <- get_equating_table(boot, direction = "1to2")
  expect_true(is.data.frame(tab))
  expect_true(all(c("log_theta", "rounded", "expected",
                    "ci_lower", "ci_upper", "see") %in% names(tab)))
})

test_that("get_equating_table rejects non-bootstrap input", {
  expect_error(get_equating_table(list()), "leunbach_bootstrap object")
})

test_that("bootstrap errors when n_cores missing under parallel=TRUE", {
  skip_if_not_installed("mirai")
  d <- make_paired_scores(n = 200)
  fit <- leunbach_ipf(d, max_score1 = 6, max_score2 = 5)
  expect_error(
    leunbach_bootstrap(fit, nsim = 5, parallel = TRUE, n_cores = NULL),
    "specify the number of cores"
  )
})
