# Shared helpers for tests: generate small reproducible test-score data.

make_paired_scores <- function(n = 300,
                               max1 = 6,
                               max2 = 5,
                               int1 = 3,
                               int2 = 2.5,
                               slope1 = 1.5,
                               slope2 = 1.3,
                               sd1 = 0.8,
                               sd2 = 0.7,
                               seed = 1) {
  set.seed(seed)
  theta <- rnorm(n)
  test1 <- pmin(pmax(round(int1 + slope1 * theta + rnorm(n, sd = sd1)), 0), max1)
  test2 <- pmin(pmax(round(int2 + slope2 * theta + rnorm(n, sd = sd2)), 0), max2)
  data.frame(test1 = test1, test2 = test2)
}
