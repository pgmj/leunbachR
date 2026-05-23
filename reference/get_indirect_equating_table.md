# Get indirect equating table

Get indirect equating table

## Usage

``` r
get_indirect_equating_table(boot)
```

## Arguments

- boot:

  A leunbach_indirect_bootstrap object

## Value

A data frame with equating results and CIs

## Examples

``` r
# \donttest{
set.seed(123)
n <- 300
theta1 <- rnorm(n)
a <- pmin(pmax(round(3 + 1.5 * theta1 + rnorm(n, sd = 0.8)), 0), 6)
b1 <- pmin(pmax(round(2.5 + 1.3 * theta1 + rnorm(n, sd = 0.7)), 0), 5)
theta2 <- rnorm(n)
b2 <- pmin(pmax(round(2.5 + 1.3 * theta2 + rnorm(n, sd = 0.7)), 0), 5)
cc <- pmin(pmax(round(3 + 1.4 * theta2 + rnorm(n, sd = 0.8)), 0), 6)
fit_ab <- leunbach_ipf(data.frame(a, b1), max_score1 = 6, max_score2 = 5)
fit_bc <- leunbach_ipf(data.frame(b2, cc), max_score1 = 5, max_score2 = 6)
boot <- leunbach_indirect_bootstrap(fit_ab, fit_bc, nsim = 25,
                                    parallel = FALSE, seed = 1)
get_indirect_equating_table(boot)
#>   Test A expected ci_lower ci_upper rounded see pct_failed
#> 0      0     0.00     0.00     0.00       0 0.0          0
#> 1      1     1.13     0.80     1.50       1 0.2          0
#> 2      2     2.00     1.79     2.25       2 0.0          0
#> 3      3     2.98     2.75     3.20       3 0.0          0
#> 4      4     3.88     3.63     4.09       4 0.0          0
#> 5      5     4.92     4.56     5.18       5 0.0          0
#> 6      6     6.00     6.00     6.00       6 0.0          0
# }
```
