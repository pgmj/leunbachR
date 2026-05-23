# Extract equating table with bootstrap confidence intervals

Extracts equating results from a bootstrap object as a clean data frame
with log theta values, confidence intervals and standard errors.

## Usage

``` r
get_equating_table(boot, direction = c("1to2", "2to1"))
```

## Arguments

- boot:

  A leunbach_bootstrap object

- direction:

  "1to2" or "2to1"

## Value

A data frame with equating results and CIs

## Examples

``` r
# \donttest{
set.seed(123)
n <- 300
theta <- rnorm(n)
test1 <- pmin(pmax(round(3 + 1.5 * theta + rnorm(n, sd = 0.8)), 0), 6)
test2 <- pmin(pmax(round(2.5 + 1.3 * theta + rnorm(n, sd = 0.7)), 0), 5)
fit <- leunbach_ipf(data.frame(test1, test2),
                    max_score1 = 6, max_score2 = 5)
boot <- leunbach_bootstrap(fit, nsim = 25, parallel = FALSE, seed = 1)
get_equating_table(boot, direction = "1to2")
#>   Test1 log_theta rounded expected ci_lower ci_upper  see
#> 0     0   -5.0000       0     0.00     0.00     0.00 0.00
#> 1     1   -2.7954       1     0.80     0.59     1.04 0.00
#> 2     2   -1.4210       2     1.63     1.49     1.78 0.20
#> 3     3   -0.1207       3     2.54     2.40     2.68 0.48
#> 4     4    1.0766       3     3.28     3.15     3.44 0.00
#> 5     5    2.8131       4     4.12     3.99     4.34 0.00
#> 6     6    5.0000       5     5.00     5.00     5.00 0.00
# }
```
