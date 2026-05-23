# Get orbit distribution for a specific total score

Get orbit distribution for a specific total score

## Usage

``` r
get_orbit(orbits, total_score)
```

## Arguments

- orbits:

  A leunbach_orbits object

- total_score:

  The total score to examine

## Value

A data frame with the orbit distribution

## Examples

``` r
set.seed(123)
n <- 400
theta <- rnorm(n)
test1 <- pmin(pmax(round(3 + 1.5 * theta + rnorm(n, sd = 0.8)), 0), 6)
test2 <- pmin(pmax(round(2.5 + 1.3 * theta + rnorm(n, sd = 0.7)), 0), 5)
fit <- leunbach_ipf(data.frame(test1, test2),
                    max_score1 = 6, max_score2 = 5)
orb <- analyze_orbits(fit)
get_orbit(orb, total_score = 5)
#>   test1 test2 expected_pct cum_left cum_right observed
#> 1     0     5         0.00     0.00    100.00        0
#> 2     1     4         0.74     0.75    100.00        0
#> 3     2     3        32.35    33.09     99.25       16
#> 4     3     2        61.02    94.11     66.91       28
#> 5     4     1         5.85    99.97      5.89        2
#> 6     5     0         0.03   100.00      0.03        0
```
