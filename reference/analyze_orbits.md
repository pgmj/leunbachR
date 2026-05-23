# Analyze Orbits for Leunbach Model

Performs orbit analysis for the Leunbach model. For each total score
(orbit), computes the expected distribution of (Test1, Test2) score
pairs and cumulative probabilities for assessing person fit.

## Usage

``` r
analyze_orbits(fit, alpha = 0.05, verbose = FALSE)
```

## Arguments

- fit:

  A leunbach_ipf object from leunbach_ipf()

- alpha:

  Significance level for critical values (default 0.05)

- verbose:

  Print detailed output

## Value

A list of class "leunbach_orbits" containing:

- orbits: Matrix of expected percentages within each total score

- left_right: Cumulative probabilities P(X \<= x \| S = s)

- right_left: Cumulative probabilities P(X \>= x \| S = s)

- crit_left: Critical values for left tail by total score

- crit_right: Critical values for right tail by total score

- crit_values: Combined critical values (two-tailed)

- expected_critical: Expected number of persons with significant
  differences

- variance_expected: Variance of expected critical count

- n_significant: Number of observed persons with significant differences

- orbit_df: Degrees of freedom for each orbit

## Details

For each total score S = X + Y, the orbit consists of all valid (X, Y)
pairs. The expected probability of each cell within an orbit is: P(X = x
\| S = s) = (gamma_x \* delta_y) / sigma_s

Cumulative probabilities are computed in both directions:

- Left-to-right: P(X \<= x \| S = s) - tests if Test1 score is unusually
  low

- Right-to-left: P(X \>= x \| S = s) - tests if Test1 score is unusually
  high

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
print(orb)
#> Leunbach Orbit Analysis
#> =======================
#> 
#> N = 400 observations
#> Significance level:  5.0%
#> 
#> Critical levels for person fit assessment:
#> 
#>  Score  N Crit_Left Crit_Right Crit_Combined DF
#>      0 14     0.000      0.000         0.000  0
#>      1 16     0.000      0.000         0.000  1
#>      2 24     0.000      0.000         0.000  2
#>      3 43     1.887      0.000         1.887  3
#>      4 52     0.060      0.609         0.670  4
#>      5 46     0.746      0.032         0.778  5
#>      6 54     4.124      0.631         4.755  5
#>      7 57     0.877      0.065         0.942  4
#>      8 34     0.000      1.657         1.657  3
#>      9 23     0.000      0.000         0.000  2
#>     10 21     0.000      0.000         0.000  1
#>     11 16     0.000      0.000         0.000  0
#> 
#> 3 (0.8%) persons with significant differences between measurements
#> 5.2 (1.3%) expected
#> 
#> 95% Confidence interval: [0.2%, 2.4%]
#> Chi-square = 0.95, df = 1, p = 0.3298
```
