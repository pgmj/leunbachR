# Equate Scores Between Tests using Leunbach Model

Creates equating tables to convert scores from one test to another using
the estimated score parameters from the Leunbach model.

## Usage

``` r
leunbach_equate(
  fit,
  direction = c("1to2", "2to1"),
  method = c("optimize", "newton"),
  verbose = FALSE
)
```

## Arguments

- fit:

  A leunbach_ipf object from leunbach_ipf()

- direction:

  Direction of equating: "1to2" (Test1 to Test2) or "2to1" (Test2 to
  Test1)

- method:

  Optimization method for person parameter estimation: "optimize"
  (default) uses stats::optimize() with Brent's method, "newton" uses
  custom Newton-Raphson with bisection fallback

- verbose:

  Print detailed output

## Value

A list of class "leunbach_equating" containing:

- equating_table: Data frame with original scores, theta, expected
  equated scores, and rounded scores

- direction: Direction of equating

- method: Optimization method used

- fit: Original leunbach_ipf object

## Examples

``` r
# Simulate paired test score data
set.seed(123)
n <- 400
theta <- rnorm(n)
test1 <- pmin(pmax(round(3 + 1.5 * theta + rnorm(n, sd = 0.8)), 0), 6)
test2 <- pmin(pmax(round(2.5 + 1.3 * theta + rnorm(n, sd = 0.7)), 0), 5)
fit <- leunbach_ipf(data.frame(test1, test2),
                    max_score1 = 6, max_score2 = 5)

# Equate Test 1 scores onto the Test 2 metric
eq <- leunbach_equate(fit, direction = "1to2")
print(eq)
#> Leunbach Equating:  Test1 to Test2
#> Method: optimize
#> ==========================================
#> 
#>  Test1    Theta Expected_Test2 Rounded_Test2
#>      0 -5.00000           0.00             0
#>      1 -2.62714           0.87             1
#>      2 -1.24105           1.73             2
#>      3  0.04112           2.52             3
#>      4  1.42240           3.36             3
#>      5  2.77380           4.32             4
#>      6  5.00000           5.00             5
```
