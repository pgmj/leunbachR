# Parametric Bootstrap for Leunbach Model with Standard Error of Equating

Performs parametric bootstrapping to assess significance of tests and
compute standard errors of equating (SEE) with confidence intervals.
Supports parallel processing using the mirai package.

## Usage

``` r
leunbach_bootstrap(
  fit,
  nsim = 1000,
  conf_level = 0.95,
  see_type = c("rounded", "expected"),
  method = c("optimize", "newton"),
  parallel = TRUE,
  n_cores = NULL,
  verbose = FALSE,
  seed = NULL
)
```

## Arguments

- fit:

  A leunbach_ipf object from leunbach_ipf()

- nsim:

  Number of bootstrap samples (default: 1000)

- conf_level:

  Confidence level for intervals (default: 0.95)

- see_type:

  Type of SEE calculation: "rounded" uses rounded (integer) scores,
  "expected" uses continuous expected scores

- method:

  Optimization method for person parameter estimation: "optimize"
  (default) uses stats::optimize() with Brent's method, "newton" uses
  custom Newton-Raphson with bisection fallback

- parallel:

  Use parallel processing if mirai package is available (default: TRUE)

- n_cores:

  Number of cores to use for parallel processing. Default NULL uses all
  available cores minus one.

- verbose:

  Print progress messages

- seed:

  Random seed for reproducibility (optional)

## Value

A list of class "leunbach_bootstrap" containing bootstrap results and
standard errors of equating

## Examples

``` r
# \donttest{
# Simulate paired test score data
set.seed(123)
n <- 300
theta <- rnorm(n)
test1 <- pmin(pmax(round(3 + 1.5 * theta + rnorm(n, sd = 0.8)), 0), 6)
test2 <- pmin(pmax(round(2.5 + 1.3 * theta + rnorm(n, sd = 0.7)), 0), 5)
fit <- leunbach_ipf(data.frame(test1, test2),
                    max_score1 = 6, max_score2 = 5)

# Sequential bootstrap with few iterations for the example
boot <- leunbach_bootstrap(fit, nsim = 25, parallel = FALSE, seed = 1)
print(boot)
#> Leunbach Model - Parametric Bootstrap Results
#> ==============================================
#> 
#> Bootstrap samples: 25 (25 valid)
#> Processing: sequential
#> Optimization method: optimize
#> SEE type: rounded scores
#> 
#> Assessment of significance by parametric bootstrapping:
#> 
#> 1. Likelihood Ratio Test:
#>    Observed LR = 18.05 (df = 21)
#>    Asymptotic p-value:   p = 0.6460
#>    Bootstrap p-value:   p = 0.2400
#> 
#> 2. Goodman-Kruskal Gamma Test (one-sided):
#>    Observed Z = 0.02
#>    Asymptotic p-value:   p = 0.4925
#>    Bootstrap p-value:   p = 0.4400
#> 
#> Equating Test1 to Test2 (with 95% CI)
#> ============================================================================
#> 
#>                                                                 Frequency of bootstrap errors
#> Score      Theta  Rounded  Expected    95% CI          SEE      -2    -1     0    +1    +2
#> --------------------------------------------------------------------------------------------------
#>     0   -5.00000        0     0.00    [ 0.00,  0.00]   0.00     0.0   0.0 100.0   0.0   0.0
#>     1   -2.79540        1     0.80    [ 0.59,  1.04]   0.00     0.0   0.0 100.0   0.0   0.0
#>     2   -1.42096        2     1.63    [ 1.49,  1.78]   0.20     0.0   4.0  96.0   0.0   0.0
#>     3   -0.12068        3     2.54    [ 2.40,  2.68]   0.48     0.0  32.0  68.0   0.0   0.0
#>     4    1.07660        3     3.28    [ 3.15,  3.44]   0.00     0.0   0.0 100.0   0.0   0.0
#>     5    2.81314        4     4.12    [ 3.99,  4.34]   0.00     0.0   0.0 100.0   0.0   0.0
#>     6    5.00000        5     5.00    [ 5.00,  5.00]   0.00     0.0   0.0 100.0   0.0   0.0
#> --------------------------------------------------------------------------------------------------
#> Average SEE: 0.10
#> 
#> Equating Test2 to Test1 (with 95% CI)
#> ============================================================================
#> 
#>                                                                 Frequency of bootstrap errors
#> Score      Theta  Rounded  Expected    95% CI          SEE      -2    -1     0    +1    +2
#> --------------------------------------------------------------------------------------------------
#>     0   -5.00000        0     0.00    [ 0.00,  0.00]   0.00     0.0   0.0 100.0   0.0   0.0
#>     1   -2.41857        1     1.27    [ 0.94,  1.51]   0.20     0.0   0.0  96.0   4.0   0.0
#>     2   -0.90032        2     2.38    [ 2.23,  2.53]   0.28     0.0   0.0  92.0   8.0   0.0
#>     3    0.59069        4     3.61    [ 3.43,  3.79]   0.41     0.0  20.0  80.0   0.0   0.0
#>     4    2.53659        5     4.87    [ 4.63,  5.01]   0.00     0.0   0.0 100.0   0.0   0.0
#>     5    5.00000        6     6.00    [ 6.00,  6.00]   0.00     0.0   0.0 100.0   0.0   0.0
#> --------------------------------------------------------------------------------------------------
#> Average SEE: 0.15
# }
```
