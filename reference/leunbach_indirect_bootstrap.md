# Bootstrap for Indirect Equating

Performs parametric bootstrapping for indirect equating to assess
standard errors of the indirect equated scores.

## Usage

``` r
leunbach_indirect_bootstrap(
  fit_ab,
  fit_bc,
  direction_ab = c("1to2", "2to1"),
  direction_bc = c("1to2", "2to1"),
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

- fit_ab:

  A leunbach_ipf object for the A-B equating

- fit_bc:

  A leunbach_ipf object for the B-C equating

- direction_ab:

  Direction for A-B equating: "1to2" or "2to1"

- direction_bc:

  Direction for B-C equating: "1to2" or "2to1"

- nsim:

  Number of bootstrap samples (default: 1000)

- conf_level:

  Confidence level for intervals (default: 0.95)

- see_type:

  Type of SEE calculation: "rounded" or "expected"

- method:

  Optimization method: "optimize" (default) or "newton"

- parallel:

  Use parallel processing if mirai is available (default: TRUE)

- n_cores:

  Number of cores for parallel processing

- verbose:

  Print progress messages

- seed:

  Random seed for reproducibility

## Value

A list of class "leunbach_indirect_bootstrap" containing:

- indirect_eq: The observed indirect equating object

- Bootstrap results and standard errors

- Bootstrap p-values for LR and Gamma tests for both equatings

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

boot <- leunbach_indirect_bootstrap(fit_ab, fit_bc,
                                    direction_ab = "1to2",
                                    direction_bc = "1to2",
                                    nsim = 25, parallel = FALSE,
                                    seed = 1)
print(boot)
#> Leunbach Indirect Equating - Parametric Bootstrap Results
#> ==========================================================
#> 
#> Path: Test A -> Test B -> Test C
#> Bootstrap samples: 25 (25 valid)
#> Processing: sequential
#> Optimization method: optimize
#> SEE type: rounded scores
#> 
#> Assessment of significance by parametric bootstrapping:
#> 
#> Equating A-B (Test A -> Test B):
#>   1. Likelihood Ratio Test:
#>      Observed LR = 18.05 (df = 21)
#>      Asymptotic p-value:    p = 0.6460
#>      Bootstrap p-value:    p = 0.2400
#>   2. Goodman-Kruskal Gamma Test (one-sided):
#>      Observed Z = 0.02
#>      Asymptotic p-value:   p = 0.4925
#>      Bootstrap p-value:    p = 0.4400
#> 
#> Equating B-C (Test B -> Test C):
#>   1. Likelihood Ratio Test:
#>      Observed LR = 25.94 (df = 21)
#>      Asymptotic p-value:    p = 0.2086
#>      Bootstrap p-value:    p = 0.0400
#>   2. Goodman-Kruskal Gamma Test (one-sided):
#>      Observed Z = -0.07
#>      Asymptotic p-value:   p = 0.5262
#>      Bootstrap p-value:    p = 0.8800
#> 
#> Indirect Equating: Test A -> Test C (with 95% CI)
#> =======================================================================================
#> 
#>                                                          Frequency of bootstrap errors
#> Score      Rounded  Expected    95% CI          SEE      -2    -1     0    +1    +2   Failed%
#> --------------------------------------------------------------------------------------------------------
#>     0        0     0.00    [ 0.00,  0.00]   0.00     0.0   0.0 100.0   0.0   0.0     0.0%
#>     1        1     1.13    [ 0.80,  1.50]   0.20     0.0   0.0  96.0   4.0   0.0     0.0%
#>     2        2     2.00    [ 1.79,  2.25]   0.00     0.0   0.0 100.0   0.0   0.0     0.0%
#>     3        3     2.98    [ 2.75,  3.20]   0.00     0.0   0.0 100.0   0.0   0.0     0.0%
#>     4        4     3.88    [ 3.63,  4.09]   0.00     0.0   0.0 100.0   0.0   0.0     0.0%
#>     5        5     4.92    [ 4.56,  5.18]   0.00     0.0   0.0 100.0   0.0   0.0     0.0%
#>     6        6     6.00    [ 6.00,  6.00]   0.00     0.0   0.0 100.0   0.0   0.0     0.0%
#> --------------------------------------------------------------------------------------------------------
#> Average SEE: 0.03
# }
```
