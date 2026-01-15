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
