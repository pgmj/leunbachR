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
