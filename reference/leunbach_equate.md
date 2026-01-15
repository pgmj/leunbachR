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

- equating_table: Data frame with original scores, expected equated
  scores, and rounded scores

- direction: Direction of equating

- method: Optimization method used

- fit: Original leunbach_ipf object
