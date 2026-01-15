# Estimate person parameter (theta) for a given score

Estimates the person parameter theta such that the expected score equals
the observed score.

## Usage

``` r
estimate_person_parameter(
  score,
  gamma,
  score_min,
  score_max,
  method = c("optimize", "newton"),
  tol = 1e-10,
  max_iter = 500
)
```

## Arguments

- score:

  The observed score

- gamma:

  Score parameters (named vector)

- score_min:

  Minimum observed score

- score_max:

  Maximum observed score

- method:

  Optimization method: "optimize" (default) uses stats::optimize() with
  Brent's method, "newton" uses custom Newton-Raphson with bisection
  fallback

- tol:

  Convergence tolerance

- max_iter:

  Maximum iterations (only for method = "newton")

## Value

Estimated theta value
