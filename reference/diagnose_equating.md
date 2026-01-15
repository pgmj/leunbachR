# Diagnostic function for equating

Diagnoses potential issues with person parameter estimation by
displaying theta values and expected scores for each source score.

## Usage

``` r
diagnose_equating(
  fit,
  direction = c("1to2", "2to1"),
  method = c("optimize", "newton")
)
```

## Arguments

- fit:

  A leunbach_ipf object

- direction:

  Direction of equating: "1to2" or "2to1"

- method:

  Optimization method: "optimize" or "newton"
