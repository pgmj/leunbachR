# Extract equating table with bootstrap confidence intervals

Extracts equating results from a bootstrap object as a clean data frame
with log theta values, confidence intervals and standard errors.

## Usage

``` r
get_equating_table(boot, direction = c("1to2", "2to1"))
```

## Arguments

- boot:

  A leunbach_bootstrap object

- direction:

  "1to2" or "2to1"

## Value

A data frame with equating results and CIs
