# Calculate Goodman-Kruskal Gamma and test statistic

Computes Goodman & Kruskal's Gamma coefficient for both observed and
expected tables, and tests whether they differ significantly. Gamma
measures the strength of association between two ordinal variables. Uses
a one-sided test following the DIGRAM implementation.

## Usage

``` r
calculate_gamma_test(observed, expected, xmin, xmax, ymin, ymax)
```

## Arguments

- observed:

  Observed contingency table

- expected:

  Expected (fitted) contingency table

- xmin:

  Minimum observed Test 1 score

- xmax:

  Maximum observed Test 1 score

- ymin:

  Minimum observed Test 2 score

- ymax:

  Maximum observed Test 2 score

## Value

A list containing:

- gamma_observed: Gamma coefficient for observed data

- gamma_expected: Gamma coefficient for expected data

- se_gamma: Standard error of gamma

- z_statistic: Z statistic (expected - observed) / SE

- p_value: One-sided p-value from standard normal
