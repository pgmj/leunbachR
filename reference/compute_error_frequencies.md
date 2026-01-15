# Compute frequency of bootstrap equating errors

Computes the frequency distribution of equating errors (differences
between bootstrap equated scores and observed equated scores) for each
source score.

## Usage

``` r
compute_error_frequencies(boot_rounded, observed_rounded, scores)
```

## Arguments

- boot_rounded:

  Matrix of bootstrap rounded scores (nsim x n_scores)

- observed_rounded:

  Vector of observed rounded scores

- scores:

  Vector of source scores

## Value

Matrix of error frequencies (n_scores x 5) for errors -2, -1, 0, +1, +2
