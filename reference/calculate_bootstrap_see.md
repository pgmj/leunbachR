# Calculate Bootstrap Standard Error of Equating

Computes the standard error of equating from bootstrap samples using the
standard deviation of the bootstrap distribution.

## Usage

``` r
calculate_bootstrap_see(boot_matrix)
```

## Arguments

- boot_matrix:

  Matrix of bootstrap equated scores (nsim x n_scores)

## Value

Named vector of SEE values for each score
