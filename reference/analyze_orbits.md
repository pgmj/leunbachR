# Analyze Orbits for Leunbach Model

Performs orbit analysis for the Leunbach model. For each total score
(orbit), computes the expected distribution of (Test1, Test2) score
pairs and cumulative probabilities for assessing person fit.

## Usage

``` r
analyze_orbits(fit, alpha = 0.05, verbose = FALSE)
```

## Arguments

- fit:

  A leunbach_ipf object from leunbach_ipf()

- alpha:

  Significance level for critical values (default 0.05)

- verbose:

  Print detailed output

## Value

A list of class "leunbach_orbits" containing:

- orbits: Matrix of expected percentages within each total score

- left_right: Cumulative probabilities P(X ≤ x \| S = s)

- right_left: Cumulative probabilities P(X ≥ x \| S = s)

- crit_left: Critical values for left tail by total score

- crit_right: Critical values for right tail by total score

- crit_values: Combined critical values (two-tailed)

- expected_critical: Expected number of persons with significant
  differences

- variance_expected: Variance of expected critical count

- n_significant: Number of observed persons with significant differences

- orbit_df: Degrees of freedom for each orbit

## Details

For each total score S = X + Y, the orbit consists of all valid (X, Y)
pairs. The expected probability of each cell within an orbit is: P(X = x
\| S = s) = (gamma_x \* delta_y) / sigma_s

Cumulative probabilities are computed in both directions:

- Left-to-right: P(X ≤ x \| S = s) - tests if Test1 score is unusually
  low

- Right-to-left: P(X ≥ x \| S = s) - tests if Test1 score is unusually
  high
