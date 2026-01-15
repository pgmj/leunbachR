# Helper Functions for Leunbach Model

Utility functions for bootstrap and equating Generate parametric
bootstrap table

Generates a bootstrap contingency table by sampling from the fitted
model. For each total score, samples are drawn from the multinomial
distribution defined by the conditional probabilities P(X=x \| S=s).

## Usage

``` r
parametric_eq_bootstrap(
  test1_scores,
  test2_scores,
  total_scores,
  total_score_freq,
  gamma,
  delta,
  sigma,
  xmin,
  xmax,
  ymin,
  ymax
)
```

## Arguments

- test1_scores:

  Score values for Test 1

- test2_scores:

  Score values for Test 2

- total_scores:

  Total score values

- total_score_freq:

  Observed frequency for each total score

- gamma:

  Score parameters for Test 1

- delta:

  Score parameters for Test 2

- sigma:

  Total score parameters

- xmin:

  Minimum observed Test 1 score

- xmax:

  Maximum observed Test 1 score

- ymin:

  Minimum observed Test 2 score

- ymax:

  Maximum observed Test 2 score

## Value

A matrix (contingency table) of bootstrap counts
