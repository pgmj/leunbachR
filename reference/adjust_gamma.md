# Adjust gamma and delta parameters (ADJUST_GAMMA from Pascal)

Applies the identification constraints to the score parameters:

1.  Normalize so lowest observed score parameter = 1

2.  Apply log-linear adjustment so product of highest parameters is
    constrained

3.  Compute sigma (total score parameters)

## Usage

``` r
adjust_gamma(
  gamma,
  delta,
  test1_scores,
  test2_scores,
  xmin,
  xmax,
  ymin,
  ymax,
  sfra,
  stil,
  smin,
  smax,
  total_scores
)
```
