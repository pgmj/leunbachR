# Estimate Leunbach Score Parameters using Iterative Proportional Fitting

Estimates score parameters for Leunbach test equating based on the power
series distribution framework. The total score distribution in Rasch
models follows a power series distribution where score parameters are
defined by generalized symmetric functions of the item parameters.

## Usage

``` r
leunbach_ipf(
  data,
  max_score1 = NULL,
  max_score2 = NULL,
  max_iter = 1000,
  tol = 1e-10,
  verbose = FALSE
)
```

## Arguments

- data:

  A data. frame or matrix with two columns: Column 1: Sum scores from
  Test 1 Column 2: Sum scores from Test 2

- max_score1:

  Maximum possible score for Test 1 (default: max observed)

- max_score2:

  Maximum possible score for Test 2 (default: max observed)

- max_iter:

  Maximum number of iterations for IPF

- tol:

  Convergence tolerance

- verbose:

  Print iteration progress

## Value

A list of class "leunbach_ipf" containing:

- gamma: Score parameters for Test1 (generalized symmetric functions)

- delta: Score parameters for Test2 (generalized symmetric functions)

- sigma: Score parameters for total score (Test1 + Test2)

- log_gamma: Log of gamma parameters

- log_delta: Log of delta parameters

- log_sigma: Log of sigma parameters

- fitted: Fitted frequency table

- observed: Observed frequency table

- iterations: Number of iterations to convergence

- converged: Logical indicating convergence

- test1_scores: Score values for Test1

- test2_scores: Score values for Test2

- total_scores: Score values for total score

- g_sq: Likelihood ratio test statistic

- chi_sq: Pearson chi-square statistic

- df: Degrees of freedom

- p_value: P-value for likelihood ratio test

## Details

The model is based on the power series distribution (PSD) for sum scores
in Rasch models. For a test with k polytomous items, the probability of
obtaining sum score r is:

\$\$P(R = r \| \theta) = \frac{\gamma_r \theta^r}{\sum_s \gamma_s
\theta^s}\$\$

where \\\gamma_r\\ is the generalized symmetric function of order r of
the item parameters, and \\\theta = exp(\beta)\\ is the person
parameter.

The model is over-parameterized, so we apply identification constraints:

- \\\gamma\_{min}\\ = 1 (lowest observed score parameter fixed at 1)

- \\\delta\_{min}\\ = 1 (lowest observed score parameter fixed at 1)

- The highest score parameters are constrained so that their product
  equals the product of all score parameters (Rasch-style constraint)

The sigma parameters for the total score (S = X + Y) are computed as:
\$\$\sigma_s = \sum\_{x+y=s} \gamma_x \delta_y\$\$

## References

Leunbach, G. (1976). A probabilistic measurement model for assessing
whether two tests measure the same personal factor. Copenhagen: Danish
Institute for Educational Research.

Adroher, N. D., Kreiner, S., Young, C., Mills, R., & Tennant, A. (2019).
Test equating sleep scales: Applying the Leunbach’s model. BMC Medical
Research Methodology, 19(1), 141.
<https://doi.org/10.1186/s12874-019-0768-y>

Kreiner, S. (2007). Validity and objectivity: Reflections on the role
and nature of Rasch models. Nordic Psychology, 59(3), 268-298.
<https://doi.org/10.1027/1901-2276.59.3.268>

## Examples

``` r
# Simulate test score data
set.seed(123)
n <- 500
theta <- rnorm(n)
test1 <- pmin(pmax(round(3 + 1.5 * theta + rnorm(n, sd = 0.8)), 0), 6)
test2 <- pmin(pmax(round(2.5 + 1.3 * theta + rnorm(n, sd = 0.7)), 0), 5)
score_data <- data.frame(test1 = test1, test2 = test2)

# Estimate parameters (specifying max possible scores)
fit <- leunbach_ipf(score_data, max_score1 = 6, max_score2 = 5, verbose = TRUE)
#> Leunbach Test Equating - Score Parameter Estimation
#> ====================================================
#> 
#> Power Series Distribution Model with Generalized Symmetric Functions
#> 
#> Data summary:
#>   N = 500 observations
#>   Test1: scores 0 to 6 (observed:  0 to 6)
#>   Test2: scores 0 to 5 (observed:  0 to 5)
#>   Total:  scores 0 to 11 (observed: 0 to 11)
#> 
#> Estimating score parameters via IPF...
#> 
#>   Iteration 1: max param change = 3.28e+00 (gamma:  3.23e+00, delta: 3.28e+00)
#>   Iteration 2: max param change = 4.95e+00 (gamma:  4.95e+00, delta: 3.42e+00)
#>   Iteration 3: max param change = 7.06e+00 (gamma:  7.06e+00, delta: 3.66e+00)
#>   Iteration 4: max param change = 9.25e+00 (gamma:  9.25e+00, delta: 3.91e+00)
#>   Iteration 5: max param change = 1.14e+01 (gamma:  1.14e+01, delta: 4.17e+00)
#>   Iteration 6: max param change = 1.36e+01 (gamma:  1.36e+01, delta: 4.41e+00)
#>   Iteration 7: max param change = 1.57e+01 (gamma:  1.57e+01, delta: 4.61e+00)
#>   Iteration 8: max param change = 1.77e+01 (gamma:  1.77e+01, delta: 4.79e+00)
#>   Iteration 9: max param change = 1.96e+01 (gamma:  1.96e+01, delta: 4.94e+00)
#>   Iteration 10: max param change = 2.13e+01 (gamma:  2.13e+01, delta: 5.06e+00)
#>   Iteration 11: max param change = 2.29e+01 (gamma:  2.29e+01, delta: 5.15e+00)
#>   Iteration 12: max param change = 2.44e+01 (gamma:  2.44e+01, delta: 5.21e+00)
#>   Iteration 13: max param change = 2.57e+01 (gamma:  2.57e+01, delta: 5.25e+00)
#>   Iteration 14: max param change = 2.68e+01 (gamma:  2.68e+01, delta: 5.27e+00)
#>   Iteration 15: max param change = 2.78e+01 (gamma:  2.78e+01, delta: 5.27e+00)
#>   Iteration 16: max param change = 2.87e+01 (gamma:  2.87e+01, delta: 5.25e+00)
#>   Iteration 17: max param change = 2.94e+01 (gamma:  2.94e+01, delta: 5.21e+00)
#>   Iteration 18: max param change = 2.99e+01 (gamma:  2.99e+01, delta: 5.16e+00)
#>   Iteration 19: max param change = 3.03e+01 (gamma:  3.03e+01, delta: 5.10e+00)
#>   Iteration 20: max param change = 3.06e+01 (gamma:  3.06e+01, delta: 5.02e+00)
#>   Iteration 21: max param change = 3.08e+01 (gamma:  3.08e+01, delta: 4.94e+00)
#>   Iteration 22: max param change = 3.09e+01 (gamma:  3.09e+01, delta: 4.85e+00)
#>   Iteration 23: max param change = 3.09e+01 (gamma:  3.09e+01, delta: 4.75e+00)
#>   Iteration 24: max param change = 3.08e+01 (gamma:  3.08e+01, delta: 4.64e+00)
#>   Iteration 25: max param change = 3.05e+01 (gamma:  3.05e+01, delta: 4.53e+00)
#>   Iteration 26: max param change = 3.03e+01 (gamma:  3.03e+01, delta: 4.42e+00)
#>   Iteration 27: max param change = 2.99e+01 (gamma:  2.99e+01, delta: 4.30e+00)
#>   Iteration 28: max param change = 2.95e+01 (gamma:  2.95e+01, delta: 4.18e+00)
#>   Iteration 29: max param change = 2.90e+01 (gamma:  2.90e+01, delta: 4.06e+00)
#>   Iteration 30: max param change = 2.85e+01 (gamma:  2.85e+01, delta: 3.94e+00)
#>   Iteration 31: max param change = 2.80e+01 (gamma:  2.80e+01, delta: 3.82e+00)
#>   Iteration 32: max param change = 2.74e+01 (gamma:  2.74e+01, delta: 3.69e+00)
#>   Iteration 33: max param change = 2.68e+01 (gamma:  2.68e+01, delta: 3.57e+00)
#>   Iteration 34: max param change = 2.61e+01 (gamma:  2.61e+01, delta: 3.45e+00)
#>   Iteration 35: max param change = 2.55e+01 (gamma:  2.55e+01, delta: 3.33e+00)
#>   Iteration 36: max param change = 2.48e+01 (gamma:  2.48e+01, delta: 3.22e+00)
#>   Iteration 37: max param change = 2.41e+01 (gamma:  2.41e+01, delta: 3.10e+00)
#>   Iteration 38: max param change = 2.34e+01 (gamma:  2.34e+01, delta: 2.99e+00)
#>   Iteration 39: max param change = 2.27e+01 (gamma:  2.27e+01, delta: 2.88e+00)
#>   Iteration 40: max param change = 2.20e+01 (gamma:  2.20e+01, delta: 2.77e+00)
#>   Iteration 41: max param change = 2.13e+01 (gamma:  2.13e+01, delta: 2.67e+00)
#>   Iteration 42: max param change = 2.06e+01 (gamma:  2.06e+01, delta: 2.56e+00)
#>   Iteration 43: max param change = 2.00e+01 (gamma:  2.00e+01, delta: 2.46e+00)
#>   Iteration 44: max param change = 1.93e+01 (gamma:  1.93e+01, delta: 2.37e+00)
#>   Iteration 45: max param change = 1.86e+01 (gamma:  1.86e+01, delta: 2.27e+00)
#>   Iteration 46: max param change = 1.80e+01 (gamma:  1.80e+01, delta: 2.18e+00)
#>   Iteration 47: max param change = 1.73e+01 (gamma:  1.73e+01, delta: 2.09e+00)
#>   Iteration 48: max param change = 1.67e+01 (gamma:  1.67e+01, delta: 2.01e+00)
#>   Iteration 49: max param change = 1.61e+01 (gamma:  1.61e+01, delta: 1.92e+00)
#>   Iteration 50: max param change = 1.55e+01 (gamma:  1.55e+01, delta: 1.84e+00)
#>   Iteration 51: max param change = 1.49e+01 (gamma:  1.49e+01, delta: 1.77e+00)
#>   Iteration 52: max param change = 1.43e+01 (gamma:  1.43e+01, delta: 1.69e+00)
#>   Iteration 53: max param change = 1.37e+01 (gamma:  1.37e+01, delta: 1.62e+00)
#>   Iteration 54: max param change = 1.32e+01 (gamma:  1.32e+01, delta: 1.55e+00)
#>   Iteration 55: max param change = 1.27e+01 (gamma:  1.27e+01, delta: 1.48e+00)
#>   Iteration 56: max param change = 1.22e+01 (gamma:  1.22e+01, delta: 1.42e+00)
#>   Iteration 57: max param change = 1.17e+01 (gamma:  1.17e+01, delta: 1.36e+00)
#>   Iteration 58: max param change = 1.12e+01 (gamma:  1.12e+01, delta: 1.30e+00)
#>   Iteration 59: max param change = 1.07e+01 (gamma:  1.07e+01, delta: 1.24e+00)
#>   Iteration 60: max param change = 1.03e+01 (gamma:  1.03e+01, delta: 1.19e+00)
#>   Iteration 61: max param change = 9.84e+00 (gamma:  9.84e+00, delta: 1.14e+00)
#>   Iteration 62: max param change = 9.43e+00 (gamma:  9.43e+00, delta: 1.09e+00)
#>   Iteration 63: max param change = 9.03e+00 (gamma:  9.03e+00, delta: 1.04e+00)
#>   Iteration 64: max param change = 8.64e+00 (gamma:  8.64e+00, delta: 9.91e-01)
#>   Iteration 65: max param change = 8.27e+00 (gamma:  8.27e+00, delta: 9.47e-01)
#>   Iteration 66: max param change = 7.92e+00 (gamma:  7.92e+00, delta: 9.04e-01)
#>   Iteration 67: max param change = 7.57e+00 (gamma:  7.57e+00, delta: 8.64e-01)
#>   Iteration 68: max param change = 7.24e+00 (gamma:  7.24e+00, delta: 8.25e-01)
#>   Iteration 69: max param change = 6.93e+00 (gamma:  6.93e+00, delta: 7.87e-01)
#>   Iteration 70: max param change = 6.62e+00 (gamma:  6.62e+00, delta: 7.52e-01)
#>   Iteration 71: max param change = 6.33e+00 (gamma:  6.33e+00, delta: 7.17e-01)
#>   Iteration 72: max param change = 6.05e+00 (gamma:  6.05e+00, delta: 6.85e-01)
#>   Iteration 73: max param change = 5.78e+00 (gamma:  5.78e+00, delta: 6.53e-01)
#>   Iteration 74: max param change = 5.52e+00 (gamma:  5.52e+00, delta: 6.24e-01)
#>   Iteration 75: max param change = 5.28e+00 (gamma:  5.28e+00, delta: 5.95e-01)
#>   Iteration 76: max param change = 5.04e+00 (gamma:  5.04e+00, delta: 5.68e-01)
#>   Iteration 77: max param change = 4.81e+00 (gamma:  4.81e+00, delta: 5.42e-01)
#>   Iteration 78: max param change = 4.60e+00 (gamma:  4.60e+00, delta: 5.17e-01)
#>   Iteration 79: max param change = 4.39e+00 (gamma:  4.39e+00, delta: 4.93e-01)
#>   Iteration 80: max param change = 4.19e+00 (gamma:  4.19e+00, delta: 4.70e-01)
#>   Iteration 81: max param change = 4.00e+00 (gamma:  4.00e+00, delta: 4.48e-01)
#>   Iteration 82: max param change = 3.82e+00 (gamma:  3.82e+00, delta: 4.28e-01)
#>   Iteration 83: max param change = 3.64e+00 (gamma:  3.64e+00, delta: 4.08e-01)
#>   Iteration 84: max param change = 3.48e+00 (gamma:  3.48e+00, delta: 3.89e-01)
#>   Iteration 85: max param change = 3.32e+00 (gamma:  3.32e+00, delta: 3.71e-01)
#>   Iteration 86: max param change = 3.16e+00 (gamma:  3.16e+00, delta: 3.54e-01)
#>   Iteration 87: max param change = 3.02e+00 (gamma:  3.02e+00, delta: 3.37e-01)
#>   Iteration 88: max param change = 2.88e+00 (gamma:  2.88e+00, delta: 3.21e-01)
#>   Iteration 89: max param change = 2.75e+00 (gamma:  2.75e+00, delta: 3.06e-01)
#>   Iteration 90: max param change = 2.62e+00 (gamma:  2.62e+00, delta: 2.92e-01)
#>   Iteration 91: max param change = 2.50e+00 (gamma:  2.50e+00, delta: 2.78e-01)
#>   Iteration 92: max param change = 2.38e+00 (gamma:  2.38e+00, delta: 2.65e-01)
#>   Iteration 93: max param change = 2.27e+00 (gamma:  2.27e+00, delta: 2.53e-01)
#>   Iteration 94: max param change = 2.17e+00 (gamma:  2.17e+00, delta: 2.41e-01)
#>   Iteration 95: max param change = 2.07e+00 (gamma:  2.07e+00, delta: 2.30e-01)
#>   Iteration 96: max param change = 1.97e+00 (gamma:  1.97e+00, delta: 2.19e-01)
#>   Iteration 97: max param change = 1.88e+00 (gamma:  1.88e+00, delta: 2.09e-01)
#>   Iteration 98: max param change = 1.79e+00 (gamma:  1.79e+00, delta: 1.99e-01)
#>   Iteration 99: max param change = 1.71e+00 (gamma:  1.71e+00, delta: 1.90e-01)
#>   Iteration 100: max param change = 1.63e+00 (gamma:  1.63e+00, delta: 1.81e-01)
#>   Iteration 101: max param change = 1.55e+00 (gamma:  1.55e+00, delta: 1.72e-01)
#>   Iteration 102: max param change = 1.48e+00 (gamma:  1.48e+00, delta: 1.64e-01)
#>   Iteration 103: max param change = 1.41e+00 (gamma:  1.41e+00, delta: 1.56e-01)
#>   Iteration 104: max param change = 1.34e+00 (gamma:  1.34e+00, delta: 1.49e-01)
#>   Iteration 105: max param change = 1.28e+00 (gamma:  1.28e+00, delta: 1.42e-01)
#>   Iteration 106: max param change = 1.22e+00 (gamma:  1.22e+00, delta: 1.35e-01)
#>   Iteration 107: max param change = 1.16e+00 (gamma:  1.16e+00, delta: 1.29e-01)
#>   Iteration 108: max param change = 1.11e+00 (gamma:  1.11e+00, delta: 1.23e-01)
#>   Iteration 109: max param change = 1.06e+00 (gamma:  1.06e+00, delta: 1.17e-01)
#>   Iteration 110: max param change = 1.01e+00 (gamma:  1.01e+00, delta: 1.12e-01)
#>   Iteration 111: max param change = 9.61e-01 (gamma:  9.61e-01, delta: 1.06e-01)
#>   Iteration 112: max param change = 9.15e-01 (gamma:  9.15e-01, delta: 1.01e-01)
#>   Iteration 113: max param change = 8.72e-01 (gamma:  8.72e-01, delta: 9.65e-02)
#>   Iteration 114: max param change = 8.31e-01 (gamma:  8.31e-01, delta: 9.19e-02)
#>   Iteration 115: max param change = 7.92e-01 (gamma:  7.92e-01, delta: 8.76e-02)
#>   Iteration 116: max param change = 7.55e-01 (gamma:  7.55e-01, delta: 8.35e-02)
#>   Iteration 117: max param change = 7.19e-01 (gamma:  7.19e-01, delta: 7.95e-02)
#>   Iteration 118: max param change = 6.85e-01 (gamma:  6.85e-01, delta: 7.58e-02)
#>   Iteration 119: max param change = 6.53e-01 (gamma:  6.53e-01, delta: 7.22e-02)
#>   Iteration 120: max param change = 6.22e-01 (gamma:  6.22e-01, delta: 6.88e-02)
#>   Iteration 121: max param change = 5.93e-01 (gamma:  5.93e-01, delta: 6.55e-02)
#>   Iteration 122: max param change = 5.65e-01 (gamma:  5.65e-01, delta: 6.24e-02)
#>   Iteration 123: max param change = 5.38e-01 (gamma:  5.38e-01, delta: 5.95e-02)
#>   Iteration 124: max param change = 5.13e-01 (gamma:  5.13e-01, delta: 5.66e-02)
#>   Iteration 125: max param change = 4.89e-01 (gamma:  4.89e-01, delta: 5.40e-02)
#>   Iteration 126: max param change = 4.65e-01 (gamma:  4.65e-01, delta: 5.14e-02)
#>   Iteration 127: max param change = 4.43e-01 (gamma:  4.43e-01, delta: 4.90e-02)
#>   Iteration 128: max param change = 4.23e-01 (gamma:  4.23e-01, delta: 4.67e-02)
#>   Iteration 129: max param change = 4.03e-01 (gamma:  4.03e-01, delta: 4.44e-02)
#>   Iteration 130: max param change = 3.83e-01 (gamma:  3.83e-01, delta: 4.23e-02)
#>   Iteration 131: max param change = 3.65e-01 (gamma:  3.65e-01, delta: 4.03e-02)
#>   Iteration 132: max param change = 3.48e-01 (gamma:  3.48e-01, delta: 3.84e-02)
#>   Iteration 133: max param change = 3.32e-01 (gamma:  3.32e-01, delta: 3.66e-02)
#>   Iteration 134: max param change = 3.16e-01 (gamma:  3.16e-01, delta: 3.49e-02)
#>   Iteration 135: max param change = 3.01e-01 (gamma:  3.01e-01, delta: 3.32e-02)
#>   Iteration 136: max param change = 2.87e-01 (gamma:  2.87e-01, delta: 3.16e-02)
#>   Iteration 137: max param change = 2.73e-01 (gamma:  2.73e-01, delta: 3.01e-02)
#>   Iteration 138: max param change = 2.60e-01 (gamma:  2.60e-01, delta: 2.87e-02)
#>   Iteration 139: max param change = 2.48e-01 (gamma:  2.48e-01, delta: 2.74e-02)
#>   Iteration 140: max param change = 2.36e-01 (gamma:  2.36e-01, delta: 2.61e-02)
#>   Iteration 141: max param change = 2.25e-01 (gamma:  2.25e-01, delta: 2.48e-02)
#>   Iteration 142: max param change = 2.14e-01 (gamma:  2.14e-01, delta: 2.36e-02)
#>   Iteration 143: max param change = 2.04e-01 (gamma:  2.04e-01, delta: 2.25e-02)
#>   Iteration 144: max param change = 1.94e-01 (gamma:  1.94e-01, delta: 2.15e-02)
#>   Iteration 145: max param change = 1.85e-01 (gamma:  1.85e-01, delta: 2.04e-02)
#>   Iteration 146: max param change = 1.76e-01 (gamma:  1.76e-01, delta: 1.95e-02)
#>   Iteration 147: max param change = 1.68e-01 (gamma:  1.68e-01, delta: 1.85e-02)
#>   Iteration 148: max param change = 1.60e-01 (gamma:  1.60e-01, delta: 1.77e-02)
#>   Iteration 149: max param change = 1.53e-01 (gamma:  1.53e-01, delta: 1.68e-02)
#>   Iteration 150: max param change = 1.45e-01 (gamma:  1.45e-01, delta: 1.60e-02)
#>   Iteration 151: max param change = 1.38e-01 (gamma:  1.38e-01, delta: 1.53e-02)
#>   Iteration 152: max param change = 1.32e-01 (gamma:  1.32e-01, delta: 1.45e-02)
#>   Iteration 153: max param change = 1.26e-01 (gamma:  1.26e-01, delta: 1.39e-02)
#>   Iteration 154: max param change = 1.20e-01 (gamma:  1.20e-01, delta: 1.32e-02)
#>   Iteration 155: max param change = 1.14e-01 (gamma:  1.14e-01, delta: 1.26e-02)
#>   Iteration 156: max param change = 1.09e-01 (gamma:  1.09e-01, delta: 1.20e-02)
#>   Iteration 157: max param change = 1.03e-01 (gamma:  1.03e-01, delta: 1.14e-02)
#>   Iteration 158: max param change = 9.85e-02 (gamma:  9.85e-02, delta: 1.09e-02)
#>   Iteration 159: max param change = 9.38e-02 (gamma:  9.38e-02, delta: 1.03e-02)
#>   Iteration 160: max param change = 8.94e-02 (gamma:  8.94e-02, delta: 9.86e-03)
#>   Iteration 161: max param change = 8.52e-02 (gamma:  8.52e-02, delta: 9.39e-03)
#>   Iteration 162: max param change = 8.11e-02 (gamma:  8.11e-02, delta: 8.95e-03)
#>   Iteration 163: max param change = 7.73e-02 (gamma:  7.73e-02, delta: 8.52e-03)
#>   Iteration 164: max param change = 7.36e-02 (gamma:  7.36e-02, delta: 8.12e-03)
#>   Iteration 165: max param change = 7.01e-02 (gamma:  7.01e-02, delta: 7.73e-03)
#>   Iteration 166: max param change = 6.68e-02 (gamma:  6.68e-02, delta: 7.36e-03)
#>   Iteration 167: max param change = 6.36e-02 (gamma:  6.36e-02, delta: 7.01e-03)
#>   Iteration 168: max param change = 6.06e-02 (gamma:  6.06e-02, delta: 6.68e-03)
#>   Iteration 169: max param change = 5.77e-02 (gamma:  5.77e-02, delta: 6.36e-03)
#>   Iteration 170: max param change = 5.50e-02 (gamma:  5.50e-02, delta: 6.06e-03)
#>   Iteration 171: max param change = 5.24e-02 (gamma:  5.24e-02, delta: 5.78e-03)
#>   Iteration 172: max param change = 4.99e-02 (gamma:  4.99e-02, delta: 5.50e-03)
#>   Iteration 173: max param change = 4.75e-02 (gamma:  4.75e-02, delta: 5.24e-03)
#>   Iteration 174: max param change = 4.53e-02 (gamma:  4.53e-02, delta: 4.99e-03)
#>   Iteration 175: max param change = 4.31e-02 (gamma:  4.31e-02, delta: 4.75e-03)
#>   Iteration 176: max param change = 4.11e-02 (gamma:  4.11e-02, delta: 4.53e-03)
#>   Iteration 177: max param change = 3.91e-02 (gamma:  3.91e-02, delta: 4.31e-03)
#>   Iteration 178: max param change = 3.73e-02 (gamma:  3.73e-02, delta: 4.11e-03)
#>   Iteration 179: max param change = 3.55e-02 (gamma:  3.55e-02, delta: 3.91e-03)
#>   Iteration 180: max param change = 3.38e-02 (gamma:  3.38e-02, delta: 3.73e-03)
#>   Iteration 181: max param change = 3.22e-02 (gamma:  3.22e-02, delta: 3.55e-03)
#>   Iteration 182: max param change = 3.07e-02 (gamma:  3.07e-02, delta: 3.38e-03)
#>   Iteration 183: max param change = 2.92e-02 (gamma:  2.92e-02, delta: 3.22e-03)
#>   Iteration 184: max param change = 2.78e-02 (gamma:  2.78e-02, delta: 3.07e-03)
#>   Iteration 185: max param change = 2.65e-02 (gamma:  2.65e-02, delta: 2.92e-03)
#>   Iteration 186: max param change = 2.53e-02 (gamma:  2.53e-02, delta: 2.78e-03)
#>   Iteration 187: max param change = 2.41e-02 (gamma:  2.41e-02, delta: 2.65e-03)
#>   Iteration 188: max param change = 2.29e-02 (gamma:  2.29e-02, delta: 2.53e-03)
#>   Iteration 189: max param change = 2.18e-02 (gamma:  2.18e-02, delta: 2.41e-03)
#>   Iteration 190: max param change = 2.08e-02 (gamma:  2.08e-02, delta: 2.29e-03)
#>   Iteration 191: max param change = 1.98e-02 (gamma:  1.98e-02, delta: 2.18e-03)
#>   Iteration 192: max param change = 1.89e-02 (gamma:  1.89e-02, delta: 2.08e-03)
#>   Iteration 193: max param change = 1.80e-02 (gamma:  1.80e-02, delta: 1.98e-03)
#>   Iteration 194: max param change = 1.71e-02 (gamma:  1.71e-02, delta: 1.89e-03)
#>   Iteration 195: max param change = 1.63e-02 (gamma:  1.63e-02, delta: 1.80e-03)
#>   Iteration 196: max param change = 1.55e-02 (gamma:  1.55e-02, delta: 1.71e-03)
#>   Iteration 197: max param change = 1.48e-02 (gamma:  1.48e-02, delta: 1.63e-03)
#>   Iteration 198: max param change = 1.41e-02 (gamma:  1.41e-02, delta: 1.55e-03)
#>   Iteration 199: max param change = 1.34e-02 (gamma:  1.34e-02, delta: 1.48e-03)
#>   Iteration 200: max param change = 1.28e-02 (gamma:  1.28e-02, delta: 1.41e-03)
#>   Iteration 201: max param change = 1.22e-02 (gamma:  1.22e-02, delta: 1.34e-03)
#>   Iteration 202: max param change = 1.16e-02 (gamma:  1.16e-02, delta: 1.28e-03)
#>   Iteration 203: max param change = 1.10e-02 (gamma:  1.10e-02, delta: 1.22e-03)
#>   Iteration 204: max param change = 1.05e-02 (gamma:  1.05e-02, delta: 1.16e-03)
#>   Iteration 205: max param change = 1.00e-02 (gamma:  1.00e-02, delta: 1.11e-03)
#>   Iteration 206: max param change = 9.55e-03 (gamma:  9.55e-03, delta: 1.05e-03)
#>   Iteration 207: max param change = 9.10e-03 (gamma:  9.10e-03, delta: 1.00e-03)
#>   Iteration 208: max param change = 8.66e-03 (gamma:  8.66e-03, delta: 9.55e-04)
#>   Iteration 209: max param change = 8.25e-03 (gamma:  8.25e-03, delta: 9.10e-04)
#>   Iteration 210: max param change = 7.86e-03 (gamma:  7.86e-03, delta: 8.67e-04)
#>   Iteration 211: max param change = 7.49e-03 (gamma:  7.49e-03, delta: 8.25e-04)
#>   Iteration 212: max param change = 7.13e-03 (gamma:  7.13e-03, delta: 7.86e-04)
#>   Iteration 213: max param change = 6.79e-03 (gamma:  6.79e-03, delta: 7.49e-04)
#>   Iteration 214: max param change = 6.47e-03 (gamma:  6.47e-03, delta: 7.13e-04)
#>   Iteration 215: max param change = 6.16e-03 (gamma:  6.16e-03, delta: 6.80e-04)
#>   Iteration 216: max param change = 5.87e-03 (gamma:  5.87e-03, delta: 6.47e-04)
#>   Iteration 217: max param change = 5.59e-03 (gamma:  5.59e-03, delta: 6.17e-04)
#>   Iteration 218: max param change = 5.33e-03 (gamma:  5.33e-03, delta: 5.87e-04)
#>   Iteration 219: max param change = 5.07e-03 (gamma:  5.07e-03, delta: 5.59e-04)
#>   Iteration 220: max param change = 4.83e-03 (gamma:  4.83e-03, delta: 5.33e-04)
#>   Iteration 221: max param change = 4.60e-03 (gamma:  4.60e-03, delta: 5.08e-04)
#>   Iteration 222: max param change = 4.39e-03 (gamma:  4.39e-03, delta: 4.83e-04)
#>   Iteration 223: max param change = 4.18e-03 (gamma:  4.18e-03, delta: 4.60e-04)
#>   Iteration 224: max param change = 3.98e-03 (gamma:  3.98e-03, delta: 4.39e-04)
#>   Iteration 225: max param change = 3.79e-03 (gamma:  3.79e-03, delta: 4.18e-04)
#>   Iteration 226: max param change = 3.61e-03 (gamma:  3.61e-03, delta: 3.98e-04)
#>   Iteration 227: max param change = 3.44e-03 (gamma:  3.44e-03, delta: 3.79e-04)
#>   Iteration 228: max param change = 3.28e-03 (gamma:  3.28e-03, delta: 3.61e-04)
#>   Iteration 229: max param change = 3.12e-03 (gamma:  3.12e-03, delta: 3.44e-04)
#>   Iteration 230: max param change = 2.97e-03 (gamma:  2.97e-03, delta: 3.28e-04)
#>   Iteration 231: max param change = 2.83e-03 (gamma:  2.83e-03, delta: 3.12e-04)
#>   Iteration 232: max param change = 2.70e-03 (gamma:  2.70e-03, delta: 2.97e-04)
#>   Iteration 233: max param change = 2.57e-03 (gamma:  2.57e-03, delta: 2.83e-04)
#>   Iteration 234: max param change = 2.45e-03 (gamma:  2.45e-03, delta: 2.70e-04)
#>   Iteration 235: max param change = 2.33e-03 (gamma:  2.33e-03, delta: 2.57e-04)
#>   Iteration 236: max param change = 2.22e-03 (gamma:  2.22e-03, delta: 2.45e-04)
#>   Iteration 237: max param change = 2.11e-03 (gamma:  2.11e-03, delta: 2.33e-04)
#>   Iteration 238: max param change = 2.01e-03 (gamma:  2.01e-03, delta: 2.22e-04)
#>   Iteration 239: max param change = 1.92e-03 (gamma:  1.92e-03, delta: 2.11e-04)
#>   Iteration 240: max param change = 1.83e-03 (gamma:  1.83e-03, delta: 2.01e-04)
#>   Iteration 241: max param change = 1.74e-03 (gamma:  1.74e-03, delta: 1.92e-04)
#>   Iteration 242: max param change = 1.66e-03 (gamma:  1.66e-03, delta: 1.83e-04)
#>   Iteration 243: max param change = 1.58e-03 (gamma:  1.58e-03, delta: 1.74e-04)
#>   Iteration 244: max param change = 1.50e-03 (gamma:  1.50e-03, delta: 1.66e-04)
#>   Iteration 245: max param change = 1.43e-03 (gamma:  1.43e-03, delta: 1.58e-04)
#>   Iteration 246: max param change = 1.36e-03 (gamma:  1.36e-03, delta: 1.50e-04)
#>   Iteration 247: max param change = 1.30e-03 (gamma:  1.30e-03, delta: 1.43e-04)
#>   Iteration 248: max param change = 1.24e-03 (gamma:  1.24e-03, delta: 1.36e-04)
#>   Iteration 249: max param change = 1.18e-03 (gamma:  1.18e-03, delta: 1.30e-04)
#>   Iteration 250: max param change = 1.12e-03 (gamma:  1.12e-03, delta: 1.24e-04)
#>   Iteration 251: max param change = 1.07e-03 (gamma:  1.07e-03, delta: 1.18e-04)
#>   Iteration 252: max param change = 1.02e-03 (gamma:  1.02e-03, delta: 1.12e-04)
#>   Iteration 253: max param change = 9.71e-04 (gamma:  9.71e-04, delta: 1.07e-04)
#>   Iteration 254: max param change = 9.25e-04 (gamma:  9.25e-04, delta: 1.02e-04)
#>   Iteration 255: max param change = 8.81e-04 (gamma:  8.81e-04, delta: 9.71e-05)
#>   Iteration 256: max param change = 8.39e-04 (gamma:  8.39e-04, delta: 9.25e-05)
#>   Iteration 257: max param change = 7.99e-04 (gamma:  7.99e-04, delta: 8.81e-05)
#>   Iteration 258: max param change = 7.61e-04 (gamma:  7.61e-04, delta: 8.39e-05)
#>   Iteration 259: max param change = 7.25e-04 (gamma:  7.25e-04, delta: 7.99e-05)
#>   Iteration 260: max param change = 6.91e-04 (gamma:  6.91e-04, delta: 7.61e-05)
#>   Iteration 261: max param change = 6.58e-04 (gamma:  6.58e-04, delta: 7.25e-05)
#>   Iteration 262: max param change = 6.27e-04 (gamma:  6.27e-04, delta: 6.91e-05)
#>   Iteration 263: max param change = 5.97e-04 (gamma:  5.97e-04, delta: 6.58e-05)
#>   Iteration 264: max param change = 5.69e-04 (gamma:  5.69e-04, delta: 6.27e-05)
#>   Iteration 265: max param change = 5.42e-04 (gamma:  5.42e-04, delta: 5.97e-05)
#>   Iteration 266: max param change = 5.16e-04 (gamma:  5.16e-04, delta: 5.69e-05)
#>   Iteration 267: max param change = 4.91e-04 (gamma:  4.91e-04, delta: 5.42e-05)
#>   Iteration 268: max param change = 4.68e-04 (gamma:  4.68e-04, delta: 5.16e-05)
#>   Iteration 269: max param change = 4.46e-04 (gamma:  4.46e-04, delta: 4.91e-05)
#>   Iteration 270: max param change = 4.25e-04 (gamma:  4.25e-04, delta: 4.68e-05)
#>   Iteration 271: max param change = 4.04e-04 (gamma:  4.04e-04, delta: 4.46e-05)
#>   Iteration 272: max param change = 3.85e-04 (gamma:  3.85e-04, delta: 4.25e-05)
#>   Iteration 273: max param change = 3.67e-04 (gamma:  3.67e-04, delta: 4.05e-05)
#>   Iteration 274: max param change = 3.50e-04 (gamma:  3.50e-04, delta: 3.85e-05)
#>   Iteration 275: max param change = 3.33e-04 (gamma:  3.33e-04, delta: 3.67e-05)
#>   Iteration 276: max param change = 3.17e-04 (gamma:  3.17e-04, delta: 3.50e-05)
#>   Iteration 277: max param change = 3.02e-04 (gamma:  3.02e-04, delta: 3.33e-05)
#>   Iteration 278: max param change = 2.88e-04 (gamma:  2.88e-04, delta: 3.17e-05)
#>   Iteration 279: max param change = 2.74e-04 (gamma:  2.74e-04, delta: 3.02e-05)
#>   Iteration 280: max param change = 2.61e-04 (gamma:  2.61e-04, delta: 2.88e-05)
#>   Iteration 281: max param change = 2.49e-04 (gamma:  2.49e-04, delta: 2.74e-05)
#>   Iteration 282: max param change = 2.37e-04 (gamma:  2.37e-04, delta: 2.61e-05)
#>   Iteration 283: max param change = 2.26e-04 (gamma:  2.26e-04, delta: 2.49e-05)
#>   Iteration 284: max param change = 2.15e-04 (gamma:  2.15e-04, delta: 2.37e-05)
#>   Iteration 285: max param change = 2.05e-04 (gamma:  2.05e-04, delta: 2.26e-05)
#>   Iteration 286: max param change = 1.95e-04 (gamma:  1.95e-04, delta: 2.15e-05)
#>   Iteration 287: max param change = 1.86e-04 (gamma:  1.86e-04, delta: 2.05e-05)
#>   Iteration 288: max param change = 1.77e-04 (gamma:  1.77e-04, delta: 1.95e-05)
#>   Iteration 289: max param change = 1.69e-04 (gamma:  1.69e-04, delta: 1.86e-05)
#>   Iteration 290: max param change = 1.61e-04 (gamma:  1.61e-04, delta: 1.77e-05)
#>   Iteration 291: max param change = 1.53e-04 (gamma:  1.53e-04, delta: 1.69e-05)
#>   Iteration 292: max param change = 1.46e-04 (gamma:  1.46e-04, delta: 1.61e-05)
#>   Iteration 293: max param change = 1.39e-04 (gamma:  1.39e-04, delta: 1.53e-05)
#>   Iteration 294: max param change = 1.32e-04 (gamma:  1.32e-04, delta: 1.46e-05)
#>   Iteration 295: max param change = 1.26e-04 (gamma:  1.26e-04, delta: 1.39e-05)
#>   Iteration 296: max param change = 1.20e-04 (gamma:  1.20e-04, delta: 1.32e-05)
#>   Iteration 297: max param change = 1.14e-04 (gamma:  1.14e-04, delta: 1.26e-05)
#>   Iteration 298: max param change = 1.09e-04 (gamma:  1.09e-04, delta: 1.20e-05)
#>   Iteration 299: max param change = 1.04e-04 (gamma:  1.04e-04, delta: 1.14e-05)
#>   Iteration 300: max param change = 9.87e-05 (gamma:  9.87e-05, delta: 1.09e-05)
#>   Iteration 301: max param change = 9.40e-05 (gamma:  9.40e-05, delta: 1.04e-05)
#>   Iteration 302: max param change = 8.95e-05 (gamma:  8.95e-05, delta: 9.87e-06)
#>   Iteration 303: max param change = 8.53e-05 (gamma:  8.53e-05, delta: 9.40e-06)
#>   Iteration 304: max param change = 8.12e-05 (gamma:  8.12e-05, delta: 8.96e-06)
#>   Iteration 305: max param change = 7.74e-05 (gamma:  7.74e-05, delta: 8.53e-06)
#>   Iteration 306: max param change = 7.37e-05 (gamma:  7.37e-05, delta: 8.13e-06)
#>   Iteration 307: max param change = 7.02e-05 (gamma:  7.02e-05, delta: 7.74e-06)
#>   Iteration 308: max param change = 6.69e-05 (gamma:  6.69e-05, delta: 7.37e-06)
#>   Iteration 309: max param change = 6.37e-05 (gamma:  6.37e-05, delta: 7.02e-06)
#>   Iteration 310: max param change = 6.07e-05 (gamma:  6.07e-05, delta: 6.69e-06)
#>   Iteration 311: max param change = 5.78e-05 (gamma:  5.78e-05, delta: 6.37e-06)
#>   Iteration 312: max param change = 5.51e-05 (gamma:  5.51e-05, delta: 6.07e-06)
#>   Iteration 313: max param change = 5.24e-05 (gamma:  5.24e-05, delta: 5.78e-06)
#>   Iteration 314: max param change = 4.99e-05 (gamma:  4.99e-05, delta: 5.51e-06)
#>   Iteration 315: max param change = 4.76e-05 (gamma:  4.76e-05, delta: 5.25e-06)
#>   Iteration 316: max param change = 4.53e-05 (gamma:  4.53e-05, delta: 5.00e-06)
#>   Iteration 317: max param change = 4.32e-05 (gamma:  4.32e-05, delta: 4.76e-06)
#>   Iteration 318: max param change = 4.11e-05 (gamma:  4.11e-05, delta: 4.53e-06)
#>   Iteration 319: max param change = 3.92e-05 (gamma:  3.92e-05, delta: 4.32e-06)
#>   Iteration 320: max param change = 3.73e-05 (gamma:  3.73e-05, delta: 4.11e-06)
#>   Iteration 321: max param change = 3.55e-05 (gamma:  3.55e-05, delta: 3.92e-06)
#>   Iteration 322: max param change = 3.38e-05 (gamma:  3.38e-05, delta: 3.73e-06)
#>   Iteration 323: max param change = 3.22e-05 (gamma:  3.22e-05, delta: 3.55e-06)
#>   Iteration 324: max param change = 3.07e-05 (gamma:  3.07e-05, delta: 3.39e-06)
#>   Iteration 325: max param change = 2.93e-05 (gamma:  2.93e-05, delta: 3.22e-06)
#>   Iteration 326: max param change = 2.79e-05 (gamma:  2.79e-05, delta: 3.07e-06)
#>   Iteration 327: max param change = 2.65e-05 (gamma:  2.65e-05, delta: 2.93e-06)
#>   Iteration 328: max param change = 2.53e-05 (gamma:  2.53e-05, delta: 2.79e-06)
#>   Iteration 329: max param change = 2.41e-05 (gamma:  2.41e-05, delta: 2.65e-06)
#>   Iteration 330: max param change = 2.29e-05 (gamma:  2.29e-05, delta: 2.53e-06)
#>   Iteration 331: max param change = 2.18e-05 (gamma:  2.18e-05, delta: 2.41e-06)
#>   Iteration 332: max param change = 2.08e-05 (gamma:  2.08e-05, delta: 2.29e-06)
#>   Iteration 333: max param change = 1.98e-05 (gamma:  1.98e-05, delta: 2.19e-06)
#>   Iteration 334: max param change = 1.89e-05 (gamma:  1.89e-05, delta: 2.08e-06)
#>   Iteration 335: max param change = 1.80e-05 (gamma:  1.80e-05, delta: 1.98e-06)
#>   Iteration 336: max param change = 1.71e-05 (gamma:  1.71e-05, delta: 1.89e-06)
#>   Iteration 337: max param change = 1.63e-05 (gamma:  1.63e-05, delta: 1.80e-06)
#>   Iteration 338: max param change = 1.55e-05 (gamma:  1.55e-05, delta: 1.71e-06)
#>   Iteration 339: max param change = 1.48e-05 (gamma:  1.48e-05, delta: 1.63e-06)
#>   Iteration 340: max param change = 1.41e-05 (gamma:  1.41e-05, delta: 1.55e-06)
#>   Iteration 341: max param change = 1.34e-05 (gamma:  1.34e-05, delta: 1.48e-06)
#>   Iteration 342: max param change = 1.28e-05 (gamma:  1.28e-05, delta: 1.41e-06)
#>   Iteration 343: max param change = 1.22e-05 (gamma:  1.22e-05, delta: 1.34e-06)
#>   Iteration 344: max param change = 1.16e-05 (gamma:  1.16e-05, delta: 1.28e-06)
#>   Iteration 345: max param change = 1.11e-05 (gamma:  1.11e-05, delta: 1.22e-06)
#>   Iteration 346: max param change = 1.05e-05 (gamma:  1.05e-05, delta: 1.16e-06)
#>   Iteration 347: max param change = 1.00e-05 (gamma:  1.00e-05, delta: 1.11e-06)
#>   Iteration 348: max param change = 9.56e-06 (gamma:  9.56e-06, delta: 1.05e-06)
#>   Iteration 349: max param change = 9.10e-06 (gamma:  9.10e-06, delta: 1.00e-06)
#>   Iteration 350: max param change = 8.67e-06 (gamma:  8.67e-06, delta: 9.56e-07)
#>   Iteration 351: max param change = 8.26e-06 (gamma:  8.26e-06, delta: 9.10e-07)
#>   Iteration 352: max param change = 7.87e-06 (gamma:  7.87e-06, delta: 8.67e-07)
#>   Iteration 353: max param change = 7.49e-06 (gamma:  7.49e-06, delta: 8.26e-07)
#>   Iteration 354: max param change = 7.14e-06 (gamma:  7.14e-06, delta: 7.87e-07)
#>   Iteration 355: max param change = 6.80e-06 (gamma:  6.80e-06, delta: 7.49e-07)
#>   Iteration 356: max param change = 6.48e-06 (gamma:  6.48e-06, delta: 7.14e-07)
#>   Iteration 357: max param change = 6.17e-06 (gamma:  6.17e-06, delta: 6.80e-07)
#>   Iteration 358: max param change = 5.88e-06 (gamma:  5.88e-06, delta: 6.48e-07)
#>   Iteration 359: max param change = 5.60e-06 (gamma:  5.60e-06, delta: 6.17e-07)
#>   Iteration 360: max param change = 5.33e-06 (gamma:  5.33e-06, delta: 5.88e-07)
#>   Iteration 361: max param change = 5.08e-06 (gamma:  5.08e-06, delta: 5.60e-07)
#>   Iteration 362: max param change = 4.84e-06 (gamma:  4.84e-06, delta: 5.33e-07)
#>   Iteration 363: max param change = 4.61e-06 (gamma:  4.61e-06, delta: 5.08e-07)
#>   Iteration 364: max param change = 4.39e-06 (gamma:  4.39e-06, delta: 4.84e-07)
#>   Iteration 365: max param change = 4.18e-06 (gamma:  4.18e-06, delta: 4.61e-07)
#>   Iteration 366: max param change = 3.98e-06 (gamma:  3.98e-06, delta: 4.39e-07)
#>   Iteration 367: max param change = 3.79e-06 (gamma:  3.79e-06, delta: 4.18e-07)
#>   Iteration 368: max param change = 3.61e-06 (gamma:  3.61e-06, delta: 3.98e-07)
#>   Iteration 369: max param change = 3.44e-06 (gamma:  3.44e-06, delta: 3.79e-07)
#>   Iteration 370: max param change = 3.28e-06 (gamma:  3.28e-06, delta: 3.61e-07)
#>   Iteration 371: max param change = 3.12e-06 (gamma:  3.12e-06, delta: 3.44e-07)
#>   Iteration 372: max param change = 2.97e-06 (gamma:  2.97e-06, delta: 3.28e-07)
#>   Iteration 373: max param change = 2.83e-06 (gamma:  2.83e-06, delta: 3.12e-07)
#>   Iteration 374: max param change = 2.70e-06 (gamma:  2.70e-06, delta: 2.97e-07)
#>   Iteration 375: max param change = 2.57e-06 (gamma:  2.57e-06, delta: 2.83e-07)
#>   Iteration 376: max param change = 2.45e-06 (gamma:  2.45e-06, delta: 2.70e-07)
#>   Iteration 377: max param change = 2.33e-06 (gamma:  2.33e-06, delta: 2.57e-07)
#>   Iteration 378: max param change = 2.22e-06 (gamma:  2.22e-06, delta: 2.45e-07)
#>   Iteration 379: max param change = 2.12e-06 (gamma:  2.12e-06, delta: 2.33e-07)
#>   Iteration 380: max param change = 2.02e-06 (gamma:  2.02e-06, delta: 2.22e-07)
#>   Iteration 381: max param change = 1.92e-06 (gamma:  1.92e-06, delta: 2.12e-07)
#>   Iteration 382: max param change = 1.83e-06 (gamma:  1.83e-06, delta: 2.02e-07)
#>   Iteration 383: max param change = 1.74e-06 (gamma:  1.74e-06, delta: 1.92e-07)
#>   Iteration 384: max param change = 1.66e-06 (gamma:  1.66e-06, delta: 1.83e-07)
#>   Iteration 385: max param change = 1.58e-06 (gamma:  1.58e-06, delta: 1.74e-07)
#>   Iteration 386: max param change = 1.51e-06 (gamma:  1.51e-06, delta: 1.66e-07)
#>   Iteration 387: max param change = 1.43e-06 (gamma:  1.43e-06, delta: 1.58e-07)
#>   Iteration 388: max param change = 1.37e-06 (gamma:  1.37e-06, delta: 1.51e-07)
#>   Iteration 389: max param change = 1.30e-06 (gamma:  1.30e-06, delta: 1.43e-07)
#>   Iteration 390: max param change = 1.24e-06 (gamma:  1.24e-06, delta: 1.37e-07)
#>   Iteration 391: max param change = 1.18e-06 (gamma:  1.18e-06, delta: 1.30e-07)
#>   Iteration 392: max param change = 1.12e-06 (gamma:  1.12e-06, delta: 1.24e-07)
#>   Iteration 393: max param change = 1.07e-06 (gamma:  1.07e-06, delta: 1.18e-07)
#>   Iteration 394: max param change = 1.02e-06 (gamma:  1.02e-06, delta: 1.12e-07)
#>   Iteration 395: max param change = 9.71e-07 (gamma:  9.71e-07, delta: 1.07e-07)
#>   Iteration 396: max param change = 9.25e-07 (gamma:  9.25e-07, delta: 1.02e-07)
#>   Iteration 397: max param change = 8.81e-07 (gamma:  8.81e-07, delta: 9.72e-08)
#>   Iteration 398: max param change = 8.40e-07 (gamma:  8.40e-07, delta: 9.26e-08)
#>   Iteration 399: max param change = 8.00e-07 (gamma:  8.00e-07, delta: 8.82e-08)
#>   Iteration 400: max param change = 7.62e-07 (gamma:  7.62e-07, delta: 8.40e-08)
#>   Iteration 401: max param change = 7.26e-07 (gamma:  7.26e-07, delta: 8.00e-08)
#>   Iteration 402: max param change = 6.91e-07 (gamma:  6.91e-07, delta: 7.62e-08)
#>   Iteration 403: max param change = 6.58e-07 (gamma:  6.58e-07, delta: 7.26e-08)
#>   Iteration 404: max param change = 6.27e-07 (gamma:  6.27e-07, delta: 6.91e-08)
#>   Iteration 405: max param change = 5.97e-07 (gamma:  5.97e-07, delta: 6.58e-08)
#>   Iteration 406: max param change = 5.69e-07 (gamma:  5.69e-07, delta: 6.27e-08)
#>   Iteration 407: max param change = 5.42e-07 (gamma:  5.42e-07, delta: 5.97e-08)
#>   Iteration 408: max param change = 5.16e-07 (gamma:  5.16e-07, delta: 5.69e-08)
#>   Iteration 409: max param change = 4.92e-07 (gamma:  4.92e-07, delta: 5.42e-08)
#>   Iteration 410: max param change = 4.68e-07 (gamma:  4.68e-07, delta: 5.16e-08)
#>   Iteration 411: max param change = 4.46e-07 (gamma:  4.46e-07, delta: 4.92e-08)
#>   Iteration 412: max param change = 4.25e-07 (gamma:  4.25e-07, delta: 4.68e-08)
#>   Iteration 413: max param change = 4.05e-07 (gamma:  4.05e-07, delta: 4.46e-08)
#>   Iteration 414: max param change = 3.86e-07 (gamma:  3.86e-07, delta: 4.25e-08)
#>   Iteration 415: max param change = 3.67e-07 (gamma:  3.67e-07, delta: 4.05e-08)
#>   Iteration 416: max param change = 3.50e-07 (gamma:  3.50e-07, delta: 3.86e-08)
#>   Iteration 417: max param change = 3.33e-07 (gamma:  3.33e-07, delta: 3.67e-08)
#>   Iteration 418: max param change = 3.17e-07 (gamma:  3.17e-07, delta: 3.50e-08)
#>   Iteration 419: max param change = 3.02e-07 (gamma:  3.02e-07, delta: 3.33e-08)
#>   Iteration 420: max param change = 2.88e-07 (gamma:  2.88e-07, delta: 3.17e-08)
#>   Iteration 421: max param change = 2.74e-07 (gamma:  2.74e-07, delta: 3.02e-08)
#>   Iteration 422: max param change = 2.61e-07 (gamma:  2.61e-07, delta: 2.88e-08)
#>   Iteration 423: max param change = 2.49e-07 (gamma:  2.49e-07, delta: 2.74e-08)
#>   Iteration 424: max param change = 2.37e-07 (gamma:  2.37e-07, delta: 2.61e-08)
#>   Iteration 425: max param change = 2.26e-07 (gamma:  2.26e-07, delta: 2.49e-08)
#>   Iteration 426: max param change = 2.15e-07 (gamma:  2.15e-07, delta: 2.37e-08)
#>   Iteration 427: max param change = 2.05e-07 (gamma:  2.05e-07, delta: 2.26e-08)
#>   Iteration 428: max param change = 1.95e-07 (gamma:  1.95e-07, delta: 2.15e-08)
#>   Iteration 429: max param change = 1.86e-07 (gamma:  1.86e-07, delta: 2.05e-08)
#>   Iteration 430: max param change = 1.77e-07 (gamma:  1.77e-07, delta: 1.95e-08)
#>   Iteration 431: max param change = 1.69e-07 (gamma:  1.69e-07, delta: 1.86e-08)
#>   Iteration 432: max param change = 1.61e-07 (gamma:  1.61e-07, delta: 1.77e-08)
#>   Iteration 433: max param change = 1.53e-07 (gamma:  1.53e-07, delta: 1.69e-08)
#>   Iteration 434: max param change = 1.46e-07 (gamma:  1.46e-07, delta: 1.61e-08)
#>   Iteration 435: max param change = 1.39e-07 (gamma:  1.39e-07, delta: 1.53e-08)
#>   Iteration 436: max param change = 1.32e-07 (gamma:  1.32e-07, delta: 1.46e-08)
#>   Iteration 437: max param change = 1.26e-07 (gamma:  1.26e-07, delta: 1.39e-08)
#>   Iteration 438: max param change = 1.20e-07 (gamma:  1.20e-07, delta: 1.32e-08)
#>   Iteration 439: max param change = 1.14e-07 (gamma:  1.14e-07, delta: 1.26e-08)
#>   Iteration 440: max param change = 1.09e-07 (gamma:  1.09e-07, delta: 1.20e-08)
#>   Iteration 441: max param change = 1.04e-07 (gamma:  1.04e-07, delta: 1.14e-08)
#>   Iteration 442: max param change = 9.88e-08 (gamma:  9.88e-08, delta: 1.09e-08)
#>   Iteration 443: max param change = 9.41e-08 (gamma:  9.41e-08, delta: 1.04e-08)
#>   Iteration 444: max param change = 8.96e-08 (gamma:  8.96e-08, delta: 9.88e-09)
#>   Iteration 445: max param change = 8.53e-08 (gamma:  8.53e-08, delta: 9.41e-09)
#>   Iteration 446: max param change = 8.13e-08 (gamma:  8.13e-08, delta: 8.96e-09)
#>   Iteration 447: max param change = 7.74e-08 (gamma:  7.74e-08, delta: 8.54e-09)
#>   Iteration 448: max param change = 7.38e-08 (gamma:  7.38e-08, delta: 8.13e-09)
#>   Iteration 449: max param change = 7.03e-08 (gamma:  7.03e-08, delta: 7.75e-09)
#>   Iteration 450: max param change = 6.69e-08 (gamma:  6.69e-08, delta: 7.38e-09)
#>   Iteration 451: max param change = 6.37e-08 (gamma:  6.37e-08, delta: 7.03e-09)
#>   Iteration 452: max param change = 6.07e-08 (gamma:  6.07e-08, delta: 6.69e-09)
#>   Iteration 453: max param change = 5.78e-08 (gamma:  5.78e-08, delta: 6.38e-09)
#>   Iteration 454: max param change = 5.51e-08 (gamma:  5.51e-08, delta: 6.07e-09)
#>   Iteration 455: max param change = 5.25e-08 (gamma:  5.25e-08, delta: 5.78e-09)
#>   Iteration 456: max param change = 5.00e-08 (gamma:  5.00e-08, delta: 5.51e-09)
#>   Iteration 457: max param change = 4.76e-08 (gamma:  4.76e-08, delta: 5.25e-09)
#>   Iteration 458: max param change = 4.54e-08 (gamma:  4.54e-08, delta: 5.00e-09)
#>   Iteration 459: max param change = 4.32e-08 (gamma:  4.32e-08, delta: 4.76e-09)
#>   Iteration 460: max param change = 4.11e-08 (gamma:  4.11e-08, delta: 4.54e-09)
#>   Iteration 461: max param change = 3.92e-08 (gamma:  3.92e-08, delta: 4.32e-09)
#>   Iteration 462: max param change = 3.73e-08 (gamma:  3.73e-08, delta: 4.12e-09)
#>   Iteration 463: max param change = 3.56e-08 (gamma:  3.56e-08, delta: 3.92e-09)
#>   Iteration 464: max param change = 3.39e-08 (gamma:  3.39e-08, delta: 3.73e-09)
#>   Iteration 465: max param change = 3.23e-08 (gamma:  3.23e-08, delta: 3.56e-09)
#>   Iteration 466: max param change = 3.07e-08 (gamma:  3.07e-08, delta: 3.39e-09)
#>   Iteration 467: max param change = 2.93e-08 (gamma:  2.93e-08, delta: 3.23e-09)
#>   Iteration 468: max param change = 2.79e-08 (gamma:  2.79e-08, delta: 3.07e-09)
#>   Iteration 469: max param change = 2.66e-08 (gamma:  2.66e-08, delta: 2.93e-09)
#>   Iteration 470: max param change = 2.53e-08 (gamma:  2.53e-08, delta: 2.79e-09)
#>   Iteration 471: max param change = 2.41e-08 (gamma:  2.41e-08, delta: 2.66e-09)
#>   Iteration 472: max param change = 2.30e-08 (gamma:  2.30e-08, delta: 2.53e-09)
#>   Iteration 473: max param change = 2.19e-08 (gamma:  2.19e-08, delta: 2.41e-09)
#>   Iteration 474: max param change = 2.08e-08 (gamma:  2.08e-08, delta: 2.30e-09)
#>   Iteration 475: max param change = 1.98e-08 (gamma:  1.98e-08, delta: 2.19e-09)
#>   Iteration 476: max param change = 1.89e-08 (gamma:  1.89e-08, delta: 2.08e-09)
#>   Iteration 477: max param change = 1.80e-08 (gamma:  1.80e-08, delta: 1.98e-09)
#>   Iteration 478: max param change = 1.71e-08 (gamma:  1.71e-08, delta: 1.89e-09)
#>   Iteration 479: max param change = 1.63e-08 (gamma:  1.63e-08, delta: 1.80e-09)
#>   Iteration 480: max param change = 1.56e-08 (gamma:  1.56e-08, delta: 1.71e-09)
#>   Iteration 481: max param change = 1.48e-08 (gamma:  1.48e-08, delta: 1.63e-09)
#>   Iteration 482: max param change = 1.41e-08 (gamma:  1.41e-08, delta: 1.56e-09)
#>   Iteration 483: max param change = 1.34e-08 (gamma:  1.34e-08, delta: 1.48e-09)
#>   Iteration 484: max param change = 1.28e-08 (gamma:  1.28e-08, delta: 1.41e-09)
#>   Iteration 485: max param change = 1.22e-08 (gamma:  1.22e-08, delta: 1.34e-09)
#>   Iteration 486: max param change = 1.16e-08 (gamma:  1.16e-08, delta: 1.28e-09)
#>   Iteration 487: max param change = 1.11e-08 (gamma:  1.11e-08, delta: 1.22e-09)
#>   Iteration 488: max param change = 1.05e-08 (gamma:  1.05e-08, delta: 1.16e-09)
#>   Iteration 489: max param change = 1.00e-08 (gamma:  1.00e-08, delta: 1.11e-09)
#>   Iteration 490: max param change = 9.56e-09 (gamma:  9.56e-09, delta: 1.05e-09)
#>   Iteration 491: max param change = 9.11e-09 (gamma:  9.11e-09, delta: 1.00e-09)
#>   Iteration 492: max param change = 8.68e-09 (gamma:  8.68e-09, delta: 9.57e-10)
#>   Iteration 493: max param change = 8.27e-09 (gamma:  8.27e-09, delta: 9.11e-10)
#>   Iteration 494: max param change = 7.87e-09 (gamma:  7.87e-09, delta: 8.68e-10)
#>   Iteration 495: max param change = 7.50e-09 (gamma:  7.50e-09, delta: 8.27e-10)
#>   Iteration 496: max param change = 7.14e-09 (gamma:  7.14e-09, delta: 7.87e-10)
#>   Iteration 497: max param change = 6.80e-09 (gamma:  6.80e-09, delta: 7.50e-10)
#>   Iteration 498: max param change = 6.48e-09 (gamma:  6.48e-09, delta: 7.14e-10)
#>   Iteration 499: max param change = 6.17e-09 (gamma:  6.17e-09, delta: 6.80e-10)
#>   Iteration 500: max param change = 5.88e-09 (gamma:  5.88e-09, delta: 6.48e-10)
#>   Iteration 501: max param change = 5.60e-09 (gamma:  5.60e-09, delta: 6.18e-10)
#>   Iteration 502: max param change = 5.33e-09 (gamma:  5.33e-09, delta: 5.88e-10)
#>   Iteration 503: max param change = 5.08e-09 (gamma:  5.08e-09, delta: 5.60e-10)
#>   Iteration 504: max param change = 4.84e-09 (gamma:  4.84e-09, delta: 5.34e-10)
#>   Iteration 505: max param change = 4.61e-09 (gamma:  4.61e-09, delta: 5.08e-10)
#>   Iteration 506: max param change = 4.39e-09 (gamma:  4.39e-09, delta: 4.84e-10)
#>   Iteration 507: max param change = 4.18e-09 (gamma:  4.18e-09, delta: 4.61e-10)
#>   Iteration 508: max param change = 3.99e-09 (gamma:  3.99e-09, delta: 4.39e-10)
#>   Iteration 509: max param change = 3.80e-09 (gamma:  3.80e-09, delta: 4.19e-10)
#>   Iteration 510: max param change = 3.61e-09 (gamma:  3.61e-09, delta: 3.98e-10)
#>   Iteration 511: max param change = 3.44e-09 (gamma:  3.44e-09, delta: 3.80e-10)
#>   Iteration 512: max param change = 3.28e-09 (gamma:  3.28e-09, delta: 3.62e-10)
#>   Iteration 513: max param change = 3.13e-09 (gamma:  3.13e-09, delta: 3.44e-10)
#>   Iteration 514: max param change = 2.98e-09 (gamma:  2.98e-09, delta: 3.28e-10)
#>   Iteration 515: max param change = 2.83e-09 (gamma:  2.83e-09, delta: 3.12e-10)
#>   Iteration 516: max param change = 2.70e-09 (gamma:  2.70e-09, delta: 2.97e-10)
#>   Iteration 517: max param change = 2.57e-09 (gamma:  2.57e-09, delta: 2.84e-10)
#>   Iteration 518: max param change = 2.45e-09 (gamma:  2.45e-09, delta: 2.70e-10)
#>   Iteration 519: max param change = 2.33e-09 (gamma:  2.33e-09, delta: 2.57e-10)
#>   Iteration 520: max param change = 2.22e-09 (gamma:  2.22e-09, delta: 2.45e-10)
#>   Iteration 521: max param change = 2.12e-09 (gamma:  2.12e-09, delta: 2.34e-10)
#>   Iteration 522: max param change = 2.02e-09 (gamma:  2.02e-09, delta: 2.22e-10)
#>   Iteration 523: max param change = 1.92e-09 (gamma:  1.92e-09, delta: 2.12e-10)
#>   Iteration 524: max param change = 1.83e-09 (gamma:  1.83e-09, delta: 2.02e-10)
#>   Iteration 525: max param change = 1.74e-09 (gamma:  1.74e-09, delta: 1.92e-10)
#>   Iteration 526: max param change = 1.66e-09 (gamma:  1.66e-09, delta: 1.83e-10)
#>   Iteration 527: max param change = 1.58e-09 (gamma:  1.58e-09, delta: 1.74e-10)
#>   Iteration 528: max param change = 1.51e-09 (gamma:  1.51e-09, delta: 1.66e-10)
#>   Iteration 529: max param change = 1.43e-09 (gamma:  1.43e-09, delta: 1.58e-10)
#>   Iteration 530: max param change = 1.37e-09 (gamma:  1.37e-09, delta: 1.51e-10)
#>   Iteration 531: max param change = 1.30e-09 (gamma:  1.30e-09, delta: 1.43e-10)
#>   Iteration 532: max param change = 1.24e-09 (gamma:  1.24e-09, delta: 1.37e-10)
#>   Iteration 533: max param change = 1.18e-09 (gamma:  1.18e-09, delta: 1.30e-10)
#>   Iteration 534: max param change = 1.13e-09 (gamma:  1.13e-09, delta: 1.24e-10)
#>   Iteration 535: max param change = 1.07e-09 (gamma:  1.07e-09, delta: 1.18e-10)
#>   Iteration 536: max param change = 1.02e-09 (gamma:  1.02e-09, delta: 1.13e-10)
#>   Iteration 537: max param change = 9.72e-10 (gamma:  9.72e-10, delta: 1.07e-10)
#>   Iteration 538: max param change = 9.26e-10 (gamma:  9.26e-10, delta: 1.02e-10)
#>   Iteration 539: max param change = 8.81e-10 (gamma:  8.81e-10, delta: 9.70e-11)
#>   Iteration 540: max param change = 8.41e-10 (gamma:  8.41e-10, delta: 9.26e-11)
#>   Iteration 541: max param change = 8.00e-10 (gamma:  8.00e-10, delta: 8.82e-11)
#>   Iteration 542: max param change = 7.63e-10 (gamma:  7.63e-10, delta: 8.41e-11)
#>   Iteration 543: max param change = 7.26e-10 (gamma:  7.26e-10, delta: 7.99e-11)
#>   Iteration 544: max param change = 6.92e-10 (gamma:  6.92e-10, delta: 7.62e-11)
#>   Iteration 545: max param change = 6.59e-10 (gamma:  6.59e-10, delta: 7.26e-11)
#>   Iteration 546: max param change = 6.27e-10 (gamma:  6.27e-10, delta: 6.91e-11)
#>   Iteration 547: max param change = 5.97e-10 (gamma:  5.97e-10, delta: 6.60e-11)
#>   Iteration 548: max param change = 5.69e-10 (gamma:  5.69e-10, delta: 6.28e-11)
#>   Iteration 549: max param change = 5.42e-10 (gamma:  5.42e-10, delta: 5.97e-11)
#>   Iteration 550: max param change = 5.17e-10 (gamma:  5.17e-10, delta: 5.70e-11)
#>   Iteration 551: max param change = 4.92e-10 (gamma:  4.92e-10, delta: 5.41e-11)
#>   Iteration 552: max param change = 4.70e-10 (gamma:  4.70e-10, delta: 5.18e-11)
#>   Iteration 553: max param change = 4.46e-10 (gamma:  4.46e-10, delta: 4.91e-11)
#>   Iteration 554: max param change = 4.25e-10 (gamma:  4.25e-10, delta: 4.66e-11)
#>   Iteration 555: max param change = 4.05e-10 (gamma:  4.05e-10, delta: 4.48e-11)
#>   Iteration 556: max param change = 3.85e-10 (gamma:  3.85e-10, delta: 4.25e-11)
#>   Iteration 557: max param change = 3.68e-10 (gamma:  3.68e-10, delta: 4.06e-11)
#>   Iteration 558: max param change = 3.50e-10 (gamma:  3.50e-10, delta: 3.85e-11)
#>   Iteration 559: max param change = 3.33e-10 (gamma:  3.33e-10, delta: 3.66e-11)
#>   Iteration 560: max param change = 3.18e-10 (gamma:  3.18e-10, delta: 3.52e-11)
#>   Iteration 561: max param change = 3.02e-10 (gamma:  3.02e-10, delta: 3.33e-11)
#>   Iteration 562: max param change = 2.88e-10 (gamma:  2.88e-10, delta: 3.18e-11)
#>   Iteration 563: max param change = 2.75e-10 (gamma:  2.75e-10, delta: 3.04e-11)
#>   Iteration 564: max param change = 2.62e-10 (gamma:  2.62e-10, delta: 2.87e-11)
#>   Iteration 565: max param change = 2.49e-10 (gamma:  2.49e-10, delta: 2.75e-11)
#>   Iteration 566: max param change = 2.37e-10 (gamma:  2.37e-10, delta: 2.60e-11)
#>   Iteration 567: max param change = 2.26e-10 (gamma:  2.26e-10, delta: 2.52e-11)
#>   Iteration 568: max param change = 2.16e-10 (gamma:  2.16e-10, delta: 2.37e-11)
#>   Iteration 569: max param change = 2.05e-10 (gamma:  2.05e-10, delta: 2.27e-11)
#>   Iteration 570: max param change = 1.96e-10 (gamma:  1.96e-10, delta: 2.14e-11)
#>   Iteration 571: max param change = 1.86e-10 (gamma:  1.86e-10, delta: 2.04e-11)
#>   Iteration 572: max param change = 1.76e-10 (gamma:  1.76e-10, delta: 1.96e-11)
#>   Iteration 573: max param change = 1.70e-10 (gamma:  1.70e-10, delta: 1.87e-11)
#>   Iteration 574: max param change = 1.61e-10 (gamma:  1.61e-10, delta: 1.77e-11)
#>   Iteration 575: max param change = 1.54e-10 (gamma:  1.54e-10, delta: 1.69e-11)
#>   Iteration 576: max param change = 1.45e-10 (gamma:  1.45e-10, delta: 1.60e-11)
#>   Iteration 577: max param change = 1.39e-10 (gamma:  1.39e-10, delta: 1.52e-11)
#>   Iteration 578: max param change = 1.32e-10 (gamma:  1.32e-10, delta: 1.46e-11)
#>   Iteration 579: max param change = 1.26e-10 (gamma:  1.26e-10, delta: 1.39e-11)
#>   Iteration 580: max param change = 1.21e-10 (gamma:  1.21e-10, delta: 1.33e-11)
#>   Iteration 581: max param change = 1.14e-10 (gamma:  1.14e-10, delta: 1.27e-11)
#>   Iteration 582: max param change = 1.09e-10 (gamma:  1.09e-10, delta: 1.21e-11)
#>   Iteration 583: max param change = 1.04e-10 (gamma:  1.04e-10, delta: 1.15e-11)
#>   Iteration 584: max param change = 9.87e-11 (gamma:  9.87e-11, delta: 1.08e-11)
#> 
#> Converged after 584 iterations
print(fit)
#> Leunbach Score Parameter Estimation
#> ====================================
#> Power Series Distribution with Generalized Symmetric Functions
#> 
#> N = 500 observations
#> 
#> Test 1 score parameters (gamma):
#> 
#>  Score Frequency      Gamma Log_Gamma
#>      0        32    1.00000  0.000000
#>      1        49   47.50220  3.860776
#>      2        99  607.69584  6.409675
#>      3       135 1480.82436  7.300354
#>      4        96  529.56519  6.272056
#>      5        56   50.08795  3.913781
#>      6        33    0.80244 -0.220098
#> 
#> Test 2 score parameters (delta):
#> 
#>  Score Frequency      Delta Log_Delta
#>      0        39   1.000000  0.000000
#>      1        84  42.348644  3.745936
#>      2       120 221.281483  5.399436
#>      3       134 234.333759  5.456746
#>      4        73  34.813941  3.550018
#>      5        50   1.246199  0.220098
#> 
#> Total score parameters (sigma = Test1 + Test2):
#> 
#>  Score Frequency        Sigma Log_Sigma
#>      0        19      1.00000  0.000000
#>      1        22     89.85085  4.498151
#>      2        30   2840.63124  7.951782
#>      3        46  37961.61108 10.544331
#>      4        59 208878.48986 12.249508
#>      5        70 494214.10151 13.110724
#>      6        69 487527.55199 13.097102
#>      7        62 187523.16197 12.141658
#>      8        51  32196.51658 10.379614
#>      9        23   2591.74144  7.860085
#>     10        26     90.35566  4.503754
#>     11        23      1.00000  0.000000
#> 
#> Converged:  TRUE (after 584 iterations)
#> 
#> --- Goodness of Fit ---
#> 
#> 1. Likelihood Ratio Test:
#>    LR = 18.30  DF = 20  p = 0.5676
#> 
#> 2. Goodman-Kruskal Gamma Test (one-sided):
#>    Gamma (observed) = 0.7738
#>    Gamma (expected) = 0.7723
#>    SE = 1897.7131
#>    Z = -0.00  p = 0.5000
#> 
#> 3. Orbit Analysis:
#>    Use analyze_orbits() to assess person fit within total score strata
summary(fit)
#> Leunbach Score Parameter Estimation - Summary
#> ==============================================
#> Power Series Distribution with Generalized Symmetric Functions
#> 
#> N = 500 observations
#> Test 1: scores 0 to 6 (observed:  0 to 6)
#> Test 2: scores 0 to 5 (observed: 0 to 5)
#> Total:  scores 0 to 11 (observed:  0 to 11)
#> 
#> === Goodness of Fit ===
#> 
#> 1. Likelihood Ratio Test (observed vs expected counts):
#>    Likelihood ratio G² =    18.30 (df = 20, p = 0.5676)
#>    Pearson chi-square  =    17.31
#> 
#> 2. Goodman-Kruskal Gamma Test (one-sided):
#>    Tests if observed correlation exceeds expected under the model.
#>    Gamma (observed)    =   0.7738
#>    Gamma (expected)    =   0.7723
#>    Standard error      = 1897.7131
#>    Z statistic         =    -0.00 (p = 0.5000, one-sided)
#> 
#> 3. Orbit Analysis (person fit):
#>    Run analyze_orbits() separately to assess the number of cases
#>    outside 95% confidence regions of orbit distributions.
#> 
#> Converged: TRUE (after 584 iterations)
```
