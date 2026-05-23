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

## Value

Invisibly returns `NULL`; called for its side effect of printing a
diagnostic table.

## Examples

``` r
set.seed(123)
n <- 400
theta <- rnorm(n)
test1 <- pmin(pmax(round(3 + 1.5 * theta + rnorm(n, sd = 0.8)), 0), 6)
test2 <- pmin(pmax(round(2.5 + 1.3 * theta + rnorm(n, sd = 0.7)), 0), 5)
fit <- leunbach_ipf(data.frame(test1, test2),
                    max_score1 = 6, max_score2 = 5)
diagnose_equating(fit, direction = "1to2")
#> Diagnostic for equating
#> =======================
#> 
#> Method: optimize
#> Source range: 0 to 6
#> Target range:  0 to 5
#> 
#> Score  Theta (raw)    Theta (log)    Expected Target
#> -----------------------------------------------------
#>     0            NA            NA            NA
#>     1      0.072285     -2.627143          0.87
#>     2      0.289082     -1.241046          1.73
#>     3      1.041972      0.041115          2.52
#>     4      4.147048      1.422397          3.36
#>     5     16.019347      2.773797          4.32
#>     6            NA            NA            NA
```
