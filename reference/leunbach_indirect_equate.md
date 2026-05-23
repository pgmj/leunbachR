# Indirect Equating via an Anchor Test

Performs indirect equating from Test A to Test C via an anchor Test B.
This chains two direct equatings: A -\> B and B -\> C.

## Usage

``` r
leunbach_indirect_equate(
  fit_ab,
  fit_bc,
  direction_ab = c("1to2", "2to1"),
  direction_bc = c("1to2", "2to1"),
  method = c("optimize", "newton"),
  verbose = FALSE
)
```

## Arguments

- fit_ab:

  A leunbach_ipf object for the A-B equating (Tests A and B)

- fit_bc:

  A leunbach_ipf object for the B-C equating (Tests B and C)

- direction_ab:

  Direction for A-B equating: "1to2" or "2to1"

- direction_bc:

  Direction for B-C equating: "1to2" or "2to1"

- method:

  Optimization method: "optimize" (default) or "newton"

- verbose:

  Print detailed output

## Value

A list of class "leunbach_indirect" containing:

- equating_table: Data frame with source scores, expected equated
  scores, and rounded scores

- eq_ab: Direct equating object for A -\> B

- eq_bc: Direct equating object for B -\> C

- fit_ab: Original leunbach_ipf object for A-B

- fit_bc: Original leunbach_ipf object for B-C

## Details

Indirect equating works by chaining two direct equatings:

1.  For a score x on Test A, find the expected score on Test B
    (typically non-integer)

2.  Find expected Test C scores for the integer B scores below and above

3.  Interpolate to get the expected Test C score for the non-integer B
    score

4.  Round to get the equated integer score

## Examples

``` r
# Simulate scores for three tests, with B serving as the anchor.
# Group 1 takes tests A and B; group 2 takes tests B and C.
set.seed(123)
n <- 400
theta1 <- rnorm(n)
a <- pmin(pmax(round(3 + 1.5 * theta1 + rnorm(n, sd = 0.8)), 0), 6)
b1 <- pmin(pmax(round(2.5 + 1.3 * theta1 + rnorm(n, sd = 0.7)), 0), 5)

theta2 <- rnorm(n)
b2 <- pmin(pmax(round(2.5 + 1.3 * theta2 + rnorm(n, sd = 0.7)), 0), 5)
cc <- pmin(pmax(round(3 + 1.4 * theta2 + rnorm(n, sd = 0.8)), 0), 6)

fit_ab <- leunbach_ipf(data.frame(a, b1), max_score1 = 6, max_score2 = 5)
fit_bc <- leunbach_ipf(data.frame(b2, cc), max_score1 = 5, max_score2 = 6)

# Indirect equating: A -> B -> C
indirect <- leunbach_indirect_equate(fit_ab, fit_bc,
                                     direction_ab = "1to2",
                                     direction_bc = "1to2")
print(indirect)
#> Leunbach Indirect Equating
#> ==========================
#> 
#> Path: Test A -> Test B -> Test C
#> Method: optimize
#> 
#> Source (Test A) range:  0 to 6
#> Anchor (Test B) range: 0 to 5
#> Target (Test C) range: 0 to 6
#> 
#>  Score_Test A Expected_Test C Rounded_Test C
#>             0            0.00              0
#>             1            1.19              1
#>             2            2.09              2
#>             3            2.89              3
#>             4            3.76              4
#>             5            4.92              5
#>             6            6.00              6
```
