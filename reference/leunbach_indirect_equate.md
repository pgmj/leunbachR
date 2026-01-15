# Indirect Equating via an Anchor Test

Performs indirect equating from Test A to Test C via an anchor Test B.
This chains two direct equatings: A → B and B → C.

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

- eq_ab: Direct equating object for A → B

- eq_bc: Direct equating object for B → C

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
# Fit models for A-B and B-C
fit_ab <- leunbach_ipf(data_ab)
#> Error: object 'data_ab' not found
fit_bc <- leunbach_ipf(data_bc)
#> Error: object 'data_bc' not found

# Indirect equating:  Test1 of fit_ab → Test2 of fit_ab → Test2 of fit_bc
indirect <- leunbach_indirect_equate(fit_ab, fit_bc, 
                                      direction_ab = "1to2", 
                                      direction_bc = "1to2")
#> Error: object 'fit_ab' not found
print(indirect)
#> Error: object 'indirect' not found
```
