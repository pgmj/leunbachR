# Convert contingency table to data frame

Converts a contingency table (matrix of counts) to a data frame with one
row per observation, suitable for input to leunbach_ipf().

## Usage

``` r
table_to_data(tab, test1_scores, test2_scores)
```

## Arguments

- tab:

  Contingency table (matrix)

- test1_scores:

  Score values for Test 1

- test2_scores:

  Score values for Test 2

## Value

A data frame with columns test1 and test2
