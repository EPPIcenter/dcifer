# Logarithms of `reval`

Calculates logarithms of `reval` and `1 - reval`, as well as other
associated quantities.

## Usage

``` r
logReval(reval, M = NULL, neval = NULL, equalr = FALSE)
```

## Arguments

- reval:

  a matrix representing a grid of (r₁, ..., r_(M)) combinations, over
  which the likelihood will be calculated. Each column is a single
  combination.

- M:

  the number of related pairs of strains.

- neval:

  the number of relatedness values/combinations to evaluate over.

- equalr:

  a logical value. If `TRUE`, the same level of relatedness is assumed
  for M pairs of strains (r₁ = ... = r_(M)).

## Value

A list of length 5 that contains `log(reval)`, `log(1 - reval)`, the
number of `reval = 1` for each column, the number of `0 < reval < 1` for
each column, and `sum(log(1 - reval[reval < 1]))` for each column.

## Details

For `equalr = TRUE` relatedness estimation, `reval` should be a
`1 x neval` matrix.

## Examples

``` r
reval <- generateReval(M = 2, nr = 1e2)
logr  <- logReval(reval, M = 2, equalr = FALSE)

reval <- generateReval(M = 1, nr = 1e3)
logr3  <- logReval(reval, M = 3, equalr = TRUE)
logr1  <- logReval(reval, M = 1)
all(logr3$sum1r == logr1$sum1r*3)
#> [1] TRUE
```
