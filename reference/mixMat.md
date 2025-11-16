# Combine two matrices

Creates a single matrix from two lower or upper triangular matrices.

## Usage

``` r
mixMat(mat1, mat2, lower = TRUE)
```

## Arguments

- mat1, mat2:

  square matrices of the same dimensions.

- lower:

  a logical value indicating if the resulting matrix should be comprised
  of lower triangular matrices.

## Value

A square matrix.

## See also

[`plotRel`](https://eppicenter.github.io/dcifer/reference/plotRel.md)
and
[`plotColorbar`](https://eppicenter.github.io/dcifer/reference/plotColorbar.md)
for plotting relatedness estimates.

## Examples

``` r
x <- matrix(1:16, 4)
y <- matrix(101:116, 4)
x[upper.tri(x, diag = TRUE)] <- NA
mixMat(x, y)
#>      [,1] [,2] [,3] [,4]
#> [1,]   NA  102  103  104
#> [2,]    2   NA  107  108
#> [3,]    3    7   NA  112
#> [4,]    4    8   12   NA
```
