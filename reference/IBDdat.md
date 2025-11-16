# Pairwise Relatedness

Provides pairwise relatedness estimates within a dataset or between two
datasets along with optional p-values and confidence intervals (CI).

## Usage

``` r
ibdDat(
  dsmp,
  coi,
  afreq,
  dsmp2 = NULL,
  coi2 = NULL,
  pval = TRUE,
  confint = FALSE,
  rnull = 0,
  side = c("right", "left", "two-sided"),
  alpha = 0.05,
  nr = 1000,
  reval = NULL
)
```

## Arguments

- dsmp:

  a list with each element corresponding to one sample.

- coi:

  a vector containing complexity of infection for each sample.

- afreq:

  a list of allele frequencies. Each element of the list corresponds to
  a locus.

- dsmp2:

  a list representing a second dataset.

- coi2:

  a vector with complexities of infection for a second dataset.

- pval:

  a logical value specifying if p-values should be returned.

- confint:

  a logical value specifying if confidence intervals should be returned.

- rnull:

  a null value of relatedness parameter for hypothesis testing (needed
  if `pval = TRUE`).

- side:

  a character string specifying if a one-sided (`"right"` or `"left"`)
  or a two-sided (`"two-sided"`) hypothesis test should be performed
  (needed if `pval = TRUE`). Set to `"right"` if `rnul = 0` and to
  `"left"` if `rnull = 1`.

- alpha:

  significance level for a 1 - α confidence region.

- nr:

  an integer specifying precision of the estimate: resolution of a grid
  of parameter values (\[0, 1\] divided into `nr` equal intervals), over
  which the likelihood will be calculated. Ignored if non-null `reval`
  is provided.

- reval:

  a vector or a single-row matrix. A grid of parameter values, over
  which the likelihood will be calculated.

## Value

A matrix if `pval` and `confint` are `FALSE` and 3-dimensional arrays
otherwise. The matrices are lower triangular if distances are calculated
within a dataset. For a 3-dimensional array, stacked matrices contain
relatedness estimates, p-values, and endpoints of confidence intervals
(if requested).

## Details

For this function, M is set to 1. If `confint = FALSE`, Newton's method
is used to find the estimates, otherwise the likelihood is calculated
for a grid of parameter values.

## See also

[`ibdPair`](https://eppicenter.github.io/dcifer/reference/IBDpair.md)
for genetic relatedness between two samples and
[`ibdEstM`](https://eppicenter.github.io/dcifer/reference/ibdEstM.md)
for estimating the number of related pairs of strains.

## Examples

``` r
coi   <- getCOI(dsmp, lrank = 2)           # estimate COI
afreq <- calcAfreq(dsmp, coi, tol = 1e-5)  # estimate allele frequencies

# subset of samples for faster processing
i1 <- 1:15     # from Maputo
i2 <- 31:40    # from Inhambane
isub <- c(i1, i2)

# matrix is returned
dres1 <- ibdDat(dsmp[isub], coi[isub], afreq, pval = FALSE)
dim(dres1)
#> [1] 25 25

# test a null hypothesis H0: r = 0, change precision
dres2 <- ibdDat(dsmp[isub], coi[isub], afreq, pval = TRUE, rnull = 0,
                nr = 1e2)
dim(dres2)
#> [1] 25 25  2

# test H0: r <= 0.2, include 99% confidence intervals
dres3 <- ibdDat(dsmp[isub], coi[isub], afreq, pval = TRUE, confint = TRUE,
                rnull = 0.2, side = "right", alpha = 0.01)
dres3[2, 1, ]
#>  estimate   p_value  CI_lower  CI_upper 
#> 0.0000000 0.9970094 0.0000000 0.1810000 

# pairwise relatedness between two datasets, H0: r = 0
drbetween <- ibdDat(dsmp[i1], coi[i1], afreq,
                    dsmp2 = dsmp[i2], coi2 = coi[i2])
dim(drbetween)
#> [1] 15 10  2
drbetween[1, 2, ]
#>   estimate    p_value 
#> 0.05595357 0.19681441 
sum(is.na(drbetween[, , 1]))
#> [1] 0
```
