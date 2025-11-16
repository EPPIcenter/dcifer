# Estimate Relatedness and a Number of Related Strains

Estimates the number of related pairs of strains between two infections
along with corresponding relatedness estimates and optional inference.

## Usage

``` r
ibdEstM(
  pair,
  coi,
  afreq,
  Mmax = 6,
  pval = FALSE,
  confreg = FALSE,
  llik = FALSE,
  rnull = 0,
  side = c("right", "left", "two-sided"),
  alpha = 0.05,
  equalr = FALSE,
  freqlog = FALSE,
  nrs = c(1000, 100, 32, 16, 12, 10),
  revals = NULL,
  tol0 = 1e-09,
  logrs = NULL,
  nevals = NULL,
  nloc = NULL
)
```

## Arguments

- pair:

  a list of length two containing data for a pair of samples.

- coi:

  a vector containing complexity of infection for each sample.

- afreq:

  a list of allele frequencies. Each element of the list corresponds to
  a locus.

- Mmax:

  a maximum number of related pairs of strains to evaluate over. If
  greater than `min(coi)`, will be set to `min(coi)`.

- pval, confreg, llik:

  logical values specifying if p-value, confidence region, and
  log-likelihood for a range of \\r\\ values should be returned.

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

- equalr:

  a logical value. If `TRUE`, the same level of relatedness is assumed
  for M pairs of strains (r₁ = ... = r_(M)).

- freqlog:

  a logical value indicating if `afreq` is on the log scale.

- nrs:

  an integer vector where `i`'th element correspons to M = i and
  indicates precision of the estimate (resolution of a grid of parameter
  values). Ignored if non-null `revals` is provided.

- revals:

  a list where `i`'th element corresponds to M = i and is a matrix
  representing a grid of parameter values (a matrix where each column
  represents a single (r₁, ..., r_(M)) combination).

- tol0:

  a tolerance value for an estimate to be considered zero.

- logrs:

  a list where `i`'th element corresponds to M = i and is a list as
  returned by
  [`logReval`](https://eppicenter.github.io/dcifer/reference/logReval.md).

- nevals:

  a vector where `i`'th element corresponds to M = i and provides the
  number of relatedness values/combinations to evaluate over.

- nloc:

  the number of loci.

## Value

A named list if multiple output logical values are `TRUE` - or a vector
if only `rhat = TRUE`. The output includes:

- a relatedness estimate (numeric vector of length corresponding to the
  estimated number of related pairs);

- a p-value if `pval = TRUE`;

- parameter values from the grid in `revals` that are within the
  confidence region if `confreg = TRUE`;

- log-likelihood values for the parameter grid in `revals` if
  `llik = TRUE`.

## See also

[`ibdPair`](https://eppicenter.github.io/dcifer/reference/IBDpair.md)
for estimates of relatedness between two samples and
[`ibdDat`](https://eppicenter.github.io/dcifer/reference/IBDdat.md) for
pairwise relatedness estimates within a dataset or between two datasets.

## Examples

``` r
coi   <- getCOI(dsmp, lrank = 2)           # estimate COI
afreq <- calcAfreq(dsmp, coi, tol = 1e-5)  # estimate allele frequencies

# two samples
ipair <- c(21, 17)
# for higher COI: c(33, 5): COI = 5-6; c(37, 20): 4-3, c(41, 50): 5-4

Mmax  <- min(coi[ipair])
# choose resolution of the grid for different M
nrs   <- c(1e3, 1e2, 32, 16, 12, 10)[1:Mmax]
revals <- mapply(generateReval, 1:Mmax, nr = nrs)

(res1 <- ibdEstM(dsmp[ipair], coi[ipair], afreq, Mmax = Mmax, equalr = FALSE,
                 reval = revals))
#> [1] 0.31 1.00
(res2 <- ibdEstM(dsmp[ipair], coi[ipair], afreq, Mmax = Mmax, equalr = TRUE))
#> [1] 0.673 0.673
# number of related pairs of strains (M')
sum(res1 > 0)
#> [1] 2
sum(res2 > 0)  # can be 0's
#> [1] 2
```
