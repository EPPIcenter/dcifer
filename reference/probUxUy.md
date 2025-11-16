# Likelihood for U_(x), U_(y)

Calculates log-likelihood for a pair of samples at a single locus.

## Usage

``` r
probUxUy(
  Ux,
  Uy,
  nx,
  ny,
  probs,
  M,
  logj,
  factj,
  equalr = FALSE,
  mnewton = TRUE,
  reval = NULL,
  logr = NULL,
  neval = NULL
)
```

## Arguments

- Ux, Uy:

  sets of unique alleles for two samples at a given locus. Vectors of
  indices corresponding to ordered probabilities in `probs`.

- nx, ny:

  complexity of infection for two samples. Vectors of length 1.

- probs:

  a vector of population allele frequencies (on a log scale) at a given
  locus. It is not checked if frequencies on a regular scale sum to 1.

- M:

  the number of related pairs of strains.

- logj, factj:

  numeric vectors containing precalculated logarithms and factorials.

- equalr:

  a logical value. If `TRUE`, the same level of relatedness is assumed
  for M pairs of strains (r₁ = ... = r_(M)).

- mnewton:

  a logical value. If `TRUE`, the coefficients for using Newton's method
  will be calculated.

- reval:

  a matrix representing a grid of (r₁, ..., r_(M)) combinations, over
  which the likelihood will be calculated. Each column is a single
  combination.

- logr:

  a list of length 5 as returned by
  [`logReval`](https://eppicenter.github.io/dcifer/reference/logReval.md).

- neval:

  the number of relatedness values/combinations to evaluate over.

## Value

- If `mnewton = TRUE`, a vector of length 2 containing coefficients for
  fast likelihood calculation;

- If `mnewton = FALSE`, a vector of length `neval` containing
  log-likelihoods for a range of parameter values.

## Examples

``` r
Ux <- c(1, 3, 7)                       # detected alleles at locus t
Uy <- c(2, 7)
coi <- c(5, 6)
aft <- runif(7)                        # allele frequencies for locus t
aft <- log(aft/sum(aft))

logj  <- log(1:max(coi))
factj <- lgamma(0:max(coi) + 1)

# M = 2, equalr = FALSE
M <- 2
reval <- generateReval(M, nr = 1e2)
logr  <- logReval(reval, M = M)
llikt <- probUxUy(Ux, Uy, coi[1], coi[2], aft, M, logj, factj,
                  equalr = FALSE, logr = logr, neval = ncol(reval))
length(llikt)
#> [1] 5151

# M = 2, equalr = TRUE
reval <- matrix(seq(0, 1, 0.001), 1)
logr  <- logReval(reval, M = M, equalr = TRUE)
llikt <- probUxUy(Ux, Uy, coi[1], coi[2], aft, M, logj, factj,
                  equalr = TRUE, logr = logr, neval = ncol(reval))

# M = 1, mnewton = FALSE
M <- 1
reval <- matrix(seq(0, 1, 0.001), 1)
logr  <- logReval(reval, M = M)
llikt <- probUxUy(Ux, Uy, coi[1], coi[2], aft, M, logj, factj,
                  mnewton = FALSE, reval = reval, logr = logr,
                  neval = ncol(reval))

# M = 1, mnewton = TRUE
probUxUy(Ux, Uy, coi[1], coi[2], aft, M, logj, factj, mnewton = TRUE)
#> [1] 2.017199e-06 2.254247e-06
```
