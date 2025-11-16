# Calculate Allele Frequencies

Calculates population allele frequencies from data, adjusting for COI.

## Usage

``` r
calcAfreq(dsmp, coi, tol = 1e-04, qstart = 0.5)
```

## Arguments

- dsmp:

  a list with each element corresponding to one sample.

- coi:

  a vector containing complexity of infection for each sample.

- tol:

  convergence tolerance for frequency estimates.

- qstart:

  a starting value for frequencies.

## Value

A list of allele frequencies, where each element is a numeric vector
containing frequencies for a single locus.

## Examples

``` r
coi   <- getCOI(dsmp, lrank = 2)           # estimate COI first
afreq <- calcAfreq(dsmp, coi, tol = 1e-5)
```
