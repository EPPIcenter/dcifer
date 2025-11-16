# Calculate COI

Calculates complexity of infection for a list of samples, using the
number of detected alleles.

## Usage

``` r
getCOI(dsmp, lrank = 2)
```

## Arguments

- dsmp:

  a list with each element corresponding to one sample.

- lrank:

  the rank of the locus that will determine a sample's COI (loci are
  ranked by the number of detected alleles).

## Value

A vector with estimated COI for each sample.

## Examples

``` r
coi <- getCOI(dsmp, lrank = 2)
```
