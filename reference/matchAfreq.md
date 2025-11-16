# Match Samples and Allele Frequencies

Checks if the list containing sample data is conformable to provided
population allele frequencies and reformats it if needed.

## Usage

``` r
matchAfreq(dsmp, afreq, minfreq = 1e-04)
```

## Arguments

- dsmp:

  a list with each element corresponding to one sample.

- afreq:

  a list of allele frequencies. Each element of the list corresponds to
  a locus.

- minfreq:

  an allele frequency to assign to alleles that are present in `dsmp`
  but not in `afreq`.

## Value

A named list of length 2 containing updated and matching `dsmp` and
`afreq`.

## Details

The function reorders loci and alleles in `dsmp` to match those in
`afreq` and inserts alleles into `dsmp` if they are present in `afreq`
and not in `dsmp`; doesn't handle cases when alleles are present in
`dsmp` but not in `afreq`. Allele names are required for this procedure.

## See also

[`readDat`](https://eppicenter.github.io/dcifer/reference/format.md) and
[`readAfreq`](https://eppicenter.github.io/dcifer/reference/format.md)
for reading in and reformating data.

## Examples

``` r
afile  <- system.file("extdata", "MozAfreq.csv", package = "dcifer")
afreq  <- readAfreq(afile, lvar = "locus", avar = "allele", fvar = "freq")
da_upd <- matchAfreq(dsmp, afreq, minfreq = 1e-3)
dsmp2  <- da_upd$dsmp
afreq2 <- da_upd$afreq
```
