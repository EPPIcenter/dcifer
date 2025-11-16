# Read and Reformat Data

`readDat` and `readAfreq` read data and population allele frequencies
from `csv` files and reformat them for further processing. `formatDat`
and `formatAfreq` reformat corresponding data frames. Original data are
assumed to be in a long format, with one row per allele.

## Usage

``` r
readDat(sfile, svar, lvar, avar, ...)

formatDat(dlong, svar, lvar, avar)

readAfreq(afile, lvar, avar, fvar, ...)

formatAfreq(aflong, lvar, avar, fvar)
```

## Arguments

- sfile:

  the name of the file containing sample data.

- svar:

  the name of the variable for sample ID.

- lvar:

  the name of the variable for locus/marker.

- avar:

  the name of the variable for allele/haplotype.

- ...:

  additional arguments for
  [`read.csv()`](https://rdrr.io/r/utils/read.table.html).

- dlong:

  a data frame containing sample data.

- afile:

  the name of the file containing population allele frequencies.

- fvar:

  the name of the variable for population allele frequiency.

- aflong:

  a data frame containing population allele frequencies.

## Value

For `readDat` and `formatDat`, a list with elements corresponding to
samples. Each element of the list is itself a list of binary vectors,
one vector for each locus. For `readAfreq` and `formatAfreq`, a list
with elements corresponding to loci. The frequencies at each locus are
normalized and sum to 1. Samples, loci, and alleles are ordered by their
IDs/names.

## See also

[`matchAfreq`](https://eppicenter.github.io/dcifer/reference/matchAfreq.md)
for making sure that the lists containing sample data and provided
population allele frequencies are matching.

## Examples

``` r
sfile <- system.file("extdata", "MozParagon.csv", package = "dcifer")
dsmp  <- readDat(sfile, svar = "sampleID", lvar = "locus", avar = "allele")

# OR, if the dataset is provided as an R data frame, e.g.
dlong <- read.csv(sfile)
# reformat only:
dsmp <- formatDat(dlong, svar = "sampleID", lvar = "locus", avar = "allele")

afile <- system.file("extdata", "MozAfreq.csv", package = "dcifer")
afreq <- readAfreq(afile, lvar = "locus", avar = "allele", fvar = "freq")

# OR, if allele frequencies are provided as an R data frame, e.g.
aflong <- read.csv(afile)
# reformat only:
afreq <- formatAfreq(aflong, lvar = "locus", avar = "allele", fvar = "freq")

da_upd <- matchAfreq(dsmp, afreq)
dsmp2  <- da_upd$dsmp
afreq2 <- da_upd$afreq
```
