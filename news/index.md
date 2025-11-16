# Changelog

## dcifer 1.5.0

CRAN release: 2025-11-11

- New: a function
  [`mixMat()`](https://eppicenter.github.io/dcifer/reference/mixMat.md)
  that combines two triangular matrices into one square matrix (useful
  for plotting).
- Changes to
  [`plotRel()`](https://eppicenter.github.io/dcifer/reference/plotRel.md)
  function:
- Highlighted (statistically significant) entries are now specified by a
  logical square matrix
- These entries can be either outlined in a different color or
  “amplified” (made larger than other entries), which is useful for
  large datasets

## dcifer 1.4.0

- Added: one-sided (one-tailed) tests for functions
  [`ibdPair()`](https://eppicenter.github.io/dcifer/reference/IBDpair.md),
  [`ibdDat()`](https://eppicenter.github.io/dcifer/reference/IBDdat.md),
  and
  [`ibdEstM()`](https://eppicenter.github.io/dcifer/reference/ibdEstM.md).
  The argument `side` indicates if a one- or two-sided test is performed
  and the side for the one-sided test.
- [`matchAfreq()`](https://eppicenter.github.io/dcifer/reference/matchAfreq.md)
  function now provides both updated `dsmp` and `afreq` (outputting a
  list of two corresponding elements). The new version allows for
  “extra” alleles in both `dsmp` and `afreq`.
- Minor updates/fixes:
  - NA returned when there are no loci with data in both samples
  - names for COI match sample ID’s

## dcifer 1.3.0

- Function
  [`calcAfreq()`](https://eppicenter.github.io/dcifer/reference/calcAfreq.md)
  adjusted to exclude missing data (per ID at a locus)

## dcifer 1.2.1

CRAN release: 2023-12-12

- Single locus data now allowed.

## dcifer 1.2.0

CRAN release: 2022-10-31

- New: functions
  [`formatDat()`](https://eppicenter.github.io/dcifer/reference/format.md)
  and
  [`formatAfreq()`](https://eppicenter.github.io/dcifer/reference/format.md)
  that reformat data frames containing sample data and population allele
  frequencies in long formats.

## dcifer 1.1.1

CRAN release: 2022-08-10

- Fix:
  [`logReval()`](https://eppicenter.github.io/dcifer/reference/logReval.md)
  now covers all relevant edge cases.

## dcifer 1.1.0

CRAN release: 2022-08-06

- Special case added: loci with no shared alleles (more efficient);
- [`logReval()`](https://eppicenter.github.io/dcifer/reference/logReval.md)
  now has an argument `equalr` and performs additional checks. Can
  provide `equalr = TRUE` and M \> 1; the result is reflected in the
  `sum1r` element of the output.

## dcifer 1.0.1

CRAN release: 2022-07-15

First public release.
