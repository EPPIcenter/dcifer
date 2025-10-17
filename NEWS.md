# dcifer 1.4.0

- Added: one-sided (one-tailed) tests for functions `ibdPair()`, `ibdDat()`, and `ibdEstM()`. The argument `side` indicates if a one- or two-sided test is performed and the side for the one-sided test. 
- Minor updates/fixes: 
  * NA returned when there are no loci with data in both samples
  * names for COI match sample ID's

# dcifer 1.3.0

- Function `calcAfreq()` adjusted to exclude missing data (per ID at a locus)

# dcifer 1.2.1

- Single locus data now allowed. 

# dcifer 1.2.0

- New: functions `formatDat()` and `formatAfreq()` that reformat data frames containing sample data and population allele frequencies in long formats.

# dcifer 1.1.1

- Fix: `logReval()` now covers all relevant edge cases. 

# dcifer 1.1.0

- Special case added: loci with no shared alleles (more efficient);
- `logReval()` now has an argument `equalr` and performs additional checks. Can provide `equalr = TRUE` and M > 1; the result is reflected in the `sum1r` element of the output. 

# dcifer 1.0.1

First public release.
