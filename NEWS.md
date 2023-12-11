
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
