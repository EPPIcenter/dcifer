#' Read and Reformat Data
#'
#' @param sfile,afile name of the files containing sample data and, optionally,
#'   population allele frequencies.
#' @param svar,lvar,avar,fvar variable names for sample ID, marker/locus,
#'   allele/haplotype, and population allele frequency.
#' @param ... additional arguments for \code{read.csv()}.

readDat <- function(sfile, afile = NULL, svar = "SampleID", lvar = "MarkerID",
                    avar = "Haplotype", fvar = "frequency", ...) {
  dat <- utils::read.csv(sfile, ...)
  if (is.null(afile)) {
    afreq <- by(dat, dat[[lvar]], function(df) table(df[[avar]])/nrow(df))
  } else {
    afreq <- utils::read.csv(afile, ...)
    afreq <- lapply(split(afreq, afreq[[lvar]]),
                    function(df) {
                      x <- stats::setNames(df[[fvar]], df[[avar]])
                      x/sum(x)
                    })
  }

  smps <- lapply(split(dat, dat[[svar]]),
                 function(df) split(df[[avar]], df[[lvar]]))
  nsmp <- length(smps)

  dsmp <- stats::setNames(vector("list", nsmp), names(smps))
  smp0 <- lapply(afreq, function(x) {x[] <- 0; x})
  for (ismp in 1:nsmp) {
    smpi <- smp0
    for (iloc in 1:length(afreq)) {
      aall <- names(afreq[[iloc]])
      asmp <- smps[[ismp]][[names(afreq)[iloc]]]
      smpi[[iloc]][match(asmp, aall)] <- 1
    }
    dsmp[[ismp]] <- smpi
  }
  return(list(dsmp = dsmp, afreq = afreq))
}

#' Rescale Allele Frequency
#' @description Handles small values of estimated population allele frequencies to protect against potential false positives driving the result.
#'
#' @param afreq a list with population allele frequencies.
#' @param dsmp  a list of lists with samples as top level elements.
#' @param rmrare a logical value. If \code{TRUE}, rare alleles will be removed from both \code{afreq} and \code{dsmp}; if \code{FALSE}, their frequencies will be increased to \code{lbound}.
#' @param lbound a bound for considering an allele "rare".
#' @return a list with updated and rescaled \code{afreq} and \code{dsmp} (updated if \code{rmrare = TRUE}).
#'
rescaleFreq <- function(afreq, dsmp = NULL, rmrare = FALSE, lbound = 0.01) {
  ismall <- lapply(afreq, function(x) which(x < lbound))
  if (rmrare) {
    afreq <- mapply(function(x, i) {x <- x[-i]; x/sum(x)}, afreq, ismall)
    for (ismp in 1:length(dsmp)) {
      dsmp[[ismp]] <- mapply(function(smp, i) smp[-i], dsmp[[ismp]], ismall)
    }
  } else {
    afreq <- mapply(function(x, i, lbound) {x[i] <- lbound; x/sum(x)},
                    afreq, ismall, lbound = lbound)
  }
  return(list(afreq = afreq, dsmp = dsmp))
}

