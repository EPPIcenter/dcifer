#' Read and Reformat Data
#'
#' @param sfile,afile names of the files containing sample data and, optionally,
#'   population allele frequencies.
#' @param svar,lvar,avar,fvar variable names for sample ID, marker/locus,
#'   allele/haplotype, and population allele frequency.
#' @param ... additional arguments for \code{read.csv()}.
#' @return A list with samples and population allele frequencies.
#' @export

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

#' Generate a grid of parameter values to evaluate over
#'
#' @param M an integer.
#' @param rval \eqn{{r}} values for the grid. Takes precedence over \code{nr}.
#' @param nr an integer. If \code{rval} is not provided, it will be generated
#'   using \code{0}, \code{1}, and \code{nr - 1} values between them.
#' @return A a matrix with \code{M} rows and \code{nr + 1} or
#'   \code{length(rval)} columns.
#' @export
#'
generateReval <- function(M, rval = NULL, nr = NULL) {
  if (is.null(rval)) {
    if (is.null(nr)) {
      return(NA)
    }
    rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
  }
  if (M == 1) {
    return(matrix(rval, 1))
  }
  reval <- as.matrix(expand.grid(rep(list(rval), M)))
  for (k in 1:nrow(reval)) {         # faster than apply()
    reval[k, ] <- sort(reval[k, ])
  }
  return(t(unique(reval)))
}

# get combinations of alleles when length(Ux) > coix and upcoi = FALSE
getComb <- function(Ux, coix) {
  if (length(Ux) <= coix) {
    return(matrix(Ux, ncol = 1))
  }
  return(utils::combn(Ux, coix))
}
