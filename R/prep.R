
#' Read and Reformat Data
#'
#' \code{readDat} and \code{readAfreq} read data and population allele
#' frequencies from \code{csv} files and reformat them for further processing.
#' \code{formatDat} and \code{formatAfreq} reformat corresponding data frames.
#' Original data are assumed to be in a long format, with one row per allele.
#'
#' @name  format
#' @param sfile  the name of the file containing sample data.
#' @param dlong  a data frame containing sample data.
#' @param afile  the name of the file containing population allele frequencies.
#' @param aflong a data frame containing population allele frequencies.
#' @param svar   the name of the variable for sample ID.
#' @param lvar   the name of the variable for locus/marker.
#' @param avar   the name of the variable for allele/haplotype.
#' @param fvar   the name of the variable for population allele frequiency.
#' @param ... additional arguments for \code{read.csv()}.
#'
#' @return For \code{readDat} and \code{formatDat}, a list with elements
#'   corresponding to samples. Each element of the list is itself a list of
#'   binary vectors, one vector for each locus. For \code{readAfreq} and
#'   \code{formatAfreq}, a list with elements corresponding to loci. The
#'   frequencies at each locus are normalized and sum to 1. Samples, loci, and
#'   alleles are ordered by their IDs/names.
#'
#' @examples
#' sfile <- system.file("extdata", "MozParagon.csv", package = "dcifer")
#' dsmp  <- readDat(sfile, svar = "sampleID", lvar = "locus", avar = "allele")
#'
#' # OR, if the dataset is provided as an R data frame, e.g.
#' dlong <- read.csv(sfile)
#' # reformat only:
#' dsmp <- formatDat(dlong, svar = "sampleID", lvar = "locus", avar = "allele")
#'
#' afile  <- system.file("extdata", "MozAfreq.csv", package = "dcifer")
#' afreq2 <- readAfreq(afile, lvar = "locus", avar = "allele", fvar = "freq")
#'
#' # OR, if allele frequencies are provided as an R data frame, e.g.
#' aflong <- read.csv(afile)
#' # reformat only:
#' afreq2 <- formatAfreq(aflong, lvar = "locus", avar = "allele", fvar = "freq")
#'
#' dsmp2  <- matchAfreq(dsmp, afreq2)
#'
#' @seealso \code{\link{matchAfreq}} for making sure that the list containing
#'   sample data is conformable to provided population allele frequencies.

#' @rdname format
#' @export
#'
readDat <- function(sfile, svar, lvar, avar, ...) {
  dlong <- utils::read.csv(sfile, ...)
  return(formatDat(dlong, svar, lvar, avar))
}

#' @rdname format
#' @export
#'
formatDat <- function(dlong, svar, lvar, avar) {
  anames <- by(dlong, dlong[[lvar]], function(df) sort(unique(df[[avar]])))
  smps <- lapply(split(dlong, dlong[[svar]]),
                 function(df) split(df[[avar]], df[[lvar]]))
  nsmp <- length(smps)

  dsmp <- stats::setNames(vector("list", nsmp), names(smps))
  smp0 <- lapply(anames, function(vchar) {
    x <- rep(0, length(vchar))
    names(x) <- vchar
    x
  })
  for (ismp in 1:nsmp) {
    smpi <- smp0
    for (iloc in 1:length(anames)) {
      asmp <- smps[[ismp]][[names(anames)[iloc]]]
      smpi[[iloc]][match(asmp, anames[[iloc]])] <- 1
    }
    dsmp[[ismp]] <- smpi
  }
  return(dsmp)
}

#' @rdname format
#' @export
#'
readAfreq <- function(afile, lvar, avar, fvar, ...) {
  aflong <- utils::read.csv(afile, ...)
  return(formatAfreq(aflong, lvar, avar, fvar))
}

#' @rdname format
#' @export
#'
formatAfreq <- function(aflong, lvar, avar, fvar) {
  lapply(split(aflong, aflong[[lvar]]), function(df) {
    x <- stats::setNames(df[[fvar]], df[[avar]])
    x/sum(x)
  })
}

#' Match Samples and Allele Frequencies
#'
#' Checks if the list containing sample data is conformable to provided
#' population allele frequencies and reformats it if needed.
#'
#' @details The function reorders loci and alleles in \code{dsmp} to match those
#'   in \code{afreq} and inserts alleles into \code{dsmp} if they are present in
#'   \code{afreq} and not in \code{dsmp}; doesn't handle cases when alleles are
#'   present in \code{dsmp} but not in \code{afreq}. Allele names are required
#'   for this procedure.
#' @inheritParams ibdDat
#' @inheritParams ibdPair
#' @return A list of the same length as \code{dsmp}, with each element matching
#'   the lengths and the names of \code{afreq} and its elements.
#'
#' @examples
#' afile  <- system.file("extdata", "MozAfreq.csv", package = "dcifer")
#' afreq2 <- readAfreq(afile, lvar = "locus", avar = "allele", fvar = "freq")
#' dsmp2  <- matchAfreq(dsmp, afreq2)
#'
#' @seealso \code{\link{readDat}} and \code{\link{readAfreq}} for reading in and
#'   reformating data.
#' @export
#'
matchAfreq <- function(dsmp, afreq) {
  if (!all(sapply(afreq, length) == sapply(dsmp[[1]], length)) ||
      !all(names(unlist(afreq)) == names(unlist(dsmp[[1]])))) {
    loci <- names(afreq)
    K    <- sapply(afreq, length)
    nloc <- length(K)
    for (ismp in 1:length(dsmp)) {
      dsmpi <- stats::setNames(dsmp[[ismp]], loci)
      for (iloc in 1:nloc) {
        old <- dsmp[[ismp]][[loci[iloc]]]
        upd <- afreq[[iloc]]
        upd[] <- 0
        upd[names(old[old == 1])] <- 1
        if (length(upd) > K[iloc]) {
          stop("extra alleles in dsmp")
        }
        dsmpi[[iloc]] <- upd
      }
      dsmp[[ismp]] <- dsmpi
    }
  }
  return(dsmp)
}

#' Grid of Parameter Values
#'
#' Generates a grid of parameter values to evaluate over.
#'
#' @param M an integer.
#' @param rval a numeric vector with parameter values for the grid. If not
#'   provided, it will be generated by dividing \ifelse{html}{\out{[0,
#'   1]}}{\eqn{[0, 1]}} into \code{nr} equal intervals.
#' @param nr an integer. Ignored if \code{rval} is provided.
#'
#' @return A a matrix with \code{M} rows and \code{nr + 1} or
#' \code{length(rval)} columns.
#'
#' @examples
#' reval <- generateReval(M = 2, nr = 1e2)
#' dim(reval)
#'
#' @export
#'
generateReval <- function(M, rval = NULL, nr = NULL) {
  if (is.null(rval)) {
    rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
  }
  if (M == 1) {
    return(matrix(rval, 1))
  }
  reval <- as.matrix(expand.grid(rep(list(rval), M)))
  colnames(reval) <- NULL
  for (k in 1:nrow(reval)) {         # faster than apply()
    reval[k, ] <- sort(reval[k, ])
  }
  return(t(unique(reval)))
}

#' Calculate COI
#'
#' Calculates complexity of infection for a list of samples, using the number of
#' detected alleles.
#'
#' @inheritParams ibdDat
#' @param lrank the rank of the locus that will determine a sample's COI (loci
#'   are ranked by the number of detected alleles).
#' @return a vector with estimated COI for each sample.
#' @examples
#' coi <- getCOI(dsmp, lrank = 2)
#' @export
#'
getCOI <- function(dsmp, lrank = 2) {
  sapply(dsmp, function(smp) sort(sapply(smp, sum), decreasing = TRUE)[lrank])
}

#' Calculate Allele Frequencies
#'
#' Calculates population allele frequencies from data, adjusting for COI.
#'
#' @inheritParams ibdDat
#' @param tol    convergence tolerance for frequency estimates.
#' @param qstart a starting value for frequencies.
#' @return A list of allele frequencies, where each element is a numeric vector
#'   containing frequencies for a single locus.
#' @examples
#' coi   <- getCOI(dsmp, lrank = 2)           # estimate COI first
#' afreq <- calcAfreq(dsmp, coi, tol = 1e-5)
#' @export
#'
calcAfreq <- function(dsmp, coi, tol = 1e-4, qstart = 0.5) {
  alist <- dsmp[[1]]
  K     <- sapply(alist, length)
  nloc  <- length(K)
  dloc <- unlist(dsmp, recursive = FALSE)
  dloc <- apply(matrix(dloc, nloc), 1, as.list)
  off  <- 10^(-ceiling(log(sum(coi), 10)))
  for (iloc in 1:nloc) {
    b <- do.call(rbind, dloc[[iloc]])
    for (iall in 1:K[iloc]) {
      q <- qNewton(b[, iall], coi, tol = tol, qstart = qstart, off = off)
      alist[[iloc]][iall] <- 1 - q
    }
  }
  sumaf <- sapply(alist, sum)
  return(mapply(function(v, sum) v/sum, alist, sumaf, SIMPLIFY = FALSE))
}

qNewton <- function(b, coi, qstart = 0.5, tol = 1e-4, toldrv = 1e-3,
                    off = 1e-3) {
  i1 <- as.logical(b)
  n1 <- coi[i1]
  sum0 <- sum(coi[!i1])
  q    <- 2
  qdrv <- 100
  qnew <- qstart
  while(abs(qnew - q) > tol || abs(qdrv[1]) > toldrv) {
    q     <- qnew
    qdrv  <- qDrvs(q, n1, sum0, drv2 = TRUE)
    qnew  <- q - qdrv[1]/qdrv[2]
    if (qnew <= 0) {
      if (qDrvs(off, n1, sum0, drv2 = FALSE) < 0) {
        return(0)
      }
      qnew <- off
    } else if (qnew >= 1) {
      if (qDrvs(1, n1, sum0, drv2 = FALSE) > 0) {
        return(1)
      }
      qnew <- 1 - off
    }
  }
  return(qnew)
}

qDrvs <- function(q, n1, sum0, drv2 = FALSE) {
  qn1   <- q^(n1 - 1)
  oneqn <- 1 - q*qn1  # 1 - q^n1
  qdrv1 <- -sum(n1*qn1/oneqn) + sum0/q
  if (!drv2) {
    return(qdrv1)
  }
  qdrv2 <- -sum(n1*((n1 - 1)*qn1/q*oneqn + qn1^2*n1)/oneqn^2) - sum0/q^2
  return(c(qdrv1, qdrv2))
}
