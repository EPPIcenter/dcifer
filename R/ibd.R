
#' Genetic Distance for a Pair of Samples
#'
#' @description Calculates likelihood over a range of values for relatedness
#'   parameter \eqn{{r}} and provides an estimate for a pair of samples along
#'   with confidence regions based on a likelihood ratio test. For data with
#'   genotyping errors where estimated COI is less than |Ux|, various options
#'   are offered.
#'
#' @details Handling of irregular cases: - Allele with population frequency of
#'   \eqn{0} is present: locus is skipped (does not contribute any information).
#'   - Number of unique alleles at a locus is greater than COI: if \code{upcoi =
#'   TRUE}, COI will be increased for that locus only; otherwise likelihoods for
#'   all the subsets of cardinality COI are calculated and geometric mean is
#'   taken.
#'
#' @param pair   a list of length two containing two samples.
#' @param coi    a vector indicating complexity of infection for each sample.
#' @param afreq  a list of allele frequencies. Each element of the list
#'   corresponds to a locus.
#' @param nm     the number of related pairs of strains.
#' @param nr     an integer value for the resolution of the grid (\eqn{nr - 1}
#'   values between 0 and 1), over which the likelihood will be calculated.
#'   Ignored if non-null \code{reval} is provided.
#' @param rval  \eqn{{r}} values for the grid or for evaluation when
#'   \code{equalr} is \code{TRUE}. If \code{NULL}, will be evenly spaced between
#'   0 and 1 and interval \eqn{1/nr}.
#' @param reval  the grid of \eqn{{r}} combinations, over which the likelihood
#'   will be calculated. A matrix where each column represents a single
#' combination.
#' @param logr a list as returned by \code{logReval} with logs of \code{reval}
#'   and other quantities.
#' @param equalr a logical value. If \code{TRUE}, the same values of \eqn{r} are
#'   assumed for all \code{nm} pairs of related strains.
#' @param out    a character string for the type of results to be returned. If
#'   \code{"mle"}, an estimate is returned, if \code{"llik"} - a vector of
#'   log-likelihood for each combination.
#' @param alpha significance level.
#' @param freqlog a logical value indicating if \code{afreq} is on the log
#'   scale.
#' @param neval the number of relatedness values/combinations to evaluate over.
#' @param nloc  the number of loci.
#' @param upcoi method to handle cases when \eqn{COI < |Ux|} (see @details).
#'
#' @return
#'   * If \code{out = "mle"}: a vector of length 1 if \code{equalr = TRUE}
#'   or of length \code{nm} otherwise, containing estimated \eqn{{r}} (or vector
#'   /matrix if not just first).
#'   * If \code{out = "llike"}, a vector of log-likelihood values for each
#'   evaluated combination.
#'   * If \code{out = "all"}, a list containing an estimate, log-likelihood,
#'   maximum log-likelihood, and \eqn{{r}} values corresponding to the
#'   acceptance region determined by the significance level.
#'
#' @seealso \code{\link{ibdEstM}} for estimating the number of related pairs of
#'   strains and \code{\link{ibdDat}} for processing multi-sample data.
#' @export
#' @useDynLib dcifer

ibdPair <- function(pair, coi, afreq, nm, nr = 1e2, rval = NULL, reval = NULL,
                    logr = NULL, equalr = FALSE, out = "mle", alpha = 0.05,
                    freqlog = FALSE, neval = NULL, nloc = NULL, upcoi = TRUE) {
  if (is.null(nloc)) {
    nloc <- length(afreq)
  }
  if (is.null(logr)) {
    if (is.null(rval) && is.null(reval)) {
      rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
    }
    nrmat <- ifelse(equalr, 1, nm)
    if (is.null(reval)) {
      reval <- generateReval(nrmat, rval = rval)
    } else if (nrow(reval) != nrmat) {
      stop("reval doesn't match nm")
    }
    neval <- ncol(reval)
    logr <- logReval(reval, nm = nrmat, neval = neval)
  } else {                            # logr has to match reval
    if (is.null(reval)) {
      if (equalr && !is.null(rval)) {
        reval <- generateReval(1, rval = rval)
      } else {
        stop("Please provide reval")  # if equalr is TRUE, can provide rval
      }
    }
    if (is.null(neval)) {
      neval <- length(logr[[3]])
    }
  }

  if (!freqlog) {
    afreq <- lapply(afreq, log)
  }

  llik <- rep(0, neval)
  for (t in 1:nloc) {
    Ux <- which(as.logical(pair[[1]][[t]]))  # in case of integer vector
    Uy <- which(as.logical(pair[[2]][[t]]))

    if (length(Ux) == 0 ||                             # NA or all 0's
        length(Uy) == 0 ||                             # NA or all 0's
        any(afreq[[t]][unique(c(Ux, Uy))] == -Inf)) {  # likelihood = 0
      next
    }

    if (upcoi) {
      coix <- max(coi[1], length(Ux))
      coiy <- max(coi[2], length(Uy))
      llikt <- probUxUy(Ux, Uy, coix, coiy, afreq[[t]], logr, nm, neval, equalr)
    } else {
      Uxcomb <- getComb(Ux, coi[1])
      Uycomb <- getComb(Uy, coi[2])
      ncx <- ncol(Uxcomb)
      ncy <- ncol(Uycomb)
      llikt <- matrix(0, ncx*ncy, neval)
      icomb <- 1
      for (icx in 1:ncx) {
        for (icy in 1:ncy) {
          llikt[icomb, ] <- probUxUy(Uxcomb[, icx], Uycomb[, icy], coi[1],
                                     coi[2], afreq[[t]], logr, nm, neval,
                                     equalr)
          icomb <- icomb + 1
        }
      }
      llikt <- colMeans(llikt, na.rm = TRUE)  # or FALSE?
    }
    llik <- llik + llikt
  }

  if (tolower(out) == "llik") {
    return(llik)
  }
  imax <- which(llik == max(llik))
  est <- rowMeans(reval[, imax, drop = FALSE])  # or reval[, imax[1]] for 1st
  if (tolower(out == "mle")) {
    return(est)
  }

  qchi <- stats::qchisq(1 - alpha, df = ifelse(equalr, 1, nrow(reval)))
  cutoff  <- max(llik) - qchi/2
  itop    <- llik >= cutoff
  rtop <- reval[, itop, drop = equalr]
  return(list(mle = est, llik = llik, maxllik = llik[imax[1]], rtop = rtop))
}

#' Pairwise Genetic Distance
#'
#' @description Pairwise estimation of a relatedness parameter for polyclonal
#'   multiallelic samples.
#'
#' @details To be added
#'
#' @param dsmp a list containing sample genotypes.
#' @param rnull null value for the relatedness parameter.
#' @param ... additional arguments for \code{ibdPair}.
#' @inheritParams ibdPair
#'
#' @return A relatedness lower triangular matrix or 3-dimensional array, which
#'   in addition to MLE provides confidence intervals and a p-value for a
#'   specified null.
#'
#' @seealso \code{\link{ibdPair}} for genetic relatedness between two samples
#'   with an option of returning log-likelihood and \code{\link{ibdEstM}} for
#'   estimating the number of related pairs of strains.
#' @export

ibdDat <- function(dsmp, coi, afreq, nr = 1e2, rval = NULL, reval = NULL,
                   equalr = FALSE, out = "all", rnull = 0, alpha = 0.05, ...) {
  if (is.null(reval)) {
    if (is.null(rval)) {
      rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
    }
    reval <- generateReval(1, rval = rval)
  }
  logr  <- logReval(reval, nm = 1)
  neval <- length(rval)
  inull <- which.min(abs(rval - rnull))
  afreq <- lapply(afreq, log)
  nloc  <- length(afreq)
  nsmp  <- length(dsmp)

  if (out == "mle") {
    res <- matrix(NA, nsmp, nsmp, dimnames = list(names(dsmp), names(dsmp)))
  } else {
    res <- array(NA, dim = c(nsmp, nsmp, 4),
                 dimnames = list(names(dsmp), names(dsmp),
                                 c("MLE", "CI lower", "CI upper", "p-value")))
  }
  for (ix in 2:nsmp) {
    for (iy in 1:(ix - 1)) {
      rxy <- ibdPair(dsmp[c(ix, iy)], coi[c(ix, iy)], afreq, nm = 1,
                     reval = reval, logr = logr, out = out, alpha = alpha,
                     freqlog = TRUE, neval = neval, nloc = nloc)#, ...)
      if (out == "mle") {
        res[ix, iy] <- rxy
      } else {
        res[ix, iy, 1] <- rxy$mle
        res[ix, iy, 2:3] <- range(rxy$rtop)
        res[ix, iy, 4] <- 1 - stats::pchisq(2*(rxy$maxllik - rxy$llik[inull]),
                                            df = 1)
      }
    }
  }
  return(res)
}

#' Estimate Relatedness and Number of Related Strains
#'
#' @description Estimates multiple relatedness parameters when the number $M$ of
#'   related pairs of strains between two samples is unknown.
#'
#' @inheritParams ibdPair
#' @param nmmax maximum number of related pairs of strains to evaluate over.
#' @param nrs   an integer vector of values for grid resolution, one value for
#'   each \code{M} (or a single value when \code{equalr = TRUE}).
#' @param revals a list where each element is an \code{reval} matrix for a
#'   corresponding \code{M}.
#' @param logrs a list where each element is a \code{logr} list for a
#'   corresponding \code{M}.
#' @param tol tolerance value to decide if a parameter estimate is close enough
#'   to zero.
#' @return
#'   * If \code{out = "mle"}: a vector containing estimated \eqn{{r}}; the
#'   length of the vector is equal to estimated number of related pairs of
#'   strains between the two samples.
#'   * If \code{out = "llike"}, a vector of log-likelihood values for each
#'   evaluated combination.
#'   * If \code{out = "all"}, a list containing an estimate, log-likelihood,
#'   maximum log-likelihood, and \eqn{{r}} values corresponding to the
#'   acceptance region determined by the significance level.
#'
#' @seealso \code{\link{ibdPair}} for genetic relatedness between two samples
#'   with an option of returning log-likelihood and \code{\link{ibdDat}} for
#'   processing multi-sample data
#' @export

ibdEstM <- function(pair, coi, afreq, nmmax = 6,
                    nrs = c(1e3, 1e2, 32, 16, 12, 10), rval = NULL,
                    revals = NULL, logrs = NULL, equalr = FALSE,
                    out = "mle", alpha = 0.05, freqlog = FALSE, nloc = NULL,
                    upcoi = TRUE, tol = 1e-9) {
  nmmax <- min(coi, nmmax)
  if (is.null(nloc)) {
    nloc <- length(afreq)
  }
  if (!freqlog) {
    afreq <- lapply(afreq, log)
  }

  if (equalr) {
    if (is.null(revals)) {
      reval <- generateReval(1, rval = rval, nr = nrs[1])
    } else {
      reval <- revals[[1]]
    }
    resall  <- list()
    llikall <- numeric(nmmax)
    for (M in 1:nmmax) {
      res <- ibdPair(pair, coi, afreq, M, reval = reval, logr = logrs[[1]],
                     equalr = TRUE, out = "all", alpha = alpha, freqlog = TRUE,
                     nloc = nloc, upcoi = upcoi)
      resall[[M]] <- res
      llikall[M]  <- res$maxllik
    }
    imax   <- which.max(llikall)
    resmax <- resall[[imax]]
    res <- switch(tolower(out), "mle" = resmax["mle"], "llik" = resmax["llik"],
                  "all" = resmax)
    res$mle <- rep(resmax$mle, imax)
  } else {
    M <- 0
    mle <- 1
    while(all(mle > tol) && M < nmmax) {
      M <- M + 1
      if (is.null(revals)) {
        reval <- generateReval(M, nr = nrs[M])
      } else {
        reval <- revals[[M]]
      }
      res <- ibdPair(pair, coi, afreq, M, reval = reval, logr = logrs[[M]],
                     equalr = FALSE, out = out, alpha = alpha, freqlog = TRUE,
                     nloc = nloc, upcoi = upcoi)
      if (inherits(res, "list")) {
        mle <- res$mle
      } else {
        mle <- res
      }
    }
  }
  return(res)
}
