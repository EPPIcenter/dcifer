
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
#' @param reval  a grid of \eqn{{r}} combinations, over which the likelihood
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
#'   or of length \code{nm} otherwise, containing the relatedness estimate.
#'   * If \code{out = "llik"}, a vector of log-likelihood values for each
#'   evaluated combination.
#'   * If \code{out = "pval"}, a list of the estimate and the p-value.
#'   * If \code{out = "all"}, a list containing an estimate, log-likelihood,
#'   maximum log-likelihood, and \eqn{{r}} values corresponding to the
#'   acceptance region determined by the significance level.
#'
#' @seealso \code{\link{ibdEstM}} for estimating the number of related pairs of
#'   strains and \code{\link{ibdDat}} for processing multi-sample data.
#' @export
#' @useDynLib dcifer

ibdPair <- function(pair, coi, afreq, nm, nr = 1e2, reval = NULL, logr = NULL,
                    equalr = FALSE, out = "mle", rnull = 0, inull = NULL,
                    alpha = 0.05, freqlog = FALSE, neval = NULL, nloc = NULL,
                    upcoi = TRUE, mnewton = TRUE, tol = 1e-3) {
                    #*** mnewton/tol is temp!!
  if (is.null(nloc)) {
    nloc <- length(afreq)
  }
  if (!freqlog) {
    afreq <- lapply(afreq, log)
  }
  logj   <- log(1:max(coi))            # starts with 1
  factj  <- lgamma(0:max(coi) + 1)     # starts with 0
  if (FALSE) {                         #*** uncomment for final version
  mnewton <- FALSE
  if (nm == 1 && out %in% c("mle", "pval")) {
    mnewton <- TRUE
  } }

  npar <- ifelse(equalr, 1, nm)
  if (!mnewton) {
    if (is.null(reval)) {
      reval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
    }
    if (!inherits(reval, "matrix")) {
      reval <- generateReval(npar, rval = reval)
    }
    if (!equalr && nm > 1 && nrow(reval) != nm) {
      stop("reval doesn't match nm")
    }
    if (is.null(neval)) {
      neval <- ncol(reval)
    }
    if (nm > 1) {
      if (is.null(logr)) {
        logReval(reval, nm = npar, neval = neval)
      }
    }
    llik <- rep(0, neval)
  } else {
    p01 <- matrix(0, 2, nloc)
  }

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
      rt <- probUxUy(Ux, Uy, coix, coiy, afreq[[t]], nm, logj, factj, logr,
                     reval, neval, mnewton, equalr)
    } else {  #*** remove later? Means for llikt or p01 not exactly justified
      Uxcomb <- getComb(Ux, coi[1])
      Uycomb <- getComb(Uy, coi[2])
      ncx <- ncol(Uxcomb)
      ncy <- ncol(Uycomb)
      rt <- matrix(0, ncx*ncy, neval)
      icomb <- 1
      for (icx in 1:ncx) {
        for (icy in 1:ncy) {
          rt[icomb, ] <- probUxUy(Uxcomb[, icx], Uycomb[, icy], coi[1], coi[2],
                                  afreq[[t]], nm, logj, factj, logr, reval,
                                  neval, mnewton, equalr)
          icomb <- icomb + 1
        }
      }
      rt <- colMeans(rt, na.rm = TRUE)
    }
    if (mnewton) {
      p01[, t] <- rt
    } else {
      llik <- llik + rt
    }
  }

  if (tolower(out) == "llik") {
    return(llik)
  }

  if (mnewton) {
    C <- p01[2, ]/p01[1, ] - 1
    rhat <- mleNewton(C, tol = tol)[1]  #*** for later make tol = 1/nr, rm [1]
  } else {
    imax <- which.max(llik)              # which(llik == max(llik))
    rhat <- reval[, imax]                # rowMeans(reval[, imax, drop = FALSE])
  }
  if (tolower(out) == "mle") {
    return(rhat)
  }

  if (mnewton) {
    lrs <- lrsP01(rhat, rnull, p01)
  } else {
    if (is.null(inull)) {
      if (equalr) {
        inull <- which.min(abs(reval - rnull))
      } else {      #*** needed? Probably won't be used much
        inull <- which.min(colSums(abs(reval - sort(rnull))))
      }
    }
    lrs <- 2*(llik[imax] - llik[inull])
  }
  adj  <- (rnull %in% c(0, 1)) + 1     # adjustment for one-sided test
  pval <- (1 - stats::pchisq(lrs, df = npar))/adj
  #*** BTW in the paper didn't address df for chisq when equalr = TRUE
  if (mnewton || tolower(out) == "pval") {
    return(list(mle = rhat, pval = pval))
  }

  #*** for ibdDat don't need this calculation to repeat, but not a big deal
  #*       (7.5 sec for 1e6)
  qchi   <- stats::qchisq(1 - alpha, df = npar)
  cutoff <- max(llik) - qchi/2
  itop   <- llik >= cutoff
  rtop   <- reval[, itop, drop = equalr]
  return(list(mle = rhat, pval = pval, llik = llik, maxllik = llik[imax],
              rtop = rtop))
  #*** note will have to get maxllik differently for mnewton = TRUE!!!
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

ibdDat <- function(dsmp, coi, afreq, nr = 1e3, reval = NULL, out = "all",
                   rnull = 0, alpha = 0.05, mnewton = TRUE, ...) {
  if (FALSE) {                     #*** uncomment later, rm mnewton from above
  mnewton <- TRUE
  if (tolower(out) == "all") {
    mnewton <- FALSE
  } }                               #*** end uncomment

  if (!mnewton) {
    if (is.null(reval)) {
      reval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
    }
    neval <- length(reval)
    inull <- which.min(abs(reval - rnull))
  } else {
    inull <- neval <- NULL
  }
  afreq <- lapply(afreq, log)
  nsmp  <- length(dsmp)
  nloc  <- length(afreq)

  if (tolower(out) == "mle") {
    res <- matrix(NA, nsmp, nsmp, dimnames = list(names(dsmp), names(dsmp)))
  } else if (tolower(out) == "pval") {
    res <- array(NA, dim = c(nsmp, nsmp, 2),
                 dimnames = list(names(dsmp), names(dsmp), c("MLE", "p-value")))
  } else {
    res <- array(NA, dim = c(nsmp, nsmp, 4),
                 dimnames = list(names(dsmp), names(dsmp),
                                 c("MLE", "CI lower", "CI upper", "p-value")))
  }
  for (ix in 2:nsmp) {
    for (iy in 1:(ix - 1)) {
      rxy <- ibdPair(dsmp[c(ix, iy)], coi[c(ix, iy)], afreq, nm = 1,
                     reval = reval, out = out, rnull = rnull, inull = inull,
                     alpha = alpha, freqlog = TRUE, neval = neval, nloc = nloc,
                     mnewton = mnewton, ...) #*** rm mnewton later???
      if (tolower(out) == "mle") {
        res[ix, iy] <- rxy
      } else {
        res[ix, iy, c("MLE", "p-value")] <- c(rxy$mle, rxy$pval)
        if (tolower(out) == "all" && !mnewton) {
          res[ix, iy, c("CI lower", "CI upper")] <- range(rxy$rtop)
        }
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
#' @param tol0 tolerance value to decide if a parameter estimate is close enough
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
                    nrs = c(1e3, 1e2, 32, 16, 12, 10), revals = NULL,
                    logrs = NULL, equalr = FALSE, out = "mle", rnull = 0,
                    alpha = 0.05, freqlog = FALSE, nloc = NULL, upcoi = TRUE,
                    tol0 = 1e-9, v2 = FALSE) {  #*** tol0 for rhat = 0, v2 temp
  nmmax <- min(coi, nmmax)
  if (is.null(nloc)) {
    nloc <- length(afreq)
  }
  if (!freqlog) {
    afreq <- lapply(afreq, log)
  }

  if (equalr) {
    if (is.null(revals)) {
      reval <- generateReval(1, nr = nrs[1])
    } else {
      reval <- revals[[1]]
    }
    neval <- length(reval)
    if (is.null(logrs)) {
      logrs <- list(logReval(reval, nm = 1, neval = neval))
    }
    if (!v2) {                  #*** temp: which method is faster
    resall  <- list()
    llikall <- numeric(nmmax)
    for (M in 1:nmmax) {
      res <- ibdPair(pair, coi, afreq, M, reval = reval, logr = logrs[[1]],
                     equalr = TRUE, out = "all", rnull = rnull, alpha = alpha,
                     freqlog = TRUE, neval = neval, nloc = nloc, upcoi = upcoi,
                     mnewton = FALSE) #*** mnewt out
      #*** might be faster to return llik, then recalculate at M' - check
      resall[[M]] <- res
      llikall[M]  <- res$maxllik
    }
    imax   <- which.max(llikall)
    res <- switch(tolower(out), "mle" = rep(resall[[imax]]$mle, imax),
                  "llik" = res[[imax]]$llik, "all" = resall[[imax]])
    } else {
      llikall <- numeric(nmmax)
      for (M in 1:nmmax) {
        res <- ibdPair(pair, coi, afreq, M, reval = reval, logr = logrs[[1]],
                       equalr = TRUE, out = "llik", rnull = rnull, alpha = alpha,
                       freqlog = TRUE, neval = neval, nloc = nloc, upcoi = upcoi,
                       mnewton = FALSE) #*** mnewt out
        llikall[M]  <- max(res)
      }
      imax   <- which.max(llikall)
      res <- ibdPair(pair, coi, afreq, imax, reval = reval, logr = logrs[[1]],
                     equalr = TRUE, out = out, rnull = rnull, alpha = alpha,
                     freqlog = TRUE, neval = neval, nloc = nloc, upcoi = upcoi,
                     mnewton = FALSE) #*** mnewt out
    }
  } else {
    M <- 0
    mle <- 1
    while(all(mle > tol0) && M < nmmax) {
      M <- M + 1
      if (is.null(revals)) {
        reval <- generateReval(M, nr = nrs[M])
      } else {
        reval <- revals[[M]]
      }
      if (is.null(logrs)) {
        logrs[[M]] <- logReval(reval, nm = M)
      }
      res <- ibdPair(pair, coi, afreq, M, reval = reval, logr = logrs[[M]],
                     equalr = FALSE, out = out, rnull = rnull, alpha = alpha,
                     freqlog = TRUE, nloc = nloc, upcoi = upcoi,
                     mnewton = FALSE) #*** mnewt out
      if (inherits(res, "list")) {
        mle <- res$mle
      } else {
        mle <- res
      }
    }
  }
  return(res)
}
