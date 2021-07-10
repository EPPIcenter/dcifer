
#' Genetic Distance for a Pair of Samples
#'
#' @description Calculates likelihood over a range of values for relatedness
#'   parameter \eqn{{r}} and provides an estimate for a pair of samples along
#'   with confidence regions based on a likelihood ratio test. For data with
#'   genotyping errors where estimated COI is less than |Ux|, various options
#'   are offered.
#' @param pair   a list of length two containing two samples.
#' @param afreq  a list of allele frequencies. Each element of the list
#'   corresponds to a locus.
#' @param coi    a vector indicating complexity of infection for each sample.
#' @param nr     an integer value for the resolution of the grid (\eqn{nr - 1}
#'   values between 0 and 1), over which the likelihood will be calculated.
#'   Ignored if non-null \code{reval} is provided.
#' @param nm     the number of related pairs of strains.
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
#'
#' @return
#'   * If \code{out = "mle"}: a vector of length 1 if \code{equalr = TRUE}
#'   or of length \code{nm} otherwise, containing estimated \eqn{{r}} (or vector
#'   /matrix if not just first).
#'   * If \code{out = "llike"}, a vector of log-likelihood values for each
#'   evaluated combination.
#'   * If \code{out = "all"}, a list containing an estimate, log-likelihood,
#'   maximum log-likelihood, \eqn{{r}} values corresponding to the acceptance
#'   region determined by the significance level, and the size of that region
#'   (as a proportion of all evaluated).
#'
#' @seealso \code{\link{ibdDat}} for processing multi-sample data in various
#'   formats and estimating pairwise relatedness.
#' @export
#' @useDynLib dcifer

ibdPair <- function(pair, afreq, coi, nr = 1e2, nm = min(coi), rval = NULL,
                    reval = NULL, logr = NULL, equalr = FALSE, out = "mle",
                    alpha = 0.05, freqlog = FALSE, neval = NULL, nloc = NULL) {
  if (is.null(nloc)) {
    nloc <- length(afreq)
  }
  if (is.null(nloc)) {
    nloc <- length(afreq)
  }
  if (is.null(logr)) {
    if (is.null(rval) && is.null(reval)) {
      rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
    }
    if (is.null(reval)) {
      reval <- generateReval(ifelse(equalr, 1, nm), rval = rval)
    }
    neval <- ncol(reval)
    logr <- logReval(reval, nm = ifelse(equalr, 1, nm), neval = neval)
  } else {
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

    if (length(Ux) == 0 ||                   # NA or all 0's
        length(Uy) == 0 ||                   # NA or all 0's
        any(afreq[[t]][unique(c(Ux, Uy))] == 0)) {  # likelihood = 0
      next                                   # likelihood = 0 or no data
    }
    coix <- max(coi[1], length(Ux))
    coiy <- max(coi[2], length(Uy))
    llikt <- probUxUy(Ux, Uy, coix, coiy, afreq[[t]], logr, nm, neval, equalr)
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
  proptop <- sum(itop)/neval
  rtop <- reval[, itop, drop = equalr]
  return(list(mle = est, llik = llik, maxllik = llik[imax[1]], rtop = rtop,
              proptop = proptop))
}

#' Genetic distance for data with genotyping errors
#'
#' @description Distance when estimated COI is less than |Ux|: likelihoods for
#'   all subsets of size = COI of Ux are calculated; geometric mean taken.
#'
#' @inherit ibdPair return params
# #' @export
#'

ibdPair2 <- function(pair, afreq, coi, nr = 1e2, nm = min(coi), rval = NULL,
                   reval = NULL, logr = NULL, equalr = FALSE, out = "mle") {
  nloc <- length(afreq)

  if (is.null(rval)) {
    rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
  }
  if (equalr) {
    neval <- length(rval)
  } else {
    if (is.null(reval)) {
      reval <- generateReval(nm, rval)
    }
    neval <- ncol(reval)
  }

  llik <- rep(0, neval)
  for (t in 1:nloc) {
    Ux <- which(as.logical(pair[[1]][[t]]))  # in case of integer vector
    Uy <- which(as.logical(pair[[2]][[t]]))  # sorted
    if (length(Ux) == 0 ||                   # NA or all 0's
        length(Uy) == 0 ||                   # NA or all 0's
        any(afreq[[t]][unique(c(Ux, Uy))] <= 0)) {  #*** < check 0 separately?
      next                                   # likelihood = 0 or no data
    }
    Uxcomb <- getComb(Ux, coi[1])
    Uycomb <- getComb(Uy, coi[2])
    ncx <- ncol(Uxcomb)
    ncy <- ncol(Uycomb)
    llikt <- matrix(0, ncx*ncy, neval)
    icomb <- 1
    for (icx in 1:ncx) {
      for (icy in 1:ncy) {
        llikt[icomb, ] <- probUxUy(Uxcomb[, icx], Uycomb[, icy], coi[1], coi[2],
                                   afreq[[t]], logr, nm, neval, equalr)
        icomb <- icomb + 1
      }
    }
    llikt <- colMeans(llikt, na.rm = TRUE)  # or FALSE?
    if (llikt[1] == -Inf) {
      if (all(llikt == -Inf)) {
        next
      } else {
        iinf <- which(llikt == -Inf)
        writeLines(paste("\n0 prob for r combinations",  # not added if 1st!!
                         paste(iinf, collapse = " ")))
      }
    } else {
      llik <- llik + llikt
    }
  }

  if (out == "llik") {
    return(llik)
  } else {
    imax <- which(llik == max(llik))
    if (equalr) {
      return(mean(rval[imax]))  # [1] first only; reval[imax] - all; mean()
    } else {
      return(rowMeans(reval[, imax, drop = FALSE]))
      # reval[, imax[1]]
    }
  }
}

#' Genetic Distance
#'
#' @description Pairwise estimation of relatedness parameters for polyclonal
#'   multiallelic samples.
#'
#' @details To be added: data formats that can be processed, formats to be
#'   returned.
#'
#' @param dat data containing sample genotypes. Various formats allowed.
#' @param nmmax the maximum number of related strains in a pair of samples.
#' @param split the allele separator character string if genotypes are provided
#'   as character strings.
#' @param FUNnm potentially a function to select \code{nm} for each pair of
#'   samples.
#' @param ...   additional arguments for FUNnm. #*** take out if FUNnm is out
#' @inheritParams ibdPair
#'
#' @return A lower triangular distance/relatedness matrix in a type of shaped
#'   list. Each pairwise distance contains a scalar, a vector, or a matrix of
#'   \eqn{{r}} estimates. Simple matrix (of type \code{double}) can be returned
#'   if all the elements are of length 1.
#'
#' @seealso \code{\link{ibdPair}} for genetic relatedness between two samples
#'   with an option of returning log-likelihood.
#' @export

ibdDat <- function(dat, afreq, coi, nmmax, nr = 1e2, rval = NULL, reval = NULL,
                   FUNnm = NULL, equalr = FALSE,
                   split = "[[:space:][:punct:]]+", ...) {
  if (inherits(dat, "matrix")) {
    nsmp <- nrow(dat)
    nloc <- ncol(dat)
    if (typeof(dat) == "list") {
      format <- "matlist"
    } else {
      format <- "matstr"
    }
  } else {
    nsmp <- length(dat)
    nloc <- length(dat[[1]])
    format <- "listlist"
  }

  if (is.null(rval)) {
    rval <- round(seq(0, 1, 1/nr), ceiling(log(nr, 10)))
  }
  if (equalr) {
    neval <- length(rval)
  } else {
    if (is.null(reval)) {
      reval <- generateRevalList(1:nmmax, rval)
      logr  <- mapply(logReval, reval, 1:nmmax)
    }
  }
  res  <- matrix(list(NA), nsmp, nsmp)

  for (ix in 1:nsmp) {
    for (iy in 1:(ix)) {
      nm <- min(coi[ix], coi[iy], nmmax)       # or FUNnm() - or AFTER ix == iy
      if (ix == iy) {                          # diagonal
        res[[ix, ix]] <- rep(1, ifelse(equalr, 1, min(coi[ix], nmmax)))
        next
      }
      llik <- rep(0, ifelse(equalr, neval, ncol(reval[[nm]])))
      for (t in 1:nloc) {
        if (format == "matstr") {
          alleles <- names(afreq[[t]])
          if (is.null(alleles))
            stop("Please provide names for allele frequencies")
          Ux  <- unlist(strsplit(dat[ix, t], split = split))
          Ux  <- sort(match(Ux, alleles))
          Uy  <- unlist(strsplit(dat[iy, t], split = split))
          Uy  <- sort(match(Uy, alleles))
        } else if (format == "matlist") {
          Ux <- which(as.logical(dat[[ix, t]]))
          Uy <- which(as.logical(dat[[iy, t]]))
        } else if (format == "listlist") {
          Ux <- which(as.logical(dat[[ix]][[t]]))
          Uy <- which(as.logical(dat[[iy]][[t]]))
        }
        if (length(Ux) == 0 || is.na(Ux[1]) ||   # NA or all 0's
            length(Uy) == 0 || is.na(Uy[1]) ||   # NA or all 0's
            any(afreq[[t]][unique(c(Ux, Uy))] <= 0)) {  #*** < check 0 separately?
          next                                   # likelihood = 0 or no data
        }
        llikt <- probUxUy(Ux, Uy, coi[ix], coi[iy], afreq[[t]], logr[[nm]],
                          nm, neval, equalr)
        if (any(llikt == -Inf)) {
          writeLines(paste("samples", ix, "and", iy, "- 0 prob for some r
                           combinations at locus", t))
          if (all(llikt == -Inf)) {       # -Inf is added otherwise
            next
          }
          llik <- llik + llikt            # lik = 0 if -Inf at any locus
        }
      }
      imax <- which(llik == max(llik))
      if (equalr) {
        res[[ix, iy]] <- mean(rval[imax])  # [1] first; rval[imax] all; mean()
      } else {
        res[[ix, iy]] <- rowMeans(reval[[nm]][, imax, drop = FALSE])
        # reval[[nm]][, imax[1]] first
      }
    }
  }
  if (all(sapply(res, length) == 1)) {
    res <- matrix(unlist(res), nrow(res), ncol(res))  # return simple matrix
  }
  return(res)
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
generateReval <- function(M, rval = NA, nr = NA) {
  if (is.na(rval[1])) {
    if (is.na(nr)) {
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

#' Generate a grid of parameter values for multiple values of \code{M}
#'
#' @param Ms an integer vector.
#' @param rvals a list of the length \code{max(Mv)}. Can also be a vector like
#'   \code{rval}; in that case the same vector will be used for all the values
#'   of \code{Mv}.
#' @param nrs an integer vector of the length \code{max(Mv)}.
#' @return A list of the length \code{max(Mv)} with elements corresponding to
#'   the values of \code{Mv}. Each element is a matrix with \code{M} rows and
#'   \code{nr + 1} or \code{length(rval)} columns.
#' @export
#' @rdname generateReval
#'
generateRevalList <- function(Ms, rvals = NA, nrs = NA) {
  Mmax <- max(Ms)
  revals <- as.list(rep(NA, Mmax))
  if (is.na(rvals[1]))          rvals <- rep(list(NA   ), Mmax)
  if (!inherits(rvals, "list")) rvals <- rep(list(rvals), Mmax)
  if (is.na(nrs[1]))            nrs   <- rep(NA, Mmax)
  for (M in Ms) {
    revals[[M]] <- generateReval(M, rvals[[M]], nrs[M])
  }
  return(revals)
}

# get combinations of alleles when length(Ux) > coix (used in ibdPair2)
getComb <- function(Ux, coix) {
  if (length(Ux) <= coix) {
    return(matrix(Ux, ncol = 1))
  }
  return(utils::combn(Ux, coix))
}



