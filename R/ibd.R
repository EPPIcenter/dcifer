
#' Relatedness Between Two Samples
#'
#' @description Provides estimates of relatedness between a pair of samples
#'   along with an optional support curve and inference.
#'
#' @details Handling of irregular cases:
#'   - Allele with population frequency of 0 is present: locus is skipped
#'     (does not contribute any information).
#'   - Number of unique alleles at a locus is greater than COI: COI will be
#'     increased for that locus only.
#'
#' @param pair   a list of length two containing data for a pair of samples.
#' @param coi    a vector containing complexity of infection for each sample.
#' @param afreq  a list of allele frequencies. Each element of the list
#'   corresponds to a locus.
#' @param M      the number of related pairs of strains.
#' @param rhat,pval,confreg,llik,maxllik  logical values specifying if
#'   relatedness estimate, p-value, confidence region, log-likelihood for a
#'   range of \eqn{r} values, and maximum log-likelihood should be returned.
#' @param rnull  a null value of relatedness parameter for hypothesis testing
#'   (needed if \code{pval = TRUE}).
#' @param alpha  significance level for a \ifelse{html}{\out{1 -
#'   &alpha;}}{\eqn{1 - \alpha}} confidence region.
#' @param equalr a logical value. If \code{TRUE}, the same level of relatedness
#'   is assumed for M pairs of strains
#'   (\ifelse{html}{\out{r<sub>1</sub>}}{\eqn{r_1}} = ... =
#'    \ifelse{html}{\out{r<sub>M</sub>}}{\eqn{r_M}}).
#' @param mnewton a logical value. If \code{TRUE}, Newton's method, adapted for
#'   a bounded parameter space, will be used to find MLE. If \code{mnewton} is
#'   not specified, it will be set to \code{TRUE} if \code{M = 1}, \code{confreg
#'   = FALSE}, and \code{llik = FALSE}.
#' @param freqlog a logical value indicating if \code{afreq} is on the log
#'   scale.
#' @param nr      an integer specifying precision of the estimate: resolution of
#'   a grid of parameter values (\ifelse{html}{\out{[0, 1]}}{\eqn{[0, 1]}}
#'   divided into \code{nr} equal intervals), over which the likelihood will be
#'   calculated. Ignored if non-null \code{reval} is provided.
#' @param reval  a matrix representing a grid of
#'   \ifelse{html}{\out{(r<sub>1</sub>, ..., r<sub>M</sub>)}}{\eqn{(r_1, ...,
#'   r_M)}} combinations, over which the likelihood will be calculated. Each
#'   column is a single combination.
#' @param tol    tolerance for calculating an estimate if \code{mnewton = TRUE}.
#'   Set to \code{1/nr} if not provided.
#' @param logr   a list as returned by \code{logReval} with logs of \code{reval}
#'   and other quantities.
#' @param neval  the number of relatedness values/combinations to evaluate over.
#' @param inull  an index of the value/column of \code{reval} that is closest to
#'   \code{rnull}.
#' @param nloc  the number of loci.
#'
#' @return A named list if multiple output logical values are \code{TRUE} - or a
#'   vector if only \code{rhat = TRUE} or \code{llik = TRUE}. Depending on
#'   these logical values, the following quantities are included:
#'   * If \code{rhat = TRUE}, a relatedness estimate (a vector of length 1 if
#'     \code{equalr = TRUE} or of length M if \code{equalr = FALSE});
#'   * If \code{pval = TRUE}, a p-value;
#'   * If \code{confreg = TRUE}, relatedness parameter values from the grid
#'     \code{reval} that are within \ifelse{html}{\out{1 - &alpha;}}{\eqn{1 -
#'     \alpha}} confidence region;
#'   * If \code{llik = TRUE}, log-likelihood values for relatedness parameter
#'     grid (provided in \code{reval} or determined by \code{nr});
#'   * If \code{maxllik = TRUE}, maximum log-likelihood.
#'
#' @examples
#' coi   <- getCOI(dsmp, lrank = 2)
#' afreq <- calcAfreq(dsmp, coi, tol = 1e-5)
#'
#' # two samples
#' ipair <- c(21, 17)
#' pair <- dsmp[ipair]
#' coip <-  coi[ipair]
#' M    <- 2
#'
#' res1 <- ibdPair(pair, coip, afreq, M = M, confreg = TRUE, alpha = 0.05,
#'                 equalr = FALSE, reval = revals[[M]])
#' res2 <- ibdPair(pair, coip, afreq, M = M, llik = TRUE,
#'                 equalr = TRUE, reval = revals[[1]])
#' res1$rhat
#' rep(res2$rhat, M)
#'
#' # plot confidence region
#' creg <- cbind(res1$confreg, res1$confreg[2:1, ])
#' plot(creg[1, ], creg[2, ], xlim = c(0, 1), ylim = c(0, 1), pch = 15,
#'      cex = 0.6, col = "cadetblue3", xlab = expression(hat(r)[1]),
#'      ylab = expression(hat(r)[2]))
#' points(res1$rhat, rev(res1$rhat), pch = 16)
#'
#' # plot log-likelihood
#' plot(revals[[1]], res2$llik, type = "l", xlab = "r", ylab = "log-likelihood")
#'
#' ipair <- c(41, 50)
#' pair <- dsmp[ipair]
#' coip <-  coi[ipair]
#'
#' # rtotal at different values of M with and without equality constraint
#' Mmax <- min(coip)
#' for (M in 1:Mmax) {
#'   print(paste0("M = ", M))
#'   print(c(sum(ibdPair(pair, coip, afreq, M = M, pval = FALSE,
#'                       equalr = FALSE, reval = revals[[M]])),
#'           ibdPair(pair, coip, afreq, M = M, pval = FALSE, equalr = TRUE)*M))
#'   cat("\n")
#' }
#'
#' # M = 1
#' # log-likelihood for specific r values
#' ibdPair(pair, coip, afreq, M = 1, rhat = FALSE, pval = FALSE, llik = TRUE,
#'         reval = c(0, 0.15, 0.38, 1))
#' # grid vs Newton's method
#' system.time(
#'   ibdPair(pair, coip, afreq, M = 1, mnewton = TRUE,  tol = 1e-5))
#' system.time(
#'   ibdPair(pair, coip, afreq, M = 1, mnewton = FALSE, nr  = 1e5))
#'
#' @seealso \code{\link{ibdEstM}} for estimating the number of related pairs of
#'   strains and \code{\link{ibdDat}} for processing multi-sample data.
#' @export
#' @useDynLib dcifer

ibdPair <- function(pair, coi, afreq, M, rhat = TRUE, pval = FALSE,
                    confreg = FALSE, llik = FALSE, maxllik = FALSE, rnull = 0,
                    alpha = 0.05, equalr = FALSE, mnewton = NULL,
                    freqlog = FALSE, nr = 1e3, reval = NULL, tol = NULL,
                    logr = NULL, neval = NULL, inull = NULL, nloc = NULL) {
  if (M > (min(coi))) {
    stop("number of related strains greater than min(coi)")
  }
  if (is.null(nloc)) {
    nloc <- length(afreq)
  }
  if (!freqlog) {
    afreq <- lapply(afreq, log)
  }
  maxj <- max(sapply(pair[[1]], sum), sapply(pair[[2]], sum), max(coi))
  logj   <- log(1:maxj)                # starts with 1
  factj  <- lgamma(0:maxj + 1)         # starts with 0
  if (is.null(mnewton)) {
    mnewton <- FALSE
    if (M == 1 && !confreg && !llik) {
      mnewton <- TRUE
      if (is.null(tol)) {
        tol <- 1/nr
      }
    }
  }

  npar <- if (equalr) 1 else M
  if (!mnewton) {
    if (!inherits(reval, "matrix")) {
      reval <- generateReval(npar, rval = reval, nr = nr)
    } else if (nrow(reval) != npar) {
      stop("reval doesn't match M")
    }
    if (is.null(neval)) {
      neval <- ncol(reval)
    }
    if (is.null(logr)) {
      logr <- logReval(reval, M = M, neval = neval, equalr = equalr)
    }
    llikr <- rep(0, neval)
  } else {
    p01 <- matrix(NA, 2, nloc)
  }

  for (t in 1:nloc) {
    Ux <- which(as.logical(pair[[1]][[t]]))  # in case of integer vector
    Uy <- which(as.logical(pair[[2]][[t]]))
    if (length(Ux) == 0 ||                             # NA or all 0's
        length(Uy) == 0 ||                             # NA or all 0's
        any(afreq[[t]][unique(c(Ux, Uy))] == -Inf)) {  # likelihood = 0
      next
    }
    coix <- max(coi[1], length(Ux))
    coiy <- max(coi[2], length(Uy))
    rt <- probUxUy(Ux, Uy, coix, coiy, afreq[[t]], M, logj, factj, equalr,
                   mnewton, reval, logr, neval)
    if (mnewton) {
      p01[, t] <- rt
    } else {
      llikr <- llikr + rt
    }
  }

  if (llik && !rhat && !pval && !confreg && !maxllik) {
    return(llikr)
  }
  if (mnewton) {
    p01 <- p01[, !is.na(p01[1, ])]
    C   <- p01[2, ]/p01[1, ] - 1
    rhat <- rNewton(C, tol = tol, off = tol)
  } else {
    imax <- which.max(llikr)           # which(llikr == max(llikr))
    rhat <- reval[, imax]              # rowMeans(reval[, imax, drop = FALSE])
  }
  if (!pval && !confreg && !llik && !maxllik) {
    return(rhat)
  }

  res <- list(rhat = rhat)
  if (pval) {
    if (mnewton) {
      lrs <- lrsP01(rhat, rnull, p01)
    } else {
      if (is.null(inull)) {
        if (equalr) {
          inull <- which.min(abs(reval - rnull))
        } else {
          inull <- which.min(colSums(abs(reval - sort(rnull))))
        }
      }
      lrs <- 2*(llikr[imax] - llikr[inull])
    }
    adj  <- (rnull %in% c(0, 1)) + 1   # adjustment for one-sided test
    res$pval <- (1 - stats::pchisq(lrs, df = npar))/adj
  }

  if (confreg) {
    qchi   <- stats::qchisq(1 - alpha, df = npar)
    cutoff <- max(llikr) - qchi/2
    itop   <- llikr >= cutoff
    res$confreg <- reval[, itop, drop = equalr]
  }

  if (llik) {
    res$llik <- llikr
  }
  if (maxllik) {                       # for ibdEstM()
    if (mnewton) {
      C <- p01[2, ]/p01[1, ] - 1
      rhat <- rNewton(C, tol = tol, off = tol)
      res$maxllik <- sum(log(rhat*(p01[2, ] - p01[1, ]) + p01[1, ]))
    } else {
      res$maxllik <- max(llikr)
    }
  }
  return(res)
}

#' Pairwise Relatedness
#'
#' @description Provides pairwise relatedness estimates within a dataset or
#'   between two datasets along with optional p-values and confidence intervals
#'   (CI).
#'
#' @details For this function, \ifelse{html}{\out{M}}{\eqn{M}} is set to 1. If
#'   \code{confint = FALSE}, Newton's method is used to find the estimates,
#'   otherwise the likelihood is calculated for a grid of parameter values.
#'
#' @param dsmp     a list with each element corresponding to one sample.
#' @param dsmp2    a list representing a second dataset.
#' @param coi2     a vector with complexities of infection for a second dataset.
#' @param pval     a logical value specifying if p-values should be returned.
#' @param confint  a logical value specifying if confidence intervals should be
#'   returned.
#' @param reval    a vector or a single-row matrix. A grid of parameter values,
#'   over which the likelihood will be calculated.
#' @inheritParams ibdPair
#'
#' @return A matrix if \code{pval} and \code{confint} are \code{FALSE} and
#'   3-dimensional arrays otherwise. The matrices are lower triangular if
#'   distances are calculated within a dataset. For a 3-dimensional array,
#'   stacked matrices contain relatedness estimates, p-values, and endpoints of
#'   confidence intervals (if requested).
#'
#' @examples
#' coi   <- getCOI(dsmp, lrank = 2)           # estimate COI
#' afreq <- calcAfreq(dsmp, coi, tol = 1e-5)  # estimate allele frequencies
#'
#' # subset of samples for faster processing
#' i1 <- 1:15     # from Maputo
#' i2 <- 31:40    # from Inhambane
#' isub <- c(i1, i2)
#'
#' # matrix is returned
#' dres1 <- ibdDat(dsmp[isub], coi[isub], afreq, pval = FALSE)
#' dim(dres1)
#'
#' # test a null hypothesis H0: r = 0, change precision
#' dres2 <- ibdDat(dsmp[isub], coi[isub], afreq, pval = TRUE, rnull = 0,
#'                 nr = 1e2)
#' dim(dres2)
#'
#' # test H0: r = 0.2, include 99% confidence intervals
#' dres3 <- ibdDat(dsmp[isub], coi[isub], afreq, pval = TRUE, confint = TRUE,
#'                 rnull = 0.2, alpha = 0.01)
#' dres3[2, 1, ]
#'
#' # pairwise relatedness between two datasets, H0: r = 0
#' drbetween <- ibdDat(dsmp[i1], coi[i1], afreq,
#'                     dsmp2 = dsmp[i2], coi2 = coi[i2])
#' dim(drbetween)
#' drbetween[1, 2, ]
#' sum(is.na(drbetween[, , 1]))
#'
#' @seealso \code{\link{ibdPair}} for genetic relatedness between two samples
#'   and \code{\link{ibdEstM}} for estimating the number of related pairs of
#'   strains.
#' @export

ibdDat <- function(dsmp, coi, afreq, dsmp2 = NULL, coi2 = NULL, pval = TRUE,
                   confint = FALSE, rnull = 0, alpha = 0.05, nr = 1e3,
                   reval = NULL) {
  dwithin <- is.null(dsmp2)
  if (confint) {
    mnewton <- FALSE
    tol     <- NULL
  } else {
    mnewton <- TRUE
    tol     <- 1/nr
  }

  if (!mnewton) {
    if (!inherits(reval, "matrix")) {
      reval <- generateReval(M = 1, rval = reval, nr = nr)
    }
    neval <- ncol(reval)
    logr  <- logReval(reval, M = 1, neval = neval)
  } else {
    neval <- NULL
  }
  #***
#  neval  <- if (mnewton)          NULL else ncol(reval)
  inull  <- if (mnewton || !pval) NULL else which.min(abs(reval - rnull))
  afreq  <- lapply(afreq, log)
  nloc   <- length(afreq)
  nsmp   <- length(dsmp)
  snames <- names(dsmp)
  if (dwithin) {
    dsmp2 <- dsmp
    coi2  <- coi
  }
  nsmp2   <- length(dsmp2)
  snames2 <- names(dsmp2)

  dim3names <- c("estimate", "p_value", "CI_lower", "CI_upper")
  if (!pval && !confint) {
    res <- matrix(NA, nsmp, nsmp2, dimnames = list(snames, snames2))
  } else if (pval && !confint) {
    res <- array(NA, dim = c(nsmp, nsmp2, 2),
                 dimnames = list(snames, snames2, dim3names[c(1:2)]))
  } else if (!pval && confint) {
    res <- array(NA, dim = c(nsmp, nsmp2, 3),
                 dimnames = list(snames, snames2, dim3names[c(1, 3:4)]))
  } else {
    res <- array(NA, dim = c(nsmp, nsmp2, 4),
                 dimnames = list(snames, snames2, dim3names))
  }

  xstart <- if (dwithin) 2 else 1
  for (ix in xstart:nsmp) {
    yend  <- if (dwithin) ix - 1 else nsmp2
    for (iy in 1:yend) {
      rxy <- ibdPair(list(dsmp[[ix]], dsmp2[[iy]]), c(coi[ix], coi2[iy]), afreq,
                     M = 1, pval = pval, confreg = confint, rnull = rnull,
                     alpha = alpha, mnewton = mnewton, freqlog = TRUE,
                     reval = reval, tol = tol, logr = logr, neval = neval,
                     inull = inull, nloc = nloc)
      if (!pval && !confint) {
        res[ix, iy] <- rxy
      } else {
        res[ix, iy, "estimate"] <- rxy$rhat
      }
      if (pval) {
        res[ix, iy, "p_value"]  <- rxy$pval
      }
      if (confint) {
        res[ix, iy, c("CI_lower", "CI_upper")] <- range(rxy$confreg)
      }
    }
  }
  return(res)
}

#' Estimate Relatedness and a Number of Related Strains
#'
#' @description Estimates the number of related pairs of strains between two
#'   infections along with corresponding relatedness estimates and optional
#'   inference.
#'
#' @inheritParams ibdPair
#' @param Mmax  a maximum number of related pairs of strains to evaluate over.
#'   If greater than \code{min(coi)}, will be set to \code{min(coi)}.
#' @param pval,confreg,llik logical values specifying if p-value, confidence
#'   region, and log-likelihood for a range of \eqn{r} values should be
#'   returned.
#' @param nrs   an integer vector where \code{i}'th element correspons to
#'   \ifelse{html}{\out{M = i}}{\eqn{M = i}} and indicates precision of the
#'   estimate (resolution of a grid of parameter values). Ignored if non-null
#'   \code{revals} is provided.
#' @param revals a list where \code{i}'th element corresponds to
#'   \ifelse{html}{\out{M = i}}{\eqn{M = i}} and is a matrix representing a grid
#'   of parameter values (a matrix where each column represents a single
#'   \ifelse{html}{\out{(r<sub>1</sub>, ..., r<sub>M</sub>)}}{\eqn{(r_1, ...,
#'   r_M)}} combination).
#' @param tol0   a tolerance value for an estimate to be considered zero.
#' @param logrs  a list where \code{i}'th element corresponds to
#'   \ifelse{html}{\out{M = i}}{\eqn{M = i}} and is a list as returned by
#'   \code{\link{logReval}}.
#' @param nevals a vector where \code{i}'th element corresponds to
#'   \ifelse{html}{\out{M = i}}{\eqn{M = i}} and provides the number of
#'   relatedness values/combinations to evaluate over.
#'
#' @return A named list if multiple output logical values are \code{TRUE} - or a
#'   vector if only \code{rhat = TRUE}. The output includes:
#'   * a relatedness estimate (numeric vector of length corresponding to the
#'     estimated number of related pairs);
#'   * a p-value if \code{pval = TRUE};
#'   * parameter values from the grid in \code{revals} that are within the
#'     confidence region if \code{confreg = TRUE};
#'   * log-likelihood values for the parameter grid in \code{revals} if
#'     \code{llik = TRUE}.
#'
#' @examples
#' coi   <- getCOI(dsmp, lrank = 2)           # estimate COI
#' afreq <- calcAfreq(dsmp, coi, tol = 1e-5)  # estimate allele frequencies
#'
#' # two samples
#' ipair <- c(21, 17)
#' # for higher COI: c(33, 5): COI = 5-6; c(37, 20): 4-3, c(41, 50): 5-4
#'
#' Mmax  <- min(coi[ipair])
#' # choose resolution of the grid for different M
#' nrs   <- c(1e3, 1e2, 32, 16, 12, 10)[1:Mmax]
#' revals <- mapply(generateReval, 1:Mmax, nr = nrs)
#'
#' (res1 <- ibdEstM(dsmp[ipair], coi[ipair], afreq, Mmax = Mmax, equalr = FALSE,
#'                  reval = revals))
#' (res2 <- ibdEstM(dsmp[ipair], coi[ipair], afreq, Mmax = Mmax, equalr = TRUE))
#' # number of related pairs of strains (M')
#' sum(res1 > 0)
#' sum(res2 > 0)  # can be 0's
#'
#' @seealso \code{\link{ibdPair}} for estimates of relatedness between two
#'   samples and \code{\link{ibdDat}} for pairwise relatedness estimates within
#'   a dataset or between two datasets.
#' @export
#'

ibdEstM <- function(pair, coi, afreq, Mmax = 6, pval = FALSE, confreg = FALSE,
                    llik = FALSE, rnull = 0, alpha = 0.05, equalr = FALSE,
                    freqlog = FALSE, nrs = c(1e3, 1e2, 32, 16, 12, 10),
                    revals = NULL, tol0 = 1e-9, logrs = NULL, nevals = NULL,
                    nloc = NULL) {
  Mmax <- min(coi, Mmax)
  if (is.null(nloc)) {
    nloc <- length(afreq)
  }
  if (!freqlog) {
    afreq <- lapply(afreq, log)
  }

  if (equalr) {
    reval <- if (is.null(revals)) generateReval(1, nr = nrs[1]) else revals[[1]]
    neval <- length(reval)
    logr  <- if (is.null(logrs))
      logReval(reval, M = 1, neval = neval) else logrs[[1]]
    sum1r <- logr$sum1r
    inull <- if (pval) which.min(abs(reval - rnull)) else NULL

    resall  <- list()
    llikall <- numeric(Mmax)
    for (M in 1:Mmax) {
      logr$sum1r <- sum1r*M
      res <- ibdPair(pair, coi, afreq, M, pval = pval, confreg = confreg,
                     llik = llik, maxllik = TRUE, rnull = rnull, alpha = alpha,
                     equalr = TRUE, freqlog = TRUE, nr = nrs[1], reval = reval,
                     logr = logr, neval = neval, inull = inull, nloc = nloc)
      resall[[M]] <- res
      llikall[M]  <- res$maxllik
    }
    imax <- which.max(llikall)
    res  <- resall[[imax]]
    res$rhat    <- rep(res$rhat, imax)
    res$maxllik <- NULL
    if (length(res) == 1) {
      res <- res[[1]]
    }
  } else {
    M <- 0
    rhat <- 1
    while(all(rhat > tol0) && M < Mmax) {
      M <- M + 1
      res <- ibdPair(pair, coi, afreq, M, pval = pval, confreg = confreg,
                     llik = llik, maxllik = FALSE, rnull = rnull, alpha = alpha,
                     equalr = FALSE, freqlog = TRUE, nr = nrs[M],
                     reval = revals[[M]], logr = logrs[[M]], neval = nevals[M],
                     nloc = nloc)
      if (inherits(res, "list")) {
        rhat <- res$rhat
      } else {
        rhat <- res
      }
    }
  }
  return(res)
}

# If outside bounds, check 1st deriv at bound, then start next iter at bound
rNewton <- function(C, rstart = 0.5, tol = 1e-3, toldrv = 1e-3, off = 1e-3) {
  rold <- 2
  sc   <- 100
  rnew <- rstart
  while(abs(rnew - rold) > tol || abs(sc) > toldrv) {
    rold <- rnew
    sc <- drv1(rold, C)
    scdrv <- -sum((C/(C*rold + 1))^2)  # 2nd derivative of llik
    (rnew <- rold - sc/scdrv)
    if (rnew < 0) {
      if (drv1(0, C) < 0) {
        return(0)
      }
      rnew <- 0
    } else if (rnew >= 1) {
      if (drv1(1, C) > 0) {
        return(1)
      }
      rnew <- 1 - off
    }
  }
  return(rnew)
}

drv1 <- function(r, C) sum(C/(C*r + 1))

lrsP01 <- function(rhat, rnull, p01) {
  p1mp0 <- p01[2, ] - p01[1, ]
  lnull <- sum(log(rnull*p1mp0 + p01[1, ]))
  lhat  <- sum(log(rhat *p1mp0 + p01[1, ]))
  return(2*(lhat - lnull))
}
