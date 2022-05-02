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
#*** out if that option in ibdPair() is out
getComb <- function(Ux, coix) {
  if (length(Ux) <= coix) {
    return(matrix(Ux, ncol = 1))
  }
  return(utils::combn(Ux, coix))
}

drv1 <- function(r, C) sum(C/(C*r + 1))
#drv2 <- function(r, C) -sum((C/(C*r + 1))^2)
#llik <- function(r, p01) sum(log(r*(p01[2, ] - p01[1, ]) + p01[1, ]))

# If outside bounds, check 1st deriv at bound, then start there
mleNewton <- function(C, rstart = 0.5, tol = 1e-3, upper = 0.99) {
  rold <- 2
  rnew <- rstart
  counter <- 0
  while(abs(rnew - rold) > tol) {
    rold <- rnew
    sc <- drv1(rold, C)
    scdrv <- -sum((C/(C*rold + 1))^2)
    rnew <- rold - sc/scdrv
    counter <- counter + 1
    if (rnew < 0) {
      if (drv1(0, C) < 0) {
        return(c(0, counter))
      }
      rnew <- 0
    } else if (rnew >= 1) {
      if (drv1(1, C) > 0) {
        return(c(1, counter))
      }
      rnew <- upper
    }
  }
  return(c(rnew, counter))            #*** remove counter later
}

lrsP01 <- function(rhat, rnull, p01) {
  p1mp0 <- p01[2, ] - p01[1, ]
  lnull <- sum(log(rnull*p1mp0 + p01[1, ]))
  lhat  <- sum(log(rhat *p1mp0 + p01[1, ]))
  return(2*(lhat - lnull))
}
