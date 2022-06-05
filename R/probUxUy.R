#' Likelihood for
#' \ifelse{html}{\out{U<sub>x</sub>}}{\eqn{U_x}},
#' \ifelse{html}{\out{U<sub>y</sub>}}{\eqn{U_y}}
#'
#' Calculates log-likelihood for a pair of samples at a single locus.
#'
#' @param Ux,Uy sets of unique alleles for two samples at a given locus. Vectors
#'   of indices corresponding to ordered probabilities in \code{probs}.
#' @param nx,ny complexity of infection for two samples. Vectors of length 1.
#' @param probs a vector of population allele frequencies (on a log scale) at a
#'   given locus. It is not checked if frequencies on a regular scale sum to 1.
#' @param logj,factj numeric vectors containing precalculated logarithms and
#'   factorials.
#' @param mnewton a logical value. If \code{TRUE}, the coefficients for using
#'   Newton's method will be calculated.
#' @param logr a list of length 5 as returned by \code{\link{logReval}}.
#' @inheritParams ibdPair
#'
#' @return
#'   * If \code{mnewton = TRUE}, a vector of length 2 containing
#'   coefficients for fast likelihood calculation;
#'   * If \code{mnewton = FALSE}, a vector of length \code{neval} containing
#'   log-likelihoods for a range of parameter values.
#'
#' @examples
#' Ux <- c(1, 3, 7)                       # detected alleles at locus t
#' Uy <- c(2, 7)
#' coi <- c(5, 6)
#' aft <- runif(7)                        # allele frequencies for locus t
#' aft <- log(aft/sum(aft))
#'
#' logj  <- log(1:max(coi))
#' factj <- lgamma(0:max(coi) + 1)
#'
#' # M = 2, equalr = FALSE
#' M <- 2
#' reval <- generateReval(M, nr = 1e2)
#' logr  <- logReval(reval, M = M)
#' llikt <- probUxUy(Ux, Uy, coi[1], coi[2], aft, M, logj, factj,
#'                   equalr = FALSE, logr = logr, neval = ncol(reval))
#' length(llikt)
#'
#' # M = 2, equalr = TRUE
#' reval <- matrix(seq(0, 1, 0.001), 1)
#' logr  <- logReval(reval, M = 1)
#' llikt <- probUxUy(Ux, Uy, coi[1], coi[2], aft, M, logj, factj,
#'                   equalr = TRUE, logr = logr, neval = ncol(reval))
#'
#' # M = 1, mnewton = FALSE
#' M <- 1
#' llikt <- probUxUy(Ux, Uy, coi[1], coi[2], aft, M, logj, factj,
#'                   mnewton = FALSE, reval = seq(0, 1, 1e-3), neval = 1e3 + 1)
#'
#' # M = 1, mnewton = TRUE
#' probUxUy(Ux, Uy, coi[1], coi[2], aft, M, logj, factj, mnewton = TRUE)
#'
#' @export
probUxUy <- function(Ux, Uy, nx, ny, probs, M, logj, factj, equalr = FALSE,
                     mnewton = TRUE, reval = NULL, logr = NULL, neval = NULL) {
  ixy <- which(Ux %in% Uy)
  iyx <- which(Uy %in% Ux)
  if (M == 1) {
    if (mnewton) {
      return(.Call("p0p1", as.integer(Ux), as.integer(Uy), as.integer(ixy),
                   as.integer(iyx), as.integer(nx),  as.integer(ny),
                   as.double(probs), as.double(logj), as.double(factj),
                   PACKAGE = "dcifer"))
    } else {
      return(.Call("llikM1", as.integer(Ux), as.integer(Uy), as.integer(ixy),
                   as.integer(iyx), as.integer(nx),  as.integer(ny),
                   as.double(probs), as.double(logj), as.double(factj),
                   as.double(reval), as.integer(neval), PACKAGE = "dcifer"))
    }
  }
  if (equalr) {
    return(.Call("llikEqr", as.integer(Ux), as.integer(Uy), as.integer(ixy),
                 as.integer(iyx), as.integer(nx),  as.integer(ny),
                 as.double(probs), as.double(logj), as.double(factj),
                 as.integer(M), logr, as.integer(neval), PACKAGE = "dcifer"))
  }
  return(.Call("llik", as.integer(Ux), as.integer(Uy), as.integer(ixy),
               as.integer(iyx), as.integer(nx), as.integer(ny),
               as.double(probs), as.double(logj), as.double(factj),
               as.integer(M), logr, as.integer(neval), PACKAGE = "dcifer"))
}

#' Logarithms of \code{reval}
#'
#' Calculates logarithms of \code{reval} and \code{1 - reval}, as well as other
#' associated quantities.
#'
#' @details Use \code{M = 1} if \code{equalr = TRUE}.
#'
#' @inheritParams ibdPair
#' @return A list of length 5 that contains \code{log(reval)}, \code{log(1 -
#'   reval)}, the number of \code{reval = 1} for each column, the number of
#'   \code{0 < reval < 1} for each column, and \code{sum(log(1 - reval))} for
#'   each column.
#' @examples
#' reval <- generateReval(M = 2, nr = 1e2)
#' logr  <- logReval(reval, M = 2)
#' @export
#'
logReval <- function(reval, M = NULL, neval = NULL) {
  if (is.null(M))     M     <- nrow(reval)
  if (is.null(neval)) neval <- ncol(reval)

  res <- .Call("logReval", as.double(reval), as.integer(neval), as.integer(M))
  res[[1]] <- matrix(res[[1]], M)
  res[[2]] <- matrix(res[[2]], M)
  names(res) <- c("logr", "log1r", "m1", "nmid", "sum1r")
  return(res)
}

