#' Calculate likelihood for a pair
#' \ifelse{html}{\out{U<sub>x</sub>}}{\eqn{U_x}},
#' \ifelse{html}{\out{U<sub>y</sub>}}{\eqn{U_y}}
#'
#' Calculates log-likelihood for a pair of samples at a single locus.
#'
#' @param Ux,Uy sets of unique alleles for two samples at a given locus. Vectors
#'   of indices corresponding to ordered probabilities in \code{probs}.
#' @param nx,ny complexity of infection for two samples. Vectors of length 1.
#' @param probs a vector of population allele frequencies (on a log scale) at a
#'   given locus. It is not checked if frequencies on a regular scale sum to
#'   \eqn{1}.
#' @param logr a list of length \eqn{5} as returned by \code{logReval()}.
#' @inheritParams ibdPair
#'
#' @return A vector of log-likelihood values for each \eqn{{r}}, i.e. for each
#'   column of \code{reval} matrix (or each element of \code{rval} vector if
#'   \code{equalr} is \code{TRUE}).
#'
#' @export

probUxUy <- function(Ux, Uy, nx, ny, probs, logr, nm, neval, equalr = FALSE) {
  if (nm > (min(nx, ny))) {
    stop("number of related strains greater than min(nx, ny)")
  }
  ixy <- which(Ux %in% Uy)
  iyx <- which(Uy %in% Ux)
  logj   <- log(1:max(nx, ny))         # starts with 1
  factj  <- lgamma(0:max(nx, ny) + 1)  # starts with 0
  if (equalr) {
    return(.Call("llikEqr", as.integer(Ux), as.integer(Uy),
                 as.integer(ixy), as.integer(iyx),
                 as.integer(nx),  as.integer(ny),
                 as.double(probs), as.double(logj), as.double(factj),
                 as.integer(nm), logr, as.integer(neval), PACKAGE = "dcifer"))
  }
  return(.Call("llik", as.integer(Ux), as.integer(Uy),
               as.integer(ixy), as.integer(iyx),
               as.integer(nx), as.integer(ny),
               as.double(probs), as.double(logj), as.double(factj),
               as.integer(nm), logr, as.integer(neval), PACKAGE = "dcifer"))
}

#' Calculate logs for reval and associated quantities
#'
#' @details Use \code{nm = 1} if \code{equalr = TRUE}.
#'
#' @inheritParams ibdPair
#' @return a list of length \eqn{5} that contains \eqn{log(r)}, \eqn{log(1 - r)}
#'   (matrices of the same dimensions as \code{reval}), the number of \eqn{r =
#'   1} for each column of \code{reval}, the number of \eqn{0 < r < 1} for each
#'   column, and the sum of \eqn{log(1 - r)} for each column.
#' @export

logReval <- function(reval, nm = NULL, neval = NULL) {
  if (is.null(nm))    nm    <- nrow(reval)
  if (is.null(neval)) neval <- ncol(reval)

  #  if (!is.loaded("src/logr")) dyn.load("src/logr.so")
  res <- .Call("logReval", as.double(reval), as.integer(neval), as.integer(nm))
  res[[1]] <- matrix(res[[1]], nm)
  res[[2]] <- matrix(res[[2]], nm)
  names(res) <- c("logr", "log1r", "m1", "nmid", "sum1r")
  return(res)
}

