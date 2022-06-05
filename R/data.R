#' Sample data
#'
#' Microhaplotype data from two clinics in Mozambique. Samples are sorted by location.
#'
#' @format A list of 52 elements (samples), each element is a list of 87
#'   elements (loci), which are integer vectors (alleles).
"dsmp"

#' Dcifer results
#'
#' Results of relatedness estimation.
#'
#' @format A three-dimensional array with 52 columns, 52 rows, and 4 matrices.
#'   Dimension names correspond to sample ID's (rows and columns) and types of
#'   results \code{c("estimate", "p_value", "CI_lower", "CI_upper")} (matrices).
"dres"

#' Parameter grid
#'
#' Precalculated parameter grids for a range of values of M (from 1 to 5).
#'
#' @format A list of length 5, where each element corresponds to a single value
#'   of M and is a matrix with M rows. Each column of a matrix is a
#'   \ifelse{html}{\out{r<sub>1</sub>}}{\eqn{r_1}} = ... =
#'   \ifelse{html}{\out{r<sub>M</sub>}}{\eqn{r_M}}) combination.
"revals"


