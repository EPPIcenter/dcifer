
#' Plot Relatedness Estimates
#'
#' Represents a matrix of pairwise relatedness estimates with colors
#' corresponding to the levels of relatedness. Optionally, also outlines results
#' of a hypothesis testing. The plot follows a matrix layout.
#'
#' @param r a matrix or a 3-dimensional array as returned by
#'   \code{\link{ibdDat}}.
#' @param rlim the range of values for colors. If \code{NULL} or \code{NA}, will
#'   be calculated from \code{r}.
#' @param sig a logical matrix specifying which entries of the relatedness
#'   matrix should be outlined or "amplified" (made larger). \code{sig} takes
#'   precedence over \code{alpha}.
#' @param alpha significance level for hypothesis testing; determines
#'   relatedness matrix entries to be outlined. Ignored if \code{sig} is not
#'   \code{NULL}.
#' @param col the colors for the range of relatedness values.
#' @param draw_diag a logical value specifying if diagonal cells should be
#'   distinguished from others by a separate color.
#' @param col_diag,border_diag,lwd_diag the color for the fill, the color for
#'   the border, and the line width for the border of diagonal entries. Ignored
#'   if \code{draw_diag = FALSE}.
#' @param border_sig,lwd_sig the color and the line width for outlining entries
#'   specified by \code{sig} or \code{alpha}. If \code{border_sig} is \code{NA},
#'   these entries will be "amplified", their size controlled by \code{lwd_sig}.
#' @param xlab,ylab axis labels.
#' @param add a logical value specifying if the graphics should be added to the
#'   existing plot (useful for triangular matrices).
#' @param idlab a logical value specifying if sample ID's should be displayed.
#' @param side_id an integer vector specifying plot sides for sample ID labels.
#' @param col_id,cex_id numeric vectors for the color and the size of sample ID
#'   labels.
#' @param srt_id a vector of the same length as \code{side_id} specifying
#'   rotation angles for sample ID labels. If \code{NULL}, the labels will be
#'   perpendicular to the axes.
#' @param ... other graphical parameters.
#' @return \code{NULL}; called for plotting.
#'
#' @examples
#' parstart <- par(no.readonly = TRUE)   # save starting graphical parameters
#'
#' par(mar = c(0.5, 0.5, 0.5, 0.5))
#' plotRel(dres, alpha = 0.05, draw_diag = TRUE)
#'
#' # draw log of p-values in the upper triangle
#' pmat <- t(log(dres[, , "p_value"]))
#' pmat[pmat == -Inf] <- min(pmat[is.finite(pmat)])
#' plotRel(pmat, rlim = NULL, draw_diag = TRUE, col = hcl.colors(101, "PuRd"),
#'         add = TRUE, col_diag = "slategray2", border_diag = 1)
#'
#' # symmetric matrix, outline significant in upper triangle, display sample ID
#' par(mar = c(3, 3, 0.5, 0.5))
#' dmat  <- dres[, , "estimate"]
#' dmats <- mixMat(dmat, dmat)
#' sig <- dres[, , "p_value"] <= 0.05
#' col_id <- rep(c("orchid4", "cadetblue4"), each = 26)
#' plotRel(dmats, sig = t(sig), border_sig = "magenta2", draw_diag = TRUE,
#'         idlab = TRUE, col_id = col_id)
#' abline(v = 26, h = 26, col = "gray45", lty = 5)
#'
#' # rotate sample ID labels on all sides, increase size for significant pairs
#' par(mar = c(3, 3, 3, 3))
#' sig <- dres[, , "p_value"] <= 0.01
#' plotRel(dmats, sig = mixMat(sig, sig), border_sig = NA, lwd_sig = 5,
#'         draw_diag = TRUE, idlab = TRUE, side_id = 1:4, col_id = col_id,
#'         srt_id = c(-55, 25, 65, -35))
#' par(parstart)
#'
#' @seealso \code{\link{plotColorbar}} for a colorbar and \code{\link{mixMat}}
#'   for combining square matrces.
#' @export
#'
plotRel <- function(r, rlim = c(0, 1), sig = NULL, alpha = NULL,
                    col = grDevices::hcl.colors(101, "YlGnBu", rev = TRUE),
                    draw_diag = FALSE, col_diag = "gray", border_diag = NA,
                    lwd_diag = 0.5, border_sig = "orangered2", lwd_sig = 1.5,
                    xlab = "", ylab = "", add = FALSE, idlab = FALSE,
                    side_id = c(1, 2), col_id = 1, cex_id = 0.5,
                    srt_id = NULL, ...) {
  z <- if (length(dim(r)) == 2) r else r[, , "estimate"]
  nx <- ncol(z)
  ny <- nrow(z)
  nxy  <- min(nx, ny)
  if (is.null(rlim) || is.na(rlim[1])) {
    rlim <- range(z, na.rm = TRUE)
  }

  # scale the matrix
  z <- round((z - rlim[1])/diff(rlim)*(length(col) - 1))

  if (!add) {
    mrg <- 0.01
    plot(NULL, xlim = c(0, nx) + c(-1, 1)*mrg*nx,
         ylim = rev(c(0, ny) + c(-1, 1)*mrg*ny),
         xlab = xlab, ylab = ylab, axes = FALSE, xaxs = "i", yaxs = "i", ...)
  }
  for (j in 1:nx) {
    for (i in 1:ny) {
      graphics::polygon(c(j - 1, j)[c(1, 2, 2, 1, 1)],
                        c(i - 1, i)[c(1, 1, 2, 2, 1)],
                        col = col[z[i, j] + 1], border = NA, ljoin = 1,
                        lend = 2)
    }
  }

  # outlined or amplified (order randomized)
  if (is.null(sig) && !is.null(alpha)) {
    sig <- r[, , "p_value"] <= alpha
  }
  if (!is.null(sig)) {
    jrand <- sample(1:nx)
    for (j in jrand) {
      irand <- sample(1:ny)
      for (i in irand) {
        colij <- if (is.na(border_sig)) col[z[i, j] + 1] else border_sig
        if (!is.na(sig[i, j]) && sig[i, j]) {
          graphics::lines(c(j - 1, j)[c(1, 2, 2, 1, 1)],
                          c(i - 1, i)[c(1, 1, 2, 2, 1)],
                          col = colij, lwd = lwd_sig, ljoin = 1, lend = 2)
        }
      }
    }
  }

  if (draw_diag) {                            # to overlay borders
    for (i in 1:nxy) {
      graphics::polygon(c(i - 1, i)[c(1, 2, 2, 1, 1)],  # in case transparent
                        c(i - 1, i)[c(1, 1, 2, 2, 1)], col = "white",
                        border = ifelse(is.na(border_diag), NA, "white"),
                        lwd = lwd_diag, ljoin = 1, lend = 2)
      graphics::polygon(c(i - 1, i)[c(1, 2, 2, 1, 1)],
                        c(i - 1, i)[c(1, 1, 2, 2, 1)],
                        col = col_diag, border = border_diag, lwd = lwd_diag,
                        ljoin = 1, lend = 2)
    }
  }

  if (idlab) {
    cc <- graphics::par("usr")
    cc <- (cc + c(c(-1, 1)*diff(cc[1:2]),
                  c(-1, 1)*diff(cc[3:4]))*0.00)[c(3, 1, 4, 2)]
    for (i in 1:length(side_id)) {
      if (side_id[i] %in% c(1, 3)) {
        labs <- colnames(z)
        at <- xc <- (1:nx) - 0.5
        yc <- cc[side_id[i]]
      } else {
        labs <- rownames(z)
        at <- yc <- (1:ny) - 0.5
        xc <- cc[side_id[i]]
      }
      if (is.null(srt_id) || (side_id[i] %in% c(1, 3) && srt_id[i] == 90) ||
          (side_id[i] %in% c(2, 4) && srt_id[i] == 0)) {
        graphics::mtext(labs, side_id[i], line = 0, at = at, col = col_id,
                        cex = cex_id, las = 2)
      } else {
        graphics::text(xc, yc, labs, adj = adjID(side_id[i], srt_id[i]),
                       col = col_id, cex = cex_id, srt = srt_id[i], xpd = TRUE)
      }
    }
  }
}
adjID <- function(s, srt) {
  srt <- srt %% 360
  if (s == 4 || (s == 1 && srt > 180) || (s == 3 && srt < 180)) 0 else 1
}

#' Colorbar
#'
#' Creates a colorbar for a plot.
#'
#' @param rlim  the range of values to be represented by colors.
#' @param by    increment size for tickmark locations. Ignored if \code{at} is provided.
#' @param at    a vector of tickmark locations.
#' @param horiz a logical value specifying if the colorbar should be drawn horizontally.
#' @param col   the colors for the colorbar.
#' @return \code{NULL}; called for plotting.
#' @inheritParams plotRel
#' @details The colorbar will fill the whole plotting region, which needs to be
#'   specified outside of this function to control proportions and location of
#'   the colorbar (see examples). To match the colors in the main plot,
#'   \code{rlim} values should be the same for \code{plotRel} and
#'   \code{plotColorbar}; if \code{rlim = NULL} or \code{rlim = NA} in
#'   \code{plotRel}, provide the actual range of relatedness estimates for
#'   \code{plotColorbar} (see examples).
#'
#' @examples
#' parstart <- par(no.readonly = TRUE)   # save starting graphical parameters
#'
#' # colorbar on the side of the main plot
#' layout(matrix(1:2, 1), width = c(7, 1))
#' par(mar = c(2, 0, 2, 0) + 0.1)
#' # make symmetric matrix
#' dmat <- dres[, , "estimate"]
#' dmat[upper.tri(dmat)] <- t(dmat)[upper.tri(t(dmat))]
#' isig <- which(dres[, , "p_value"] <= 0.05, arr.ind = TRUE)
#' plotRel(dmat, draw_diag = TRUE, isig = rbind(isig, isig[, 2:1]))
#' abline(v = 26, h = 26, col = "gray45", lty = 5)
#' par(mar = c(2, 1, 2, 2) + 0.1)
#' plotColorbar()
#'
#' # shorter colorbar, tick mark locations provided
#' par(mar = c(2, 0, 2, 0) + 0.1)
#' plotRel(dmat, draw_diag = TRUE, isig = rbind(isig, isig[, 2:1]))
#' par(mar = c(5, 0.5, 5, 2.5) + 0.1)
#' plotColorbar(at = c(0.0625, 0.125, 0.25, 0.5, 0.78))
#' par(parstart)
#'
#' # triangular matrix, inset horizontal colorbar
#' par(mar = c(1, 1, 1, 1))
#' plotRel(dres, rlim = NULL, draw_diag = TRUE, border_diag = 1, alpha = 0.05)
#' par(fig = c(0.3, 1.0, 0.73, 0.83), new = TRUE)
#' rlim <- range(dres[, , 1], na.rm = TRUE)
#' plotColorbar(rlim = rlim, at = c(0.2, 0.4, 0.6, 0.8), horiz = TRUE)
#' par(parstart)
#'
#' @seealso \code{\link{plotRel}} and \code{\link{mixMat}} for plotting
#'   relatedness estimates.
#' @export
#'
plotColorbar <- function(rlim = c(0, 1), by = 0.1, at = NULL, horiz = FALSE,
                         col = grDevices::hcl.colors(301, "YlGnBu",
                                                     rev = TRUE), ...) {
  nc <- length(col)
  if (is.null(at)) {
    at <- seq(rlim[1], rlim[2], by)
    at <- at[-c(1, length(at))]
  }
  labtck <- at
  at <- (at - rlim[1])/diff(rlim)
  at <- at*(nc - 1) + 0.5
  ntick <- length(at)
  ltick <- 0.2
  lwd   <- 2.5
  mrg   <- 0.01                        # to align with plotRel if needed
  lim1 <- c(0, 1)
  lim2 <- c(0, nc) + c(-1, 1)*mrg*nc
  if (horiz) {
    xlim <- lim2; ylim <- lim1
  } else {
    xlim <- lim1; ylim <- lim2
  }
  plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE,
       xaxs = "i", yaxs = "i", ...)
  for (i in 1:nc) {
    c2 <- c(i - 1, i)
    if (horiz) {
      xc <- c2;   yc <- lim1
    } else {
      xc <- lim1; yc <- c2
    }
    graphics::polygon(xc[c(1, 2, 2, 1, 1)], yc[c(1, 1, 2, 2, 1)],
                      col = col[i], border = col[i])
  }
  c0  <- rep(c(0, 1 - ltick), each = ntick)
  c1  <- rep(c(ltick, 1),     each = ntick)
  c01 <- rep(at, 2)
  cmid <-
  if (horiz) {
    x0 <- c01; x1 <- c01
    y0 <- c0;  y1 <- c1
  } else {
    x0 <- c0;  x1 <- c1
    y0 <- c01; y1 <- c01
  }
  graphics::segments(x0, y0, x1, y1, col = col[round(nc - at + 0.5)], lwd = lwd)
  if ((nc/2) %in% at) {
    imid <- which(at == nc/2)
    ii <- c(imid, imid + ntick)
    graphics::segments(x0[ii], y0[ii], x1[ii], y1[ii], col = "gray78",
                       lwd = lwd)
  }
  labtck <- rlim[1] + (at - 0.5)*diff(rlim)/(nc - 1)
  side <- if (horiz) 1 else 4
  graphics::mtext(labtck, side, line = 0.3, at = at, las = 1, cex = 0.8)
  graphics::mtext(expression(widehat(r)), 3, line = 0.2, cex = 1.0)
}

#' Combine two matrices
#'
#' Creates a single matrix from two lower or upper triangular matrices.
#'
#' @param mat1,mat2 square matrices of the same dimensions.
#' @param lower a logical value indicating if the resulting matrix should be comprised of lower triangular matrices.
#' @return A square matrix.
#'
#' @examples
#' x <- matrix(1:16, 4)
#' y <- matrix(101:116, 4)
#' x[upper.tri(x, diag = TRUE)] <- NA
#' mixMat(x, y)
#'
#' @seealso \code{\link{plotRel}} and \code{\link{plotColorbar}} for plotting
#'   relatedness estimates.
#' @export
#'
mixMat <- function(mat1, mat2, lower = TRUE) {
  itri <- if (lower) upper.tri(mat1) else lower.tri(mat1)
  mat1[itri] <- t(mat2)[itri]
  return(mat1)
}


