# Morphospace functions -------

# @export
morphospacePCA <- function(PCA, xax, yax, pos.shp, nb.shp = 24,
                           nr.shp = 6, nc.shp = 5, amp.shp = 1,
                           rotate.shp = 0,
                           flipx.shp = FALSE,
                           flipy.shp = FALSE,
                           size.shp = 1, wdw=max(.wdw()),
                           pts.shp = 60,
                           col.shp = "#00000011", border.shp = "#00000055",
                           lwd.shp = 1, plot = TRUE) {
  # we check here, though it shoudl have been before
  if (length(PCA$method)>4 | is.null(PCA$method)) {
    stop("morphospacePCA needs a $method of length <= 5")}
  # we retrive the values corresponding to the two plotted axes and the meanshape
  xy <- PCA$x[, c(xax, yax)]
  rot <- PCA$rotation[, c(xax, yax)]
  mshape <- PCA$mshape
  # we define the position of shapes
  pos <- morphospace_positions(xy, pos.shp = pos.shp, nb.shp = nb.shp,
                    nr.shp = nr.shp, nc.shp = nc.shp)
  # according to the type of morphometrics applied, we switch the method
  # and the way we plot reconstruct shapes (polygon, lines, points for Out, Opn, Ldk)
  # when the object combines different morphometric approaches (up to 4)
  # their size is divided by 2 and the shapes and set (of d) around the (x; y) coordinates of pos.shp
  method <- PCA$method
  lm <- length(method)
  if (length(rotate.shp) != lm){
    rotate.shp <- rep(rotate.shp, lm)
  }
  if (length(flipx.shp) != lm){
    flipx.shp <- rep(flipx.shp, lm)
  }
  if (length(flipy.shp) != lm){
    flipy.shp <- rep(flipy.shp, lm)
  }

  if (length(size.shp)!=lm) size.shp <- rep(size.shp[1], lm)
  size.shp.final <- (size.shp*wdw/14) / ifelse(lm<2, 1, 2)
  d <- mean(size.shp.final) / 2
  # here we define the translation x and y for every sub-morphoshape
  # and the coe to retrieve
  if (lm==1){
    dx <- 0
    dy <- 0}
  if (lm==2){ #met1 over met2 - h center
    dx <- c(0, 0)
    dy <- c(d, -d)}
  if (lm==3){ #podium arrangement
    dx <- c(0, -d, d)
    dy <- c(d, -d, -d)}
  if (lm==4){ #form top left, clockwise
    dx <- c(-d, d, -d, d)
    dy <- c(d, d, -d, -d)}
  # indices of successive coe to select
  if (lm==1){
    col.start <- 1
    col.end   <- length(mshape)
  } else {
    col.start <- cumsum(PCA$cuts) - PCA$cuts + 1
    col.end   <- cumsum(PCA$cuts)}
  # not very satisfactory...
  # hack in case of multi
  if (!plot) SHP <- list()
  for (i in seq(along=method)){
    shp <- NULL
    plot.method <- NULL
    ids <- col.start[i]:col.end[i]
    # outlines
    # efourier
    if (method[i] == "efourier") {
      shp <- PCA2shp_efourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids],
                              amp.shp = amp.shp, pts.shp = pts.shp)
      shp <- lapply(shp, coo_close)
      plot.method <- "poly"}
    # rfourier
    if (method[i] == "rfourier") {
      shp <- PCA2shp_rfourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids],
                              amp.shp = amp.shp, pts.shp = pts.shp)
      shp <- lapply(shp, coo_close)
      plot.method <- "poly"}
    # sfourier
    if (method[i] == "sfourier") {
      shp <- PCA2shp_sfourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids],
                              amp.shp = amp.shp, pts.shp = pts.shp)
      shp <- lapply(shp, coo_close)
      plot.method <- "poly"}
    # tfourier
    if (method[i] == "tfourier") {
      shp <- PCA2shp_tfourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids],
                              amp.shp = amp.shp, pts.shp = pts.shp)
      shp <- lapply(shp, coo_close)
      plot.method <- "poly"}
    ### open outlines
    # dfourier
    if (method[i] == "dfourier") {
      shp <- PCA2shp_dfourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids],
                              amp.shp = amp.shp, pts.shp = pts.shp)
      plot.method <- "lines"}
    # opoly
    if (method[i] == "opoly") {
      shp <- PCA2shp_polynomials(pos = pos, rot = rot[ids, ], mshape = mshape[ids],
                                 amp.shp = amp.shp, pts.shp = pts.shp, ortho = TRUE,
                                 baseline1 = PCA$baseline1[1:2 + (i-1)*2],
                                 baseline2 = PCA$baseline2[1:2 + (i-1)*2])
      plot.method <- "lines"}
    # npoly
    if (method[i] == "npoly") {
      shp <- PCA2shp_polynomials(pos = pos, rot = rot[ids, ], mshape = mshape[ids],
                                 amp.shp = amp.shp, pts.shp = pts.shp, ortho = FALSE,
                                 baseline1 = PCA$baseline1[1:2 + (i-1)*2],
                                 baseline2 = PCA$baseline2[1:2 + (i-1)*2])
      plot.method <- "lines"}
    ### configuration of landmarks
    if (method[i] == "procrustes") {
      shp <- PCA2shp_procrustes(pos = pos, rot = rot[ids, ],
                                mshape = mshape[ids],
                                amp.shp = amp.shp)
      plot.method <- "points"}
    ### Then...
    # we template shapes
    shp <- lapply(shp, coo_template, size = size.shp.final[i])
    # since coo_template does not center shapes but the bounding box
    shp <- lapply(shp, coo_center)
    # we rotate shapes
    shp <- lapply(shp, coo_rotate, rotate.shp[i])
    # we flip (if required)
    if (flipx.shp[i]) shp <- lapply(shp, coo_flipx)
    if (flipy.shp[i]) shp <- lapply(shp, coo_flipy)
    # we translate shapes
    if (plot) { # to translate only for morphospace PCA, not PCcontrib, etc.
      for (s in 1:length(shp)) {
        shp[[s]] <- coo_trans(shp[[s]], pos[s, 1] + dx[i], pos[s, 2] + dy[i])}
    } else {
      SHP[[i]] <- shp
    }
    # otherwise, we plot the morphospace
    if (plot) {
      if (plot.method == "poly") {
        garbage <- lapply(shp, coo_draw, col = col.shp, border = border.shp, lwd = lwd.shp,
                          points = FALSE, centroid = FALSE, first.point = FALSE)}
      if (plot.method == "lines"){
        garbage <- lapply(shp, lines, col = border.shp, lwd = lwd.shp * 2)}
      if (plot.method == "points"){
        garbage <- lapply(shp, points, col = border.shp, cex = lwd.shp*0.25, pch=20)
        if (!is.null(PCA$links)) lapply(shp, function(x) ldk_links(x, PCA$links, col="grey90"))
      }
    }
  }
  if (!plot) SHP else invisible(shp)
}

# @export
morphospaceLDA <- function(LDA, xax, yax, pos.shp, nb.shp = 24,
                           nr.shp = 6, nc.shp = 5, amp.shp = 1, size.shp = 1, pts.shp = 60,
                           col.shp = "#00000011", border.shp = "#00000055") {

  xy <- LDA$mod.pred$x[, c(xax, yax)]
  rot <- LDA$LDs[, c(xax, yax)]
  mshape <- LDA$mshape

  # we fill any removed variables with 0s
  r <- LDA$removed
  if (length(r) > 0) {
    m2 <- matrix(rep(0, length(r) * 2), nrow = length(r),
                 byrow = TRUE, dimnames = list(names(r), colnames(rot)))
    m3 <- rbind(rot, m2)
    rot <- m3[match(names(mshape), rownames(m3)), ]
  }

  # we define the position of shapes
  pos <- morphospace_positions(xy, pos.shp = pos.shp, nb.shp = nb.shp,
                    nr.shp = nr.shp, nc.shp = nc.shp)
  # according to the type of morphometrics applied, we
  # reconstruct shapes
  method <- LDA$method
  ## outlines
  if (method == "efourier") {
    shp <- LDA2shp_efourier(pos = pos, rot = rot, mshape = mshape,
                            amp.shp = amp.shp, pts.shp = pts.shp)
    cd <- TRUE
  }

  shp <- lapply(shp, coo_template, size = (size.shp*max(.wdw())/14))
  shp <- lapply(shp, coo_close)
  for (i in 1:length(shp)) {
    shp[[i]] <- coo_trans(shp[[i]], pos[i, 1], pos[i, 2])
  }
  if (cd) {
    garbage <- lapply(shp, coo_draw, col = col.shp, border = border.shp,
                      points = FALSE, centroid = FALSE, first.point = TRUE)
  } else {
    garbage <- lapply(shp, lines, col = border.shp)
  }
  invisible(shp)
}

#' Calculates nice positions on a plane for drawing shapes
#'
#' @param xy a matrix of points typically from a PCA or other multivariate method on
#' which morphospace can be calculated
#' @param pos.shp how shapes should be positionned: \code{range} of xy,
#' \code{full} extent of the plane, \code{circle} as a rosewind,
#' on \code{xy} values provided, \code{range_axes} on the range of xy
#' but on the axes, \code{full_axes} same thing but on (0.85) range of the axes.
#' You can also directly pass a matrix (or a data.frame)
#' with columns named \code{("x", "y")}.
#' @param nb.shp the total number of shapes
#' @param nr.shp the number of rows to position shapes
#' @param nc.shp the number of cols to position shapes
#' @param circle.r.shp if circle, its radius
#' @return a data.frame of positions
#' @details See \link{plot.PCA} for self-speaking examples
#' @export
morphospace_positions <- function(xy, pos.shp = c("range", "full", "circle", "xy",
                                       "range_axes", "full_axes")[1],
                       nb.shp = 12, nr.shp = 6, nc.shp = 5, circle.r.shp) {
  if (is.data.frame(pos.shp) | is.matrix(pos.shp)) {
    return(as.matrix(pos.shp))
  }
  if (pos.shp == "xy") {
    return(xy)
  }
  if (pos.shp == "circle") {
    if (missing(circle.r.shp)) {
      # mean distance from origin
      circle.r.shp <- coo_centsize(xy)
    }
    t <- seq(0, 2 * pi, len = nb.shp + 1)[-(nb.shp + 1)]
    pos <- cbind(circle.r.shp * cos(t), circle.r.shp * sin(t))
    colnames(pos) <- c("x", "y")  # pure cosmetics
    return(pos)
  }
  if (pos.shp == "range") {
    pos <- expand.grid(seq(min(xy[, 1]), max(xy[, 1]), len = nr.shp),
                       seq(min(xy[, 2]), max(xy[, 2]), len = nc.shp))
    pos <- as.matrix(pos)
    colnames(pos) <- c("x", "y")  # pure cosmetics
    return(pos)
  }
  if (pos.shp == "full") {
    w <- par("usr")
    pos <- expand.grid(seq(w[1] * 0.85, w[2] * 0.85, len = nr.shp),
                       seq(w[3] * 0.85, w[4] * 0.85, len = nc.shp))
    pos <- as.matrix(pos)
    colnames(pos) <- c("x", "y")  # pure cosmetics
    return(pos)
  }
  if (pos.shp == "range_axes") {
    pos <- matrix(c(0, 0,
                    min(xy[, 1]), 0,
                    max(xy[, 1]), 0,
                    0, min(xy[, 2]),
                    0, max(xy[, 2])),
                  byrow=TRUE, ncol=2)
    colnames(pos) <- c("x", "y")  # pure cosmetics
    return(pos)
  }
  if (pos.shp == "full_axes") {
    w <- par("usr") * 0.85
    pos <- matrix(c(0, 0,
                    w[1], 0,
                    w[2], 0,
                    0, w[3],
                    0, w[4]),
                  byrow=TRUE, ncol=2)
    colnames(pos) <- c("x", "y")  # pure cosmetics
    return(pos)
  }
  # if a non-valid method is passed
  return(xy)
}

# Domestic -----------
# stupid function
# @export
.mprod <- function(m, s) {
  res <- m
  for (i in 1:ncol(m)) {
    res[, i] <- m[, i] * s[i]
  }
  return(res)
}

# Calculates shapes from PC plane: e/r/tfourier
#
# @param pos the position on two PC axis
# @param rot the corresponding loadings
# @param mshape the meanshape
# @param amp.shp amplification factor for the shape deformation
# @param pts.shp number of points to reconstruct the shape
# @rdname PCA2shp_fourier
# @export
PCA2shp_efourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
  if (ncol(pos) != ncol(rot))
    stop("'rot' and 'pos' must have the same ncol")
  if (length(mshape) != nrow(rot))
    stop("'mshape' and ncol(rot) lengths differ")
  nb.h <- length(mshape)/4
  n <- nrow(pos)
  # we prepare the array
  res <- list()
  for (i in 1:n) {
    ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
    coe <- mshape + apply(ax.contrib, 1, sum)
    xf <- coeff_split(coe)
    coo <- efourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
    # reconstructed shapes are translated on their centroid if
    # (trans) {
    dx <- pos[i, 1] - coo_centpos(coo)[1]
    dy <- pos[i, 2] - coo_centpos(coo)[2]
    coo <- coo_trans(coo, dx, dy)
    # }
    res[[i]] <- coo
  }
  return(res)
}

# @rdname PCA2shp_fourier
# @export
PCA2shp_rfourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
  if (ncol(pos) != ncol(rot))
    stop("'rot' and 'pos' must have the same ncol")
  if (length(mshape) != nrow(rot))
    stop("'mshape' and ncol(rot) lengths differ")
  nb.h <- length(mshape)/2
  n <- nrow(pos)
  # we prepare the array
  res <- list()
  for (i in 1:n) {
    ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
    coe <- mshape + apply(ax.contrib, 1, sum)
    xf <- coeff_split(coe, cph = 2)
    coo <- rfourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
    # reconstructed shapes are translated on their centroid if
    # (trans) {
    dx <- pos[i, 1] - coo_centpos(coo)[1]
    dy <- pos[i, 2] - coo_centpos(coo)[2]
    coo <- coo_trans(coo, dx, dy)
    # }
    res[[i]] <- coo
  }
  return(res)
}

# @rdname PCA2shp_fourier
# @export
PCA2shp_sfourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
  if (ncol(pos) != ncol(rot))
    stop("'rot' and 'pos' must have the same ncol")
  if (length(mshape) != nrow(rot))
    stop("'mshape' and ncol(rot) lengths differ")
  nb.h <- length(mshape)/2
  n <- nrow(pos)
  # we prepare the array
  res <- list()
  for (i in 1:n) {
    ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
    coe <- mshape + apply(ax.contrib, 1, sum)
    xf <- coeff_split(coe, cph = 2)
    coo <- sfourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
    # reconstructed shapes are translated on their centroid if
    # (trans) {
    dx <- pos[i, 1] - coo_centpos(coo)[1]
    dy <- pos[i, 2] - coo_centpos(coo)[2]
    coo <- coo_trans(coo, dx, dy)
    # }
    res[[i]] <- coo
  }
  return(res)
}
# @rdname PCA2shp_fourier
# @export
PCA2shp_tfourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
  if (ncol(pos) != ncol(rot))
    stop("'rot' and 'pos' must have the same ncol")
  if (length(mshape) != nrow(rot))
    stop("'mshape' and ncol(rot) lengths differ")
  nb.h <- length(mshape)/2
  n <- nrow(pos)
  # we prepare the array
  res <- list()
  for (i in 1:n) {
    ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
    coe <- mshape + apply(ax.contrib, 1, sum)
    xf <- coeff_split(coe, cph = 2)
    coo <- tfourier_i(xf, nb.h = nb.h, nb.pts = pts.shp,
                      force2close = TRUE)
    # reconstructed shapes are translated on their centroid if
    # (trans) {
    dx <- pos[i, 1] - coo_centpos(coo)[1]
    dy <- pos[i, 2] - coo_centpos(coo)[2]
    coo <- coo_trans(coo, dx, dy)
    # }
    res[[i]] <- coo
  }
  return(res)
}

PCA2shp_dfourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
  if (ncol(pos) != ncol(rot))
    stop("'rot' and 'pos' must have the same ncol")
  if (length(mshape) != nrow(rot))
    stop("'mshape' and ncol(rot) lengths differ")
  nb.h <- length(mshape)/2
  n <- nrow(pos)
  # we prepare the array
  res <- list()
  for (i in 1:n) {
    ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
    coe <- mshape + apply(ax.contrib, 1, sum)
    xf <- coeff_split(coe, cph=2)
    #     coo <- dfourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
    coo <- dfourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
    # reconstructed shapes are translated on their centroid if
    # (trans) {
    dx <- pos[i, 1] - coo_centpos(coo)[1]
    dy <- pos[i, 2] - coo_centpos(coo)[2]
    coo <- coo_trans(coo, dx, dy)
    # }
    res[[i]] <- coo
  }
  return(res)
}

# Calculates shapes from PC plane: polynomials
#
# @param pos the position on two PC axis
# @param rot the corresponding loadings
# @param mshape the meanshape
# @param amp.shp amplification factor for the shape deformation
# @param pts.shp number of points to reconstruct the shape
# @param ortho logical whether working with raw or orthogonal polynomials
# @param baseline1 the (x; y) coordinates of the first baseline point
# @param baseline2 the (x; y) coordinates of the second baseline point
# @export
PCA2shp_polynomials <- function(pos, rot, mshape, amp.shp = 1,
                                pts.shp = 60, ortho, baseline1, baseline2) {
  if(ortho) {
    method_i <- opoly_i
  } else {
    method_i <- npoly_i
  }
  if (ncol(pos) != ncol(rot))
    stop("'rot' and 'pos' must have the same ncol")
  if (length(mshape) != nrow(rot))
    stop("'mshape' and ncol(rot) lengths differ")
  degree <- length(mshape)
  n <- nrow(pos)
  # an empy pol object
  pol <- list(coeff = rep(NA, degree), ortho = ortho, baseline1 = baseline1, baseline2 = baseline2)
  # we prepare the array
  res <- list()
  for (i in 1:n) {
    ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
    pol$coeff <- mshape + apply(ax.contrib, 1, sum)
    coo <- method_i(pol, nb.pts = pts.shp, reregister = TRUE)
    pol$coeff <- rep(NA, degree)
    # reconstructed shapes are translated on their centroid if
    # (trans) {
    dx <- pos[i, 1] - coo_centpos(coo)[1]
    dy <- pos[i, 2] - coo_centpos(coo)[2]
    coo <- coo_trans(coo, dx, dy)
    res[[i]] <- coo
  }
  # }
  return(res)
}

# Calculates shapes from PC plane: (aligned) landmarks
#
# @param pos the position on two PC axis
# @param rot the corresponding loadings
# @param mshape the meanshape
# @param amp.shp amplification factor for the shape deformation
# @export
PCA2shp_procrustes <- function(pos, rot, mshape, amp.shp = 1) {
  if (ncol(pos) != ncol(rot))
    stop("'rot' and 'pos' must have the same ncol")
  if (length(mshape) != nrow(rot))
    stop("'mshape' and ncol(rot) lengths differ")
  n <- nrow(pos)
  # we prepare the array
  res <- list()
  for (i in 1:n) {
    ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
    shape.i <- mshape + apply(ax.contrib, 1, sum)
    coo <- matrix(shape.i, ncol = 2, byrow = FALSE)
    # reconstructed shapes are translated on their centroid if
    # (trans) {
    dx <- pos[i, 1] - coo_centpos(coo)[1]
    dy <- pos[i, 2] - coo_centpos(coo)[2]
    coo <- coo_trans(coo, dx, dy)
    res[[i]] <- coo
  }
  # }
  return(res)
}

# todo r/t
# Calculates shapes from LD plane: e/r/tfourier
#
# @param pos the position on two PC axis
# @param rot the (unstardized) loadings
# @param mshape the meanshape
# @param amp.shp amplification factor for the shape deformation
# @param pts.shp number of points to reconstruct the shape
# @rdname LDA2shp_fourier
# @export

# with layerize, I think LDA2shp_* equals PCA2shp_*
# apropos("PCA2shp_") %>% paste(gsub("PCA", "LDA", .), "<-", .) %>% cat(sep="\n")
LDA2shp_dfourier    <- PCA2shp_dfourier
LDA2shp_efourier    <- PCA2shp_efourier
LDA2shp_polynomials <- PCA2shp_polynomials
LDA2shp_procrustes  <- PCA2shp_procrustes
LDA2shp_rfourier    <- PCA2shp_rfourier
LDA2shp_sfourier    <- PCA2shp_sfourier
LDA2shp_tfourier    <- PCA2shp_tfourier

# LDA2shp_efourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
#   if (ncol(pos) != ncol(rot))
#     stop("'rot' and 'pos' must have the same ncol")
#   if (length(mshape) != nrow(rot))
#     stop("'mshape' and ncol(rot) lengths differ")
#   nb.h <- length(mshape)/4
#   n <- nrow(pos)
#   # we prepare the array
#   res <- list()
#   for (i in 1:n) {
#     ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
#     coe <- mshape + apply(ax.contrib, 1, sum)
#     xf <- coeff_split(coe)
#     coo <- efourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
#     # reconstructed shapes are translated on their centroid if
#     # (trans) {
#     dx <- pos[i, 1] - coo_centpos(coo)[1]
#     dy <- pos[i, 2] - coo_centpos(coo)[2]
#     coo <- coo_trans(coo, dx, dy)
#     # }
#     res[[i]] <- coo
#   }
#   return(res)
# }
#
# # @rdname LDA2shp_fourier
# # @export
# LDA2shp_rfourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
#   if (ncol(pos) != ncol(rot))
#     stop("'rot' and 'pos' must have the same ncol")
#   if (length(mshape) != nrow(rot))
#     stop("'mshape' and ncol(rot) lengths differ")
#   nb.h <- length(mshape)/4
#   n <- nrow(pos)
#   # we prepare the array
#   res <- list()
#   for (i in 1:n) {
#     ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
#     coe <- mshape + apply(ax.contrib, 1, sum)
#     xf <- coeff_split(coe)
#     coo <- rfourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
#     # reconstructed shapes are translated on their centroid if
#     # (trans) {
#     dx <- pos[i, 1] - coo_centpos(coo)[1]
#     dy <- pos[i, 2] - coo_centpos(coo)[2]
#     coo <- coo_trans(coo, dx, dy)
#     # }
#     res[[i]] <- coo
#   }
#   return(res)
# }
#
# # @rdname LDA2shp_fourier
# # @export
# LDA2shp_tfourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
#   if (ncol(pos) != ncol(rot))
#     stop("'rot' and 'pos' must have the same ncol")
#   if (length(mshape) != nrow(rot))
#     stop("'mshape' and ncol(rot) lengths differ")
#   nb.h <- length(mshape)/4
#   n <- nrow(pos)
#   # we prepare the array
#   res <- list()
#   for (i in 1:n) {
#     ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
#     coe <- mshape + apply(ax.contrib, 1, sum)
#     xf <- coeff_split(coe)
#     coo <- tfourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
#     # reconstructed shapes are translated on their centroid if
#     # (trans) {
#     dx <- pos[i, 1] - coo_centpos(coo)[1]
#     dy <- pos[i, 2] - coo_centpos(coo)[2]
#     coo <- coo_trans(coo, dx, dy)
#     # }
#     res[[i]] <- coo
#   }
#   return(res)
# }
#
# LDA2shp_sfourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
#   if (ncol(pos) != ncol(rot))
#     stop("'rot' and 'pos' must have the same ncol")
#   if (length(mshape) != nrow(rot))
#     stop("'mshape' and ncol(rot) lengths differ")
#   nb.h <- length(mshape)/4
#   n <- nrow(pos)
#   # we prepare the array
#   res <- list()
#   for (i in 1:n) {
#     ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
#     coe <- mshape + apply(ax.contrib, 1, sum)
#     xf <- coeff_split(coe)
#     coo <- sfourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
#     # reconstructed shapes are translated on their centroid if
#     # (trans) {
#     dx <- pos[i, 1] - coo_centpos(coo)[1]
#     dy <- pos[i, 2] - coo_centpos(coo)[2]
#     coo <- coo_trans(coo, dx, dy)
#     # }
#     res[[i]] <- coo
#   }
#   return(res)
# }

##### end morphospaces
