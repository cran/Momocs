# 5. Coe / OutCoe / OpnCoe plotters
# ------------------------------------------
#' Boxplot of morphometric coefficients
#'
#' Explores the distribution of coefficient values.
#'
#' @param x the \link{Coe} object
#' @param ... useless here
#' @return a ggplot2 object
#' @aliases boxplot.Coe
#' @family Coe_graphics
#' @examples
#' # on OutCoe
#' bot %>% efourier(9) %>% rm_harm(1) %>% boxplot()
#'
#' data(olea)
#' op <- opoly(olea)
#' boxplot(op)
#' @export
boxplot.OutCoe <- function(x, ...){
  # we convert to a data.frame
  x$coe %>%
    .melt_mat %>%
    dplyr::mutate(coeff=factor(as.numeric(substr(key, 2, nchar(key)))),
                  harm=substr(key, 1, 1)) %>%
    ggplot() +
    aes(coeff, value) +
    geom_boxplot() +
    facet_grid(harm~., scales = "free_y")
}

#'@export
boxplot.OpnCoe <- function(x, ...){
  # xfourier case should be treated as a Out
  if (grepl("fourier", x$method)) {
    return(boxplot.OutCoe(x))
  }
  # otherwise...
  x$coe %>%
    .melt_mat %>%
    ggplot() +
    aes(key, value) +
    geom_boxplot()
}

#' Harmonic contribution to shape
#'
#' Calculates contribution of harmonics to shape. The amplitude of every coefficients
#' of a given harmonic is multiplied by the coefficients provided and the resulting
#' shapes are reconstructed and plotted. Naturally, only works on Fourier-based methods.
#' @param Coe a \code{\link{Coe}} object (either \code{OutCoe} or (soon) \code{OpnCoe})
#' @param id the id of a particular shape, otherwise working on the meanshape
#' @param harm.r range of harmonics on which to explore contributions
#' @param amp.r a vector of numeric for multiplying coefficients
#' @param main a title for the plot
#' @param xlab a title for the x-axis
#' @param ylab a title for the y-axis
#' @param ... additional parameter to pass to \code{\link{coo_draw}}
#' @rdname harm.contrib
#' @return a plot
#' @family Coe_graphics
#' @examples
#' data(bot)
#' bot.f <- efourier(bot, 12)
#' hcontrib(bot.f)
#' hcontrib(bot.f, harm.r=3:10, amp.r=1:8, col="grey20",
#'    main="A huge panel")
#' @export
hcontrib <- function(Coe, ...){
  UseMethod("hcontrib")
  }

#' @rdname harm.contrib
#' @export
hcontrib.OutCoe <- function(Coe,
                            id,
                            harm.r,
                            amp.r   = c(0, 0.5, 1, 2, 5, 10),
                            main="Harmonic contribution to shape",
                            xlab="Harmonic rank",
                            ylab="Amplification factor", ...){
  x <- Coe
  # we handle the method
  p <- pmatch(tolower(x$method), c("efourier", "rfourier", "tfourier"))
  if (is.na(p)) { warning("unvalid method. efourier is used")
  } else {
    method.i <- switch(p, efourier_i, rfourier_i, tfourier_i)}
  # we deduce the number of coefficient / harmonic, and their number
  cph <- ifelse(p==1, 4, 2)
  nb.h <- ncol(x$coe)/cph
  # we handle for missing harm.r
  if (missing(harm.r)) harm.r <- 1:ifelse(nb.h > 6, 6, nb.h)
  # if id is provided, we work on it, otherwise, on the average shape
  if (missing(id)){
    coe <- apply(x$coe, 2, mean)
    message("no 'id' provided, working on the meanshape")
  } else {
    coe <- x$coe[id, ]}
  # we prepare a xf to feed the method.i functions
  xf <- coeff_split(coe, nb.h, cph)
  # we prepare a neutral amplification factor
  mult <- rep(1, nb.h)
  # the core below
  shp <- list()
  p <- 1 # kinda dirty
  # we loop over harm.r and amp.r by just multiplying xf
  # by a given vector, all set to 1 except the harmonic that has to
  # amplified
  for (j in seq(along=harm.r)){
    for (i in seq(along=amp.r)){
      mult.loc    <- mult
      mult.loc[harm.r[j]] <- amp.r[i]
      xfi <- lapply(xf, function(x) x*mult.loc)
      shp[[p]] <- method.i(xfi)
      p <- p+1}}
  # graphics start here
  # we borrow this block to PC.contrib
  # except the expand.grid and coo_trans that needed to be "transposed"
  xs <- 1:length(harm.r) - 0.5
  ys <- rev(1:length(amp.r) - 0.5)
  plot(NA, xlim=c(0, length(harm.r)), ylim=c(0, length(amp.r)),
       asp=1, frame=FALSE, axes=FALSE,
       main=main, xlab=xlab, ylab=ylab)
  axis(1, at = xs, labels = harm.r)
  axis(2, at = ys, labels = amp.r, las=1)
  # we template the size of the shapes
  shp <- lapply(shp, coo_close)
  shp <- lapply(shp, coo_template, 0.95)
  # here we prepare and apply the translation values
  trans <- expand.grid(ys, xs)
  colnames(trans) <- c("x", "y")
  for (i in seq(along=shp)){
    shp[[i]] <- coo_trans(shp[[i]], trans[i, 2], trans[i, 1])}
  # we finally plot the shapes
  gc <- lapply(shp, coo_draw, centroid = FALSE, first.point=FALSE, ...)
  invisible(list(shp=shp, trans=trans))}

# hcontrib.Opn (dct) # todo

##### end graphics Coe
