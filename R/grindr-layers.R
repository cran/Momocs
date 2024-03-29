# .orange   <- "#ffa500"
# .grey90   <- "#e5e5e5"
# .grey60   <- "#999999"

.layerize_df <- function(x, f, axes=c(1, 2), palette=pal_qual){
  # grab the selected columns
  # not optimal - to check
  # if (!is.null(x$x))
  xy <- x$x[, axes]
  # else
  # xy <- x[, axes]
  # prepare a factor
  if (missing(f) | is.null(f)){ # no factor or NULL provided
    f <- factor(rep(1, nrow(x$x)))
    colors_groups <- rep(par("fg"), nlevels(f))
    colors_rows   <- colors_groups[f]
  } else {         # something provided, handle with fac_dispatcher
    f <- x %>% fac_dispatcher(f)

    # handle palettes
    if (missing(palette) | is.null(palette)){
      if (is.numeric(f))
        palette <- pal_seq
      else
        palette <- pal_qual
    }

    if (is.numeric(f)){
      colors_groups <- NA
      colors_rows <- f %>% .normalize()  %>%
        cut(breaks = 1e3)  %>% as.numeric() %>% `[`(palette(1e3), .)
    }
    if (is.factor(f)){
      colors_groups <- palette(nlevels(f))
      colors_rows   <- colors_groups[f]
    }
  }

  # NA handling
  nas <- which(is.na(f))
  if (length(nas)>0){
    xy <- xy[-nas, ]
    f  <- f[-nas]
    colors_groups <- colors_groups[-nas]
    colors_rows   <- colors_rows[-nas]
  }
  # return as a list and
  # add potentially useful other components
  list(xy=xy, f=f,
       colors_groups=colors_groups,
       colors_rows=colors_rows,
       palette=palette)
}


.layerize_PCA <- function(x, f, axes=c(1, 2), palette=pal_qual){
  # grab the selected columns
  xy <- x$x[, axes]
  # prepare a factor
  if (missing(f)){ # no factor provided
    f <- factor(rep(1, nrow(x$x)))
    colors_groups <- rep(par("fg"), nlevels(f))
    colors_rows   <- colors_groups[f]
  } else {         # something provided, handle with fac_dispatcher
    f <- x %>% fac_dispatcher(f)

    # handle palettes
    if (missing(palette) | is.null(palette)){
      if (is.numeric(f))
        palette <- pal_seq
      else
        palette <- pal_qual
    }

    # NA handling
    nas <- which(is.na(f))
    if (length(nas)>0){
      xy <- xy[-nas, ]
      f  <- f[-nas]
      message("some NAs were dropped")
    }

    if (is.numeric(f)){
      colors_groups <- NA
      colors_rows <- f %>% .normalize()  %>%
        cut(breaks = 1e3)  %>% as.numeric() %>% `[`(palette(1e3), .)
    }
    if (is.factor(f)){
      colors_groups <- palette(nlevels(f))
      colors_rows   <- colors_groups[f]
    }
  }

  # NA handling
  nas <- which(is.na(f))
  if (length(nas)>0){
    xy <- xy[-nas, ]
    f  <- f[-nas]
    colors_groups <- colors_groups[-nas]
    colors_rows   <- colors_rows[-nas]
  }
  # return as a list and
  # add potentially useful other components
  list(xy=xy, f=f,
       colors_groups=colors_groups, colors_rows=colors_rows,
       object="PCA", axes=axes, palette=palette,
       method=x$method, mshape=x$mshape, cuts=x$cuts,
       eig=x$eig, sdev=x$sdev, rotation=x$rotation[, axes],
       baseline1=x$baseline1, baseline2=x$baseline2,
       links=x$links)
}


.layerize_LDA <- function(x, f=x$f, axes=c(1, 2), palette=pal_qual){
  # restore LDs if some were dropped because constant or collinear
  x$LDs <- .restore_LDs(x)

  if (ncol(x$LDs) < 2){
    message("* Only two levels, so a single LD and preparing for an histogram")
    xy <- x$mod.pred$x[, 1, drop=FALSE]
    axes <- 1
  } else {
    xy <- x$mod.pred$x[, axes]
  }

  # prepare a factor
  # if (missing(f)){ # no factor provided
  #   f <- factor(rep(1, nrow(x$x)))
  #   colors_groups <- rep(par("fg"), nlevels(f))
  #   colors_rows   <- colors_groups[f]
  # } else {         # something provided, handle with fac_dispatcher
  colors_groups <- palette(nlevels(f))
  colors_rows   <- colors_groups[f]

  # NA handling
  nas <- which(is.na(f))
  if (length(nas)>0){
    xy <- xy[-nas, ]
    f  <- f[-nas]
    colors_groups <- colors_groups[-nas]
    colors_rows   <- colors_rows[-nas]
  }
  # return as a list and
  # add potentially useful other components
  list(xy=xy, f=f,
       colors_groups=colors_groups, colors_rows=colors_rows,
       object="PCA", axes=axes, palette=palette,
       method=x$method, mshape=x$mshape, cuts=x$cuts,
       eig=x$eig, sdev=x$mod$svd,
       rotation=x$LDs[, axes, drop=FALSE],
       LDs=x$LDs,
       baseline1=x$baseline1, baseline2=x$baseline2,
       links=x$links)
}

.layerize_NMDS <- function(x, f=NULL, axes=c(1, 2), palette=pal_qual){
  x0 <- list(x=x$points, fac=x$fac)
  res <- .layerize_df(x0, f=f, axes=axes, palette=palette)
  res$axes <- axes
  res
}

.layerize_MDS <- function(x, f=NULL, axes=c(1, 2), palette=pal_qual){
  x0 <- list(x=x$x, fac=x$fac)
  res <- .layerize_df(x0, f=f, axes=axes, palette=palette)
  res$axes <- axes
  res
}

# LDA(bp, ~type) %>% .layerize_LDA() -> x


### __plot_PCA__ -------------------------------------------
#' PCA plot using grindr layers
#'
#' Quickly vizualise [PCA] objects and friends and build customs plots
#' using the [layers]. See examples.
#'
#' @note This approach will replace \link{plot.PCA} (and `plot.lda` in further versions.
#' This is part of `grindr` approach that may be packaged at some point. All comments are welcome.
#'
#' @param x a [PCA] object
#' @param f factor specification to feed [fac_dispatcher]
#' @param axes \code{numeric} of length two to select PCs to use
#' (\code{c(1, 2)} by default)
#' @param palette \code{color palette} to use \code{col_summer} by default
#' @param points `logical` whether to draw this with [layer_points]
#' @param points_transp `numeric` to feed [layer_points] (default:0.25)
#' @param morphospace `logical` whether to draw this using [layer_morphospace_PCA]
#' @param morphospace_position to feed [layer_morphospace_PCA] (default: "range")
#' @param chull `logical` whether to draw this with [layer_chull]
#' @param chullfilled `logical` whether to draw this with [layer_chullfilled]
#' @param labelpoints `logical` whether to draw this with [layer_labelpoints]
#' @param labelgroups `logical` whether to draw this with [layer_labelgroups]
#' @param legend `logical` whether to draw this with [layer_legend]
#' @param title `character` if specified, fee [layer_title] (default to `""`)
#' @param center_origin `logical` whether to center origin
#' @param zoom `numeric` zoom level for the frame (default: 0.9)
#' @param eigen `logical` whether to draw this using [layer_eigen]
#' @param box `logical` whether to draw this using [layer_box]
#' @param axesnames `logical` whether to draw this using [layer_axesnames]
#' @param axesvar `logical` whether to draw this using [layer_axesvar]
#' @return a plot
#' @family grindr
#'
#' @examples
#' ### First prepare two PCA objects.
#'
#' # Some outlines with bot
#' bp <- bot %>% mutate(fake=sample(letters[1:5], 40, replace=TRUE)) %>%
#' efourier(6) %>% PCA
#' plot_PCA(bp)
#' plot_PCA(bp, ~type)
#' plot_PCA(bp, ~fake)
#'
#' # Some curves with olea
#' op <- olea %>%
#' mutate(s=coo_area(.)) %>%
#' filter(var != "Cypre") %>%
#' chop(~view) %>% opoly(5, nb.pts=90) %>%
#' combine %>% PCA
#' op$fac$s %<>% as.character() %>% as.numeric()
#'
#' op %>% plot_PCA(title="hi there!")
#'
#' ### Now we can play with layers
#' # and for instance build a custom plot
#' # it should start with plot_PCA()
#'
#' my_plot <- function(x, ...){
#'
#' x %>%
#'     plot_PCA(...) %>%
#'     layer_points %>%
#'     layer_ellipsesaxes %>%
#'     layer_rug
#' }
#'
#' # and even continue after this function
#' op %>% my_plot(~var, axes=c(1, 3)) %>%
#'     layer_title("hi there!")
#'
#' # grindr allows (almost nice) tricks like highlighting:
#'
#' # bp %>% .layerize_PCA(~fake) %>%
#' #   layer_frame %>% layer_axes() %>%
#' #   layer_morphospace_PCA() -> x
#'
#' # highlight <- function(x, ..., col_F="#CCCCCC", col_T="#FC8D62FF"){
#' #  args <- list(...)
#' #  x$colors_groups <- c(col_F, col_T)
#' #  x$colors_rows <- c(col_F, col_T)[(x$f %in% args)+1]
#' #  x
#' #  }
#' # x %>% highlight("a", "b") %>% layer_points()
#'
#' # You get the idea.
#' @export
plot_PCA <- function(x,
                     f=NULL,
                     axes=c(1, 2),
                     palette=NULL,
                     points=TRUE,
                     points_transp=1/4,
                     # morphospace
                     morphospace=TRUE,
                     morphospace_position="range",
                     # chulls
                     chull=TRUE,
                     chullfilled=FALSE,
                     # legends
                     labelpoints=FALSE,
                     labelgroups=FALSE,
                     legend=TRUE,
                     # cosmetics (mainly)
                     title="",
                     center_origin=TRUE,
                     zoom=0.9,
                     eigen=TRUE,
                     box=TRUE,
                     axesnames=TRUE, axesvar=TRUE){
  # check ------
  .check(any(class(x)=="PCA"),
         "only supported on PCA objects")

  # prepare ---------------------------
  if (missing(f) | is.null(f)){
    x %<>% .layerize_PCA(axes=axes, palette=palette)
    labelgroups <- legend <- FALSE
  } else {
    x %<>% .layerize_PCA(f, axes=axes, palette=palette)
  }

  # frame
  x %>%
    layer_frame(center_origin=center_origin, zoom = zoom) -> x
  x %>%
    layer_axes() -> x

  # this one did not work, quite inconsistently
  # so we "wait" for the plot to be created
  # thanks to Bill for pointing this
  # x %<>%
  #   layer_frame(center_origin=center_origin, zoom = zoom) %>%
  #   layer_axes()

  # cosmetics
  if (axesnames)
    x %<>% layer_axesnames(name = "PC")

  if (axesvar)
    x %<>% layer_axesvar()

  if (eigen)
    x %<>% layer_eigen()

  if (box)
    x %<>% layer_box()

  # morphospace -----------------------
  if (morphospace & !is.null(x$method))
    x %<>% layer_morphospace_PCA(position = morphospace_position)

  # data ------------------------------
  if (points)
    x %<>% layer_points(transp=points_transp)

  # groups dispersion -----------------
  if (chull & nlevels(x$f)>1)
    x %<>% layer_chull()

  if (chullfilled & nlevels(x$f)>1)
    x %<>% layer_chullfilled()

  # legends
  if (legend){
    if (missing(legend) && is.factor(x$f) && nlevels(x$f)>18){
      message("Many levels! Pass with legend=TRUE if you really want it.")
    } else {
      x %<>% layer_legend()
    }
  }

  if (labelpoints)
    x %<>% layer_labelpoints()

  if (labelgroups)
    x %<>% layer_labelgroups()

  if (title != "")
    x %<>% layer_title(title)

  # propagate
  invisible(x)
}

### __plot_LDA__ -------------------------------------------
#' LDA plot using grindr layers
#'
#' Quickly vizualise [LDA] objects and build customs plots
#' using the [layers]. See examples.
#'
#' @note This approach will replace \link{plot.LDA}.
#' This is part of `grindr` approach that may be packaged at some point. All comments are welcome.
#'
#' @param x [LDA] object
#' @param axes \code{numeric} of length two to select PCs to use
#' (\code{c(1, 2)} by default)
#' @param palette \code{color palette} to use \code{col_summer} by default
#' @param points `logical` whether to draw this with [layer_points]
#' @param points_transp `numeric` to feed [layer_points] (default:0.25)
#' @param morphospace `logical` whether to draw this using [layer_morphospace_PCA]
#' @param morphospace_position to feed [layer_morphospace_PCA] (default: "range")
#' @param chull `logical` whether to draw this with [layer_chull]
#' @param chullfilled `logical` whether to draw this with [layer_chullfilled]
#' @param labelgroups `logical` whether to draw this with [layer_labelgroups]
#' @param legend `logical` whether to draw this with [layer_legend]
#' @param title `character` if specified, fee [layer_title] (default to `""`)
#' @param center_origin `logical` whether to center origin
#' @param zoom `numeric` zoom level for the frame (default: 0.9)
#' @param eigen `logical` whether to draw this using [layer_eigen]
#' @param box `logical` whether to draw this using [layer_box]
#' @param iftwo_layer  function (no quotes) for drawing LD1 when there are two levels.
#' So far, one of [layer_histogram_2] (default) or [layer_density_2]
#' @param iftwo_split to feed `split` argument in [layer_histogram_2] or [layer_density_2]
#' @param axesnames `logical` whether to draw this using [layer_axesnames]
#' @param axesvar `logical` whether to draw this using [layer_axesvar]
#' @return a plot
#' @family grindr
#'
#' @examples
#' ### First prepare an LDA object
#'
#' # Some outlines with bot
#' bl <- bot %>%
#'       # cheap alignement before efourier
#'       coo_align() %>% coo_center %>% coo_slidedirection("left") %>%
#'       # add a fake column
#'       mutate(fake=sample(letters[1:5], 40, replace=TRUE)) %>%
#'       # EFT
#'       efourier(6, norm=FALSE) %>%
#'       # LDA
#'       LDA(~fake)
#'
#' bl %>% plot_LDA %>% layer_morphospace_LDA
#'
#' # Below inherited from plot_PCA and to adapt here.
#' #plot_PCA(bp)
#' #plot_PCA(bp, ~type)
#' #plot_PCA(bp, ~fake)
#'
#' # Some curves with olea
#' #op <- olea %>%
#' #mutate(s=coo_area(.)) %>%
#' #filter(var != "Cypre") %>%
#' #chop(~view) %>% lapply(opoly, 5, nb.pts=90) %>%
#' #combine %>% PCA
#' #op$fac$s %<>% as.character() %>% as.numeric()
#'
#' #op %>% plot_PCA(title="hi there!")
#'
#' ### Now we can play with layers
#' # and for instance build a custom plot
#' # it should start with plot_PCA()
#'
#' #my_plot <- function(x, ...){
#'
#' #x %>%
#' #     plot_PCA(...) %>%
#' #    layer_points %>%
#' #     layer_ellipsesaxes %>%
#' #    layer_rug
#' # }
#'
#' # and even continue after this function
#' # op %>% my_plot(~var, axes=c(1, 3)) %>%
#' #     layer_title("hi there!") %>%
#' #    layer_stars()
#'
#' # You get the idea.
#' @export
plot_LDA <- function(x,
                     axes=c(1, 2),
                     palette=pal_qual,
                     points=TRUE,
                     points_transp=1/4,
                     # morphospace
                     morphospace=FALSE,
                     morphospace_position="range",
                     # chulls
                     chull=TRUE,
                     chullfilled=FALSE,
                     # legends
                     labelgroups=FALSE,
                     legend=TRUE,
                     # cosmetics (mainly)
                     title="",
                     center_origin=TRUE, zoom=0.9,
                     eigen=TRUE,
                     box=TRUE,
                     iftwo_layer=layer_histogram_2,
                     iftwo_split=FALSE,
                     axesnames=TRUE, axesvar=TRUE){

  # check ------
  .check(any(class(x)=="LDA"),
         "only supported on LDA objects")

  # prepare ---------------------------
  x %<>% .layerize_LDA(axes=axes, palette=palette)

  if (length(x$axes) < 2){
    x %>% iftwo_layer(split=iftwo_split)
    return(x)
  }

  .check(all(axes <= ncol(x$LDs)),
         "axes must all be <= number of LDs")

  # frame (see plot_PCA, this no longer worked and we had to "wait"
  # for the plot to be created)
  # x %<>%
  #   layer_frame(center_origin=center_origin, zoom = zoom) %>%
  #   layer_axes()

  x %>%
    layer_frame(center_origin=center_origin, zoom = zoom) -> x
  x %>%
    layer_axes() -> x


  # cosmetics
  if (axesnames)
    x %<>% layer_axesnames(name = "LD")

  if (axesvar)
    x %<>% layer_axesvar()

  if (eigen)
    x %<>% layer_eigen()

  if (box)
    x %<>% layer_box()

  # morphospace -----------------------
  if (morphospace)
    x %<>% layer_morphospace_LDA(position = morphospace_position)

  # data ------------------------------
  if (points)
    x %<>% layer_points(transp=points_transp)

  # groups dispersion -----------------
  if (chull)
    x %<>% layer_chull()

  if (chullfilled)
    x %<>% layer_chullfilled()

  # legends
  if (legend){
    if (missing(legend) && is.factor(x$f) && nlevels(x$f)>18){
      message("Many levels! Pass with legend=TRUE if you really want it.")
    } else {
      x %<>% layer_legend()
    }
  }

  if (labelgroups)
    x %<>% layer_labelgroups()

  if (title != "")
    x %<>% layer_title(title)

  # propagate
  invisible(x)
}


# bot %>% efourier(6) %>% PCA %>% plot_PCA
### _____Layers_____ ---------------------------------------
# frame and options ----------------------------------------

#' grindr layers for multivariate plots
#'
#' Useful layers for building custom
#' mutivariate plots using the cheapbabi approach. See examples.
#'
#' @name layers
#' @rdname layers
#' @seealso grindr_drawers
#' @family grindr
#'
#' @param x a list, typically returned by \link{plot_PCA}
#' @param center_origin \code{logical} whether to center the origin (default \code{TRUE})
#' @param zoom \code{numeric} to change the zoom (default \code{0.9})
#' @return a drawing layer
#' @export
layer_frame <- function(x, center_origin = TRUE, zoom = 0.9){
  # neater par
  old <- par(mar=rep(1/8, 4))
  on.exit(par(old))
  xy <- x$xy
  if (center_origin) {
    w <- (1/zoom) * max(abs(xy))
    plot(NA, xlim = c(-w, w), ylim = c(-w, w), asp = 1,
         axes = FALSE, frame = FALSE, mar=rep(0, 4))
  }
  else {
    w <- (1/zoom) * apply(xy, 2, range)
    plot(xy, xlim = w[, 1], ylim = w[, 2],
         type = "n", asp = 1, axes = FALSE, frame = FALSE)
  }
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
#' @param col color (hexadecimal) to use for drawing components
#' @param lwd linewidth for drawing components
#' @param ... additional options to feed core functions for each layer
layer_axes <- function(x, col="#999999", lwd=1/2, ...){
  # neater par
  # old <- par(mar=rep(1/8, 4))
  # on.exit(par(old))
  # add x=0 and y=0 lines for axes
  abline(h=0, v=0, col=col, lwd=lwd, ...)
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
# cosmetics
layer_ticks <- function(x, col="#333333", cex=3/4, lwd=3/4, ...){
  # neater par
  old <- par(mar=rep(1/8, 4))
  on.exit(par(old))

  # calculate positions
  at_x <- pretty(seq(par("usr")[1], par("usr")[2], length.out=5))
  at_y <- pretty(seq(par("usr")[3], par("usr")[4], length.out=5))

  # calculate gaps
  hx <- .wdw()[1]/200
  hy <- .wdw()[2]/200

  # draw ticjs
  segments(at_x, -hy, at_x, hy, col=col)
  segments(-hx, at_y, hx, at_y, col=col)

  # add positions text
  text(at_x, -strheight(at_x) - hy, labels=at_x, adj=c(0.5, 1), cex=cex)
  text(-strwidth(at_y) - hx, at_y, labels=at_y, adj=c(1, 0.5), cex=cex)

  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
#' @param lty linetype for drawing components
#' @param grid \code{numeric} number of grid to draw
layer_grid <- function(x, col="#999999", lty=3, grid = 3, ...) {
  # neater par
  old <- par(mar=rep(1/8, 4))
  on.exit(par(old))
  m <- max(abs(par("usr")))
  g <- seq(0, m, length = grid)
  g <- c(g[-1] * -1, 0, g[-1])
  abline(h=g, v=g, col = col, lty = lty, ...)
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
#' @param border color (hexadecimal) to use to draw border
layer_box <- function(x, border="#e5e5e5", ...){
  # neater par
  old <- par(mar=rep(1/8, 4))
  on.exit(par(old))
  w <- par("usr")
  # draw the box
  rect(w[1], w[3], w[2], w[4], border=border, ...)
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
layer_fullframe <- function(x, ...){
  x %>%
    layer_frame(...) %>%
    layer_grid() %>%
    layer_axes() %>%
    layer_box() %>%
    layer_axesvar() %>%
    layer_axesnames()
}


layer_semiframe <- function(x, ...){
  x %>%
    layer_frame(...) %>%
    layer_grid() %>%
    layer_axes() %>%
    layer_box()
}

# shapes ---------------------------------------------------
#' @export
#' @rdname layers
#' @param pch to use for drawing components
#' @param cex to use for drawing components
#' @param transp transparency to use (min: 0 defaut:0 max:1)
layer_points <- function(x, pch=20, cex=4/log1p(nrow(x$xy)), transp=0, ...){
  points(x$xy, col=pal_alpha(x$colors_rows, transp=transp), pch=pch, cex=cex, ...)
  # propagate
  invisible(x)
}


# ellipses layers ------------------------------------------
#' @export
#' @rdname layers
#' @param conf \code{numeric} between 0 and 1 for confidence ellipses
#' @param alpha \code{numeric} between 0 and 1 for the transparency of components
layer_ellipses <- function(x, conf=0.5, lwd=1, alpha=0, ...) {
  # if numeric, do not apply this layer but still propagate
  if (!is.factor(x$f))
    invisible(x)
  # if conf and lwd lengths dont match, recycle
  if (length(conf)!=length(lwd))
    lwd <- rep(lwd[1], length(conf))
  # if conf and alpha lengths dont match, recycle
  if (length(conf)!=length(alpha))
    alpha <- rep(alpha[1], length(conf))
  # loop along all levels
  for (i in seq_along(levels(x$f))) {
    xy_i <- x$xy[x$f == levels(x$f)[i], ,drop=FALSE]
    # with less than 3 points in a group,
    # conf_ell would fail
    if (nrow(xy_i)>2) {
      # loop along conf levels
      for (j in seq_along(conf)){
        conf_ell(x = xy_i, conf = conf[j])$ell %>%
          coo_close %>%
          draw_lines(col=x$colors_groups[i] %>%
                       pal_alpha(alpha[j]), ...)
      }
    }
  }
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
layer_ellipsesfilled <- function(x, conf=0.5, lwd=1, alpha=0, ...) {
  # if numeric, do not apply this layer but still propagate
  if (!is.factor(x$f))
    invisible(x)
  # if conf and lwd lengths dont match, recycle
  if (length(conf)!=length(lwd))
    lwd <- rep(lwd[1], length(conf))
  # if conf and alpha lengths dont match, recycle
  if (length(conf)!=length(alpha))
    alpha <- rep(alpha[1], length(conf))

  # loop along all levels
  for (i in seq_along(levels(x$f))) {
    xy_i <- x$xy[x$f == levels(x$f)[i], ,drop=FALSE]
    # with less than 3 points in a group,
    # conf_ell would fail
    if (nrow(xy_i)>2) {
      # loop along conf levels
      for (j in seq_along(conf)){
        conf_ell(x = xy_i, conf = conf[j])$ell %>%
          coo_close %>%
          draw_polygon(fill=NA,
                       col=x$colors_groups[i] %>%
                         pal_alpha(alpha[j]), ...)
      }
    }
  }
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
layer_ellipsesaxes <- function(x, conf=0.5, lwd=1, alpha=0, ...) {
  # if numeric, do not apply this layer but still propagate
  if (!is.factor(x$f))
    invisible(x)
  # if conf and lwd lengths dont match, recycle
  if (length(conf)!=length(lwd))
    lwd <- rep(lwd[1], length(conf))
  # if conf and alpha lengths dont match, recycle
  if (length(conf)!=length(alpha))
    alpha <- rep(alpha[1], length(conf))

  # loop along all levels
  for (i in seq_along(levels(x$f))) {
    xy_i <- x$xy[x$f == levels(x$f)[i], ,drop=FALSE]
    # with less than 3 points in a group,
    # conf_ell would fail
    if (nrow(xy_i)>2) {
      # loop along conf levels
      for (j in seq_along(conf)){
        seg_i <- conf_ell(x = xy_i, conf = conf[j], nb.pts = 360)$seg
        segments(seg_i[1, 1], seg_i[1, 2],
                 seg_i[2, 1], seg_i[2, 2],
                 lwd = lwd[j],
                 col=x$colors_groups[i] %>%
                   pal_alpha(alpha[j]), lend=2, ...)
        segments(seg_i[3, 1], seg_i[3, 2],
                 seg_i[4, 1], seg_i[4, 2],
                 lwd = lwd[j],
                 col=x$colors_groups[i] %>%
                   pal_alpha(alpha[j]), lend=2, ...)

      }
    }
  }
  # propagate
  invisible(x)
}

# chull layers ---------------------------------------------
#' @export
#' @rdname layers
layer_chull <- function(x, ...){
  # if numeric, do not apply this layer but still propagate
  if (!is.factor(x$f))
    invisible(x)
  # loop along levels and draw
  for (i in seq_along(levels(x$f))) {
    coo <- x$xy[x$f == levels(x$f)[i],, drop=FALSE]
    # with less than 3 points in a group,
    # coo_chull would fail
    if (nrow(coo)>2) {
      coo %>% coo_chull() %>% coo_close %>%
        draw_lines(col=x$colors_groups[i], ...)
    }
  }
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
layer_chullfilled <- function(x, alpha=0.8, ...){
  # if numeric, do not apply this layer but still propagate
  if (!is.factor(x$f))
    invisible(x)
  # loop along levels and draw
  for (i in seq_along(levels(x$f))) {
    coo <- x$xy[x$f == levels(x$f)[i],, drop=FALSE]
    if (nrow(coo)>2) {
      # with less than 3 points in a group,
      # coo_chull would fail
      coo %>% coo_chull() %>% coo_close %>%
        draw_polygon(fill=x$colors_groups[i], col=x$colors_groups[i], transp = alpha, ...)
    }
  }
  # propagate
  invisible(x)
}

# stars layers ---------------------------------------------
#' @export
#' @rdname layers
layer_stars <- function(x, alpha=0.5, ...) {
  # if numeric, apply with a local trick on f
  if (!is.factor(x$f)){
    f <- factor(rep(1, length(x$f)))
    colors_groups <- "#000000"
  } else {
    f <- x$f
    colors_groups <- x$colors_groups
  }

  # loop along levels and draw
  for (i in seq_along(levels(f))) {
    xy_i <- x$xy[f == levels(f)[i],, drop=FALSE]
    # to prevent levels with 0
    if (nrow(xy_i) < 1)
      next()
    c_i <- coo_centpos(xy_i)
    for (j in 1:nrow(xy_i)) {
      segments(c_i[1], c_i[2],
               xy_i[j, 1], xy_i[j, 2],
               col = pal_alpha(colors_groups[i], alpha), ...)
    }
  }
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
layer_delaunay <- function(x, ...){
  # if numeric, apply with a local trick on f
  if (!is.factor(x$f)){
    f <- factor(rep(1, length(x$f)))
    colors_groups <- "#000000"
  } else {
    f <- x$f
    colors_groups <- x$colors_groups
  }

  # loop along levels and draw
  for (i in seq_along(levels(f))) {
    xy_i <- x$xy[f == levels(f)[i],, drop=FALSE]
    # with less than 3 points in a group,
    # coo_chull would fail
    if (nrow(xy_i)>2) {
      xy_i %>% geometry::delaunayn() %>%
      {rbind(.[, -1], .[, -2], .[, -3])} %>%
        `[`(-duplicated(.), ) %>%
        draw_links(xy_i, ., col=colors_groups[i], ...)
    }
  }
  # propagate
  invisible(x)
}


# density layers -------------------------------------------
#' @export
#' @rdname layers
#' @param levels_density \code{numeric} number of levels to use to feed \code{MASS::kde2d}
#' @param levels_contour \code{numeric} number of levels to use to feed \code{graphics::contour}
#' @param n \code{numeric} number of grid points to feed \code{MASS::kde2d}
#' @param density \code{logical} whether to draw density estimate
#' @param contour \code{logical} whether to draw contour lines
layer_density <- function(x, levels_density=20,
                          levels_contour=4, alpha=1/3,
                          n = 200, density=TRUE, contour=TRUE) {
  # if numeric, apply with a local trick on f
  if (!is.factor(x$f)){
    f <- factor(rep(1, length(x$f)))
    colors_groups <- "#000000"
  } else {
    f <- x$f
    colors_groups <- x$colors_groups
  }

  # neater par
  old <- par(mar=rep(1/8, 4), xpd=NA)
  on.exit(par(old))
  # loop over levels and use MASS::kde2d
  for (i in seq_along(levels(f))) {
    xy_i <- x$xy[f == levels(f)[i],, drop=FALSE]
    ki <- MASS::kde2d(xy_i[, 1], xy_i[, 2], n = n, lims = c(par("usr")))
    ki$z %<>% .normalize()
    if (density)
      graphics::image(ki$x, ki$y, ki$z, add = TRUE,
                      xlim = range(ki$x), ylim = range(ki$y),
                      col = pal_manual(colors_groups[i], transp=alpha)(levels_density))
    if (contour)
      graphics::contour(ki$x, ki$y, ki$z, add = TRUE,
                        nlevels = levels_contour,
                        col=pal_manual(colors_groups[i], transp=alpha)(levels_contour))
  }
  # propagate
  invisible(x)
}

# label layers ---------------------------------------------
#' @export
#' @rdname layers
#' @param font to feed \link{text}
#' @param abbreviate \code{logical} whether to abbreviate names
layer_labelpoints <- function(x, col=par("fg"), cex=2/3,
                              font=1, abbreviate=FALSE, ...){
  if (missing(col))
    col <- x$colors_rows

  # prepare labs
  labs <- rownames(x$xy)
  if (abbreviate)
    labs %<>% abbreviate(minlength = 1)

  # draw labels
  text(x$xy[, 1], x$xy[, 2],
       labels = labs, col = col,
       cex = cex, font = font, ...)
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
#' @param rect \code{logical} whether to draw a rectangle below names
layer_labelgroups <- function(x, col=par("fg"), cex=3/4, font=2,
                              rect=TRUE, alpha=1/4, abbreviate=FALSE, ...){
  # if numeric, do not apply this layer but still propagate
  if (!is.factor(x$f))
    return(x)
  cxys <- apply(x$xy, 2, function(.) tapply(., x$f, mean))
  # no f or single level fac
  if (!is.matrix(cxys))
    cxys %<>%
    matrix(nrow=1, ncol=2) %>%
    `row.names<-`(levels(x$f))

  # prepare labs
  labs <- rownames(cxys)
  if (abbreviate)
    labs %<>% abbreviate(minlength = 1)

  # if rect, calculate their dimensions and draw them
  p <- strheight(labs[1])
  if (rect) {
    w <- strwidth(labs, cex = cex)
    h <- strheight(labs, cex = cex)
    rect(xleft   = cxys[, 1] - w/2 - p,
         ybottom = cxys[, 2] - h/2 + p,
         xright  = cxys[, 1] + w/2 + p,
         ytop    = cxys[, 2] + h/2 + p,
         col     = pal_alpha("#FFFFFF", alpha), border = NA)
  }
  # draw labels
  text(cxys[, 1], cxys[, 2] + p,
       labels = labs, col = x$colors_groups,
       cex = cex, font = font, ...)
  # propagate
  invisible(x)
}

# meta layers ----------------------------------------------
#' @param size `numeric` as a fraction of graphical window (default: `1/200`)
#' @export
#' @rdname layers
layer_rug <- function(x, size=1/200, ...){
  # if numeric, apply with a local trick on f
  if (!is.factor(x$f)){
    f <- factor(rep(1, length(x$f)))
    colors_groups <- "#000000"
  } else {
    f <- x$f
    colors_groups <- x$colors_groups
  }
  # neater par
  old <- par(mar=rep(1/8, 4))
  on.exit(par(old))
  # tick size
  h <- size*max(.wdw())
  for (i in seq_along(levels(f))){
    rug(x$xy[f == levels(f)[i], 1], ticksize = h, side=1,
        col=colors_groups[i], lend=1, quiet=TRUE)
    rug(x$xy[f == levels(f)[i], 2], ticksize = h, side=2,
        col=colors_groups[i], lend=1, quiet=TRUE)
  }
  #propagate
  invisible(x)
}

# layers 2 groups ------------------------------------------
#' @param freq `logical`to feed` [hist] (default: `FALSE`)
#' @param breaks to feed [hist] (default: calculated on the pooled values)
#' @param split `logical` whether to split the two distributions into two plots
#' @rdname layers
#' @export
layer_histogram_2 <- function(x, freq=FALSE, breaks, split=FALSE, transp=0){
  # handles par for split or non-split
  if (split){
    op <- par(mfrow=c(2, 1), oma=rep(0, 4), mar=c(4, 4, 3, 1 ))
  } else {
    op <- par(oma=rep(0, 4), mar=c(4, 4, 3, 1 ))
    if (missing(transp))
      transp <- 0.5 # to see overlaps (if any)
  }
  on.exit(op)

  # shortcuts
  xy <- x$xy
  f <- x$f
  xl <- range(xy)*1.2
  cols <- x$colors_groups %>% pal_alpha(transp)

  # homegeneize breaks number
  if (missing(breaks))
    breaks <- graphics::hist(xy, plot=FALSE)$breaks

  # handle ylim for both cases
  hs <- split(xy, f) %>% lapply(graphics::hist, breaks=breaks, plot=FALSE)
  if (freq)
    yl <- c(0, sapply(hs, `[`, "counts") %>% do.call("c", .) %>% max())
  else
    yl <- c(0, sapply(hs, `[`, "density") %>% do.call("c", .) %>% max())

  # draw the two hists
  graphics::hist(xy[f==levels(f)[1], ], freq=freq, breaks=breaks,
                 xlim=xl, ylim=yl, xlab="LD1", ylab=NA,
                 main="",
                 col=cols[1], axes=TRUE)

  if (split)
    title(levels(f)[1])
  else
    layer_legend(x)

  graphics::hist(xy[f==levels(f)[2], ], freq=freq, breaks=breaks,
                 xlim=xl, ylim=yl, xlab="LD1", ylab=NA,
                 main="",
                 col=cols[2], axes=TRUE, add=!split)

  if (split)
    title(levels(f)[2])

  # restore layout (incompatible with on.exit - do not know why)
  graphics::layout(1)
  # propagate
  invisible(x)
}

#' @param bw to feed [density] (default: [stats::bw.nrd0])
#' @param rug `logical` whether to add [rug] (default: `TRUE`)
#' @rdname layers
#' @export
layer_density_2 <- function(x, bw, split=FALSE, rug=TRUE, transp=0){
  # handles par for split or non-split
  if (split){
    op <- par(mfrow=c(2, 1), oma=rep(0, 4), mar=c(4, 4, 3, 1 ))
    ticksize <- 1/20
  } else {
    op <- par(oma=rep(0, 4), mar=c(4, 4, 3, 1 ))
    if (missing(transp))
      transp <- 0.5
    ticksize <- 1/50
  }
  on.exit(op)

  # shortcuts
  xy <- x$xy
  f <- x$f
  cols0 <- x$colors_groups
  cols <- x$colors_groups %>% pal_alpha(transp)

  # define a bandwidth, if missing
  if (missing(bw))
    bw <- stats::bw.nrd(xy)

  # handles xl and yl for both cases
  ds <- split(xy, f) %>% lapply(stats::density, bw=bw)
  xl <- sapply(ds, `[`, "x") %>% do.call("c", .) %>% range()
  yl <- sapply(ds, `[`, "y") %>% do.call("c", .) %>% range()

  # plots the first density
  d1 <- stats::density(xy[f==levels(f)[1], ], bw=bw)
  plot(d1, xlim=xl, ylim=yl, yaxs="i",
       type = "n", xlab="LD1", bty="n", main="")
  polygon(d1, col = cols[1])
  # and the first rug if required
  if (rug)
    graphics::rug(xy[f==levels(f)[1], ], ticksize = ticksize, col=cols0[1], line=1/2)
  # add the 1st title if splitted, legend otherwise
  if (split)
    title(levels(f)[1])
  else
    layer_legend(x)


  # same for second density
  d2 <- stats::density(xy[f==levels(f)[2], ], bw=bw)
  #
  if (split){
    plot(d2,
         xlim=xl, ylim=yl, yaxs="i",
         type = "n", xlab="LD1", bty="n", main="")
    title(levels(f)[2])
  }
  graphics::polygon(d2, col = cols[2])
  # rug if any
  if(rug)
    graphics::rug(xy[f==levels(f)[2], ], ticksize = ticksize, col=cols0[2], line=1/2)

  # restore layout (not working with on.exit)
  graphics::layout(1)

  # propagate
  invisible(x)
}

# cosmetics layers -----------------------------------------
#' @export
#' @rdname layers
#' @param title to add to the plot (default \code{""})
layer_title <- function(x, title="", cex=3/4, ...) {
  # neater par
  old <- par(mar=rep(1/8, 4), xpd=NA)
  on.exit(par(old))
  text(par("usr")[1], par("usr")[4] - strheight(title),
       pos = 4, labels = title, font = 2, cex = cex, ...)
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
#' @param name to use on axes (default \code{"Axis"})
layer_axesnames <- function(x, cex=3/4, name="Axis", ...){
  # neater par
  old <- par(mar=rep(1/8, 4), xpd=NA)
  on.exit(par(old))
  # calculate dimensions
  gy <- strheight(paste(name, 1), cex = cex)/1.5
  gx <- strwidth("00.0%", cex = cex)/1.5
  # add axes names
  text(par("usr")[2] - gx, gy,
       col = "grey50", cex = cex,
       labels = paste(name, x$axes[1]), ...)
  text(-gy, par("usr")[4] - gx,
       col = "grey50", cex = cex,
       labels = paste(name, x$axes[2]), srt = 90, ...)
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
#' @param nb_max \code{numeric} number of eigen values to display (default \code{5})
layer_eigen <- function(x, nb_max=5, cex=1/2, ...){
  # if non pertinent, propagate
  if (is.null(x$sdev))
    invisible(x)
  # neater par
  old <- par(mar=rep(0, 4), xpd=NA)
  on.exit(par(old))

  # default dimensions
  range_padding = 1/50
  range_width   = 1/20
  range_height  = 1/10
  # window dimensions
  u <- par("usr")
  w <- .wdw()
  # deduce subplot bounding box
  x0 <- u[2] - w[1]*(range_width+range_padding)
  x1 <- u[2] - w[1]*range_padding
  y0 <- u[3] + w[2]*range_padding
  y1 <- u[3] + w[2]*(range_height+range_padding)
  xs_left  <- seq(x0, x1, length.out = nb_max+1)[-(nb_max+1)]
  xs_right <- xs_left + diff(xs_left[1:2])
  s_ys <- seq(y0, y1, length.out = 3)

  # handle bars selection, height and colors
  if (max(x$axes)>nb_max)
    nb_max <- max(x$axes)
  var <- x$sdev^2
  ev <- var/sum(var)
  h0 <- ev[1:nb_max]
  hb <- h0*range_height*w[2]*2 # *2 because of 0.5 below
  cols <- rep("grey98", nb_max)
  cols[x$axes] <- "grey60"

  # axis 1 ticks
  segments(x0, s_ys, x0 - range_width*w[2]/6, s_ys, lwd=1/2)
  text(x0, s_ys, seq(0, 0.5, length.out = 3), cex=cex, pos=2, adj=1)
  # draw bars
  rect(xs_left, y0, xs_right, y0+hb, col=cols, lwd=1/2)

  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
layer_axesvar <- function(x, cex=3/4, ...){
  # neater par
  old <- par(mar=rep(1/8, 4), xpd=NA)
  on.exit(par(old))
  # calculate variance
  var <- x$sdev^2
  var <- signif(100 * var/sum(var), 3)
  # and dimensions of labels on window
  gx <- strwidth("00.0%", cex = cex)/1.5
  gy <- strheight("00.0%", cex = cex)/1.5
  # add it below the axes
  text(par("usr")[2] - gx, -gy,
       col = "grey50", cex = cex,
       labels = paste0(var[x$axes[1]], "%"), ...)
  text(gy, par("usr")[4] - gx,
       col = "grey50", cex = cex,
       labels = paste0(var[x$axes[2]], "%"), srt = 90, ...)
  # propagate
  invisible(x)
}

#' @export
#' @rdname layers
#' @param probs \code{numeric} sequence to feed \code{stats::quantile}
#' and to indicate where to draw ticks and legend labels
layer_legend <- function(x, probs=seq(0, 1, 0.25), cex=3/4,  ...){
  # neater par
  old <- par(mar=rep(0, 4), xpd=NA)
  on.exit(par(old))
  # default dimensions
  range_padding = 1/30
  range_width   = 1/60
  range_height  = 1/8
  # window dimensions
  u <- par("usr")
  w <- .wdw()
  # factor case
  # simple legend
  if (is.factor(x$f)){
    wid_leg <- levels(x$f) %>% strwidth(cex=cex) %>% max
    x0 <- u[2] - w[1]*(range_padding+range_width) - wid_leg*1.5
    y0 <- u[4] - w[2]*range_padding
    legend(x0, y0, legend = levels(x$f), fill = x$colors_groups,
           border=par("fg"), bty="n", cex = cex)
  }
  # numeric case
  # redraw continuous scale manually
  if (is.numeric(x$f)){
    qf <- stats::quantile(x$f, probs = probs) %>% signif(2)
    wid_leg <- qf %>% strwidth(cex=cex) %>% max()
    # scale position
    s_y0 <- u[4] - w[2]*(range_padding+range_height)
    s_y1 <- u[4] - w[2]*range_padding
    s_x0 <- u[2] - w[1]*(range_padding+range_width) - wid_leg*1.5
    s_x1 <- u[2] - w[1]*range_padding - wid_leg
    # scale_bars positions
    s_ys <- seq(s_y0, s_y1, length.out = 100)
    s_yg <- diff(s_ys[1:2])
    # legend text position
    t_x0 <- s_x1 + w[2]*(range_padding + range_width) - wid_leg*0.9
    t_ys <- s_y0 + probs*(range_height*w[2])
    rect(s_x0, s_ys,
         s_x1, s_ys+s_yg,
         col=x$palette(100), lwd=0)
    rect(s_x0, s_ys[1],
         s_x1, (s_ys+s_yg)[length(s_ys)],
         lwd=1/2, col=NA, border=par("fg"))
    segments(s_x1, t_ys + s_yg*0.5,
             s_x1 + range_width*w[2]/4, t_ys + s_yg*0.5, lwd=1/2)
    text(t_x0, t_ys,
         qf, cex=cex, pos=4)
  }
  # propagate
  invisible(x)
}


#' # add loading vectors
#' .loadings <- function(loadings.mat, d = 1, d.lab = 1.2, col = "red") {
#'   loadings.mat <- loadings.mat * d
#'   loadings.lab <- loadings.mat * d.lab
#'   arrows(0, 0, loadings.mat[, 1], loadings.mat[, 2], angle = 20,
#'          length = 0.1, col = col)
#'   text(loadings.lab[, 1], loadings.lab[, 2], labels = rownames(loadings.lab),
#'        cex = 0.8, col = col)
#' }


