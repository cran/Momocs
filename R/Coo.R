######################################################################
#######################   Coo methods   ##############################
######################################################################

# Domestic ###########################################################
# Builder
Coo <- function(coo, ...){new(Class="Coo",
                              coo=coo,
                              names=names(coo),
                              coo.nb=length(coo),
                              coo.len=as.numeric(lapply(coo, nrow)), ...)}

#Getters
# setMethod(f = "[", signature = "Coo", definition = function(x, i, j, value){
    # switch(EXPR=i,
           # "coo"    = {return(x@coo[[j]])   },
           # "names"  = {return(x@names[j])   },
           # "fac"    = {return(x@fac[j])     },
           # "ldk"    = {return(x@ldk[j])     },
           # "scale"  = {return(x@names[j])   },
           # "coo.nb" = {return(x@coo.nb)     },
           # "coo.len"= {return(x@coo.len[j]) },
           # stop("This attribute does not exist or cannot be accessed"))
  # })

#Setters
# setReplaceMethod(f = "[", signature = "Coo", definition = function(x, i, j, value){
    # switch(EXPR=i,
           # "coo"    = {x@coo[[j]] <- value
                       # x@coo.nb  <- length(x@coo)
                       # x@coo.len <- as.numeric(lapply(x@coo, nrow))},
           # "names"  = {x@names[j] <- value
                       # names(x@coo)[j] <- value},
           # "fac"    = {x@fac   <- value },
           # "ldk"    = {x@ldk   <- value },
           # "scale"  = {x@scale <- value },
           # stop("This attribute does not exist or cannot be modified"))
    # validObject(x)
    # return(x)
  # })

# show ###############################################################
setMethod(f="show", signature="Coo", definition=function(object){
  cat("A Coo object (see ?Coo) \n")
  cat(rep("*", 30),"\n", sep="")
  
  cat("\nGeneral\n")
  cat(rep("-", 20),"\n", sep="")
  cat(" -", object@coo.nb, ifelse(object@coo.nb<2, "outline\n", "outlines\n"))
  cat(" -", round(mean(object@coo.len)), "+/-", round(sd(object@coo.len)), "coordinates per outline\n")
  if (length(object@ldk)!=0) cat(" -", length(object@ldk), "landmarks defined\n") else cat(" - No landmark defined\n")
  if (ncol(object@fac)!=0) cat(" -", ncol(object@fac), "grouping factors defined\n") else cat(" - No groups defined\n")
  
  cat("\nCoordinates: @coo\n")
  cat(rep("-", 20),"\n", sep="")
  print(object@coo[[1]][1:10,])
  cat("...\n")
  
  cat("\nNames: @names\n")
  cat(rep("-", 20),"\n", sep="")
  if(object@coo.nb>20) {
    cat(object@names[1:20])
    cat("\n...") 
  } else {
    cat(object@names)}
  cat("\n")
  
  if (ncol(object@fac)){
    cat("\nFactors: @fac\n")
    cat(rep("-", 20),"\n", sep="")
    for (i in 1:ncol(object@fac)){
      cat(names(object@fac[i]), ":", levels(object@fac[, i]), "\n")
    }}})

# plot ###############################################################
setMethod(f="plot",
          signature="Coo",
          definition=
            function(x, y, method=c("stack", "single", "panel")[1],
                     subset, center=TRUE, scale=TRUE, align=FALSE,
                     col="#70809033", border="#708090", ...) {
              coo <- x@coo
              if (!missing(subset)) coo <- coo[[subset]]
              if (method == "panel") {
                coo.list.panel(coo, cols=col, borders=border)
                return()
              }
              if (center) coo <- lapply(coo, coo.center)
              if (scale)  coo <- lapply(coo, coo.scale)
              if (align)  coo <- lapply(coo, coo.align)
              if (method == "stack") {
                coo.plot(coo[[1]], size=0.5)
                lapply(coo, polygon, col=col, border=border)}
              if (method == "single") {
                for (i in 1:length(coo)) {
                  coo.plot(coo[[i]], col=col, border=border)
                  title(main=names(coo)[i], outer=T, line=-2)
                  readline(prompt = "Press <Enter> to continue...")}}
              return()})

# Methods on coordinates #############################################
setGeneric(name= "align",   def=function(Coo, ...){standardGeneric("align")})
setMethod(f="align", signature="Coo", definition= function(Coo){
    return(Coo(lapply(Coo@coo, coo.align)))})

setGeneric(name= "center",   def=function(Coo, ...){standardGeneric("center")})
setMethod(f="center", signature="Coo", definition=function(Coo){
  return(Coo(lapply(Coo@coo, coo.center)))})

setGeneric(name= "sample",   def=function(Coo, ...){standardGeneric("sample")})
setMethod(f="sample", signature="Coo", definition=function(Coo, nb.pts=100){
  return(Coo(lapply(Coo@coo, coo.sample, nb.pts)))})

setMethod(f="scale", signature="Coo", definition=function(x){
  return(Coo(lapply(x@coo, coo.scale)))})

# Fourier computation ################################################
# compute elliptical Fourier analysis
setMethod(f="eFourier", signature="Coo", definition=
  function(Coo,
           nb.h      = 32,
           smooth.it = 0,
           normalize = TRUE,
           start     = FALSE) {
    coo <- Coo@coo
    col.n <- paste(rep(LETTERS[1:4], each = nb.h), rep(1:nb.h, times = 4), sep="")
    coe <- matrix(ncol = 4 * nb.h, nrow = length(coo), dimnames = list(names(coo), col.n))
    for (i in seq(along = coo)) {
      ef <- efourier(coo[[i]], nb.h = nb.h, smooth.it = smooth.it)
      if (normalize) ef <- efourier.norm(ef, start=start)
      if (ef$A[1] < 0) {
        ef$A <- (-ef$A)
        ef$B <- (-ef$B)
        ef$C <- (-ef$C)
        ef$D <- (-ef$D)
        ef$lnef <- (-ef$lnef)}
      coe[i, ] <- c(ef$A, ef$B, ef$C, ef$D)}
    return(Nef(coe, fac=Coo@fac))})

# harm.qual
setGeneric(name= "harm.qual",   def=function(Coo,
           id        = 1,
           smooth.it = 0,
           harm.range= c(1, 2, 4, 8, 16, 32),
           scale = FALSE,
           center = TRUE,
           align = FALSE,
           method = c("stack", "panel")[1],
           legend = TRUE,
           palette = col.summer,
           shp.col="#70809033",
           shp.border="#708090EE"){standardGeneric("harm.qual")})

setMethod(f="harm.qual", signature="Coo", definition=
  function(Coo,
           id        = 1,
           smooth.it = 0,
           harm.range= c(1, 2, 4, 8, 16, 32),
           scale = FALSE,
           center = TRUE,
           align = FALSE,
           method = c("stack", "panel")[1],
           legend = TRUE,
           palette = col.summer,
           shp.col="#70809033",
           shp.border="#708090EE"){
    # for one signle outline
    coo <- Coo@coo[[id]]
    coo <- coo.smooth(coo, smooth.it)
    if (scale)  coo <- coo.scale(coo)
    if (center) coo <- coo.center(coo)
    if (align)  coo <- coo.align(coo)
    res <- list()
    for (i in seq(along=harm.range)) {
      res[[i]] <- efourier.i(efourier(coo, nb.h=max(harm.range)), nb.h=harm.range[i])}
    # plotting
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    cols <- paste(palette(length(harm.range)), "88", sep="")
    if (method=="stack") {
      coo.plot(coo, main=Coo@names[id], col=shp.col, border=shp.col)
      for (i in seq(along=harm.range)) {lines(res[[i]], col=cols[i], lwd=2)}
      if (legend) {
        legend("topright", legend = as.character(harm.range), bty="o",
               col = cols, lty = 1, lwd=1, bg="#FFFFFFCC", cex=0.7,
               title = "Number of harmonics")}
      if (center) {points(0, 0, pch=3, col=shp.col)}
    } else {
      if (method=="panel") {
        #par(oma=c(1, 1, 3, 0))
        pos <- coo.list.panel(res, col=cols)
        if (legend) {text(x=pos[, 1], y=pos[, 2], as.character(harm.range))}
        mtext(Coo@names[id], side=3, line=0, cex=1.3, outer=TRUE)
      }}})

# harm.quant
setGeneric(name= "harm.quant", def=function(Coo,
           id        = 1:Coo@coo.nb,
           smooth.it = 0,
           harm.range= seq(8, 32, 6),
           scale = FALSE,
           center = TRUE,
           align = FALSE,
           plot = TRUE,
           legend = TRUE,
           palette = col.summer,
           lineat.y=c(1, 5, 10)){standardGeneric("harm.quant")})

setMethod(f="harm.quant", signature="Coo", definition=
  function(Coo,
           id        = 1:Coo@coo.nb,
           smooth.it = 0,
           harm.range= seq(8, 32, 6),
           scale = FALSE,
           center = TRUE,
           align = FALSE,
           plot = TRUE,
           legend = TRUE,
           palette = col.summer,
           lineat.y=c(1, 5, 10)){
    # we prepare ...
    max.h <- max(harm.range)
    nb.pts <- max.h*2
    cols <- palette(length(harm.range))
    coo.nb <- length(id) 
    RES <- array(dim=c(length(harm.range), nb.pts, coo.nb))
    for (j in seq(along=id)) {
      # we prepare         
      coo <- Coo@coo[[id[j]]]
      coo <- coo.smooth(coo, smooth.it)
      if (scale)  coo <- coo.scale(coo)
      if (center) coo <- coo.center(coo)
      if (align)  coo <- coo.align(coo)
      ef.max  <- efourier(coo, nb.h=max.h)
      coo.max <- l2m(efourier.i(ef.max, nb.h=max.h, nb.pts))
      # we prepare and fill the matrix of deviations for every harmonic
      res <- matrix(nrow=length(harm.range), ncol=nb.pts)
      for (i in seq(along=harm.range)) {
        coo.h    <- l2m(efourier.i(ef.max, harm.range[i], nb.pts))
        res[i, ] <- edm(coo.h, coo.max)}
      RES[,,j] <- res  
    }
    dev.med <- apply(RES, 1:2, median)
    dev.se  <- apply(RES, 1:2, sd)/sqrt(coo.nb)
    if (plot) {
      # for St Thomas
      plot(NA, xlim=c(1, nb.pts),   xlab="Points sampled along the outline",
           ylim=c(0, max(res)), ylab="Deviation in pixels",
           main="Deviations along the outline", yaxs="i", xaxs="i", las=1)
      abline(h=lineat.y, lty=2, col="grey80")
      if (coo.nb > 1) {
        dev.plot(dev.med, dev.se, cols=cols)
      } else {
        dev.plot(dev.med, cols=cols)}
      if (legend) {
        legend("topright", legend = as.character(harm.range), bty="o",
               col = cols, lty = 1, lwd=1, bg="#FFFFFFCC", cex=0.7,
               title = "Number of harmonics")}
    }
    rownames(dev.se) <- rownames(dev.med) <- paste("harm", harm.range, sep="")
    colnames(dev.se) <- colnames(dev.med)  <- paste("coo", 1:nb.pts, sep="")
    if (coo.nb > 1) {
      return(list(dev.med=dev.med, dev.se=dev.se))
    } else {
      return(list(dev.med=dev.med)) }
  })

# harm.pow
setGeneric(name= "harm.pow", def=function(Coo,
           id=1:Coo@coo.nb,
           probs=c(0, 0.5, 1),
           nb.h = 24,
           drop   = 1,
           smooth.it = 0,
           plot = TRUE,
           legend = TRUE,
           title="Fourier power spectrum",
           lineat.x=seq(0, nb.h, by=6),
           lineat.y=c(0.9, 0.99),
           bw=0.1){standardGeneric("harm.pow")})
setMethod(f="harm.pow", signature="Coo", definition=
  function(Coo,
           id=1:Coo@coo.nb,
           probs=c(0, 0.5, 1),
           nb.h = 24,
           drop   = 1,
           smooth.it = 0,
           plot = TRUE,
           legend = TRUE,
           title="Fourier power spectrum",
           lineat.x=seq(0, nb.h, by=6),
           lineat.y=c(0.9, 0.99),
           bw=0.1) {
    res <- matrix(nrow=length(id), ncol=(nb.h-drop))
    x <- (drop+1) : nb.h
    for (i in seq(along=id)) {
      ef  <- efourier(Coo@coo[[id[i]]], nb.h = nb.h, smooth.it = smooth.it)
      pow <- ((ef$an^2 + ef$bn^2 + ef$cn^2 + ef$dn^2)/2)[x]
      res[i, ] <-  (cumsum(pow)/sum(pow))}
    res <- apply(res, 2, quantile, probs=probs)
    rownames(res) <- paste("q", probs)
    colnames(res) <- paste("h", x, sep="")
    if (plot){
      if (length(probs)!=3) stop("probs of length != 3 are not (yet) supported")  
      plot(NA, xlim = range(x), ylim = c(min(res), 1), las=1, yaxs="i", 
           xlab = "Number of harmonics included", ylab = "Cumulative harmonic power",
           main=title, sub=paste(length(id), "outlines"), axes=FALSE)
      axis(1, at=x) ; axis(2)
      abline(v=lineat.x, lty=2, col="grey80")
      abline(h=lineat.y, lty=2, col="grey80")
      segments(x,    res[1, ], x,    res[3, ], lwd=0.5)
      segments(x-bw, res[1, ], x+bw, res[1, ], lwd=0.5)
      segments(x-bw, res[3, ], x+bw, res[3, ], lwd=0.5)
      lines(x, res[2, ], type="o", pch=20, cex=0.6) 
      if (legend) {
        legend("topright", legend = as.character(probs), bty="o", lwd=1,
               bg="#FFFFFFCC", cex=0.7,
               title = "Quantiles")}
      box()
    }
    return(res)})

# Ptolemy
setGeneric(name="Ptolemy", def=function(Coo,
           id=1,
           t=seq(0, 2*pi, length=7)[-1],
           nb.h=3,
           nb.pts=360,
           palette=col.sari,
           legend=FALSE){standardGeneric("Ptolemy")})  
setMethod(f="Ptolemy", signature="Coo", definition=
  function(Coo,
           id=1,
           t=seq(0, 2*pi, length=7)[-1],
           nb.h=3,
           nb.pts=360,
           palette=col.sari,
           legend=FALSE) {    # we prepare and deduce
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(xpd=NA)
    cols <- palette(nb.h)
    coo <- coo.center(Coo@coo[[id]])
    #k <- floor(length(coo$x)/4)
    coo.plot(coo, main=Coo@names[id])
    # now we calculate for every harmonic
    coo.ef  <- efourier(coo, nb.h)
    coo.efi <- efourier.i(coo.ef, nb.h, nb.pts)
    vect   <- matrix(nrow=nb.h, ncol=2)
    vect <- rbind(c(0, 0), vect)
    for (i in seq(along=t)) {
      for(j in 1:nb.h) {
        vect[j+1, 1] <- coo.ef$an[j] * cos(j * t[i]) + coo.ef$bn[j] * sin(j * t[i])
        vect[j+1, 2] <- coo.ef$cn[j] * cos(j * t[i]) + coo.ef$dn[j] * sin(j * t[i])}
      vs <- apply(vect, 2, cumsum)
      for (j in 1:nb.h){
        lh   <- efourier.shape(coo.ef$an[1:j], coo.ef$bn[1:j],
                               coo.ef$cn[1:j], coo.ef$dn[1:j],
                               nb.h=j, nb.pts=nb.pts, plot=FALSE)
        ellh <- efourier.shape(coo.ef$an[j], coo.ef$bn[j],
                               coo.ef$cn[j], coo.ef$dn[j],
                               nb.h=1, nb.pts=nb.pts, plot=FALSE)
        lines(lh, col=paste(cols[j], "22", sep=""), lwd=0.8)
        lines(ellh$x + vs[j, 1], ellh$y + vs[j, 2],
              col=cols[j], lwd=1)
        points(vs[j+1, 1], vs[j+1, 2], col=cols[j], cex=0.8)
        arrows(vs[j, 1], vs[j, 2], vs[j+1, 1], vs[j+1, 2],
               col=cols[j], ang=10, length=0.05, lwd=1.2)
      }
    }
    points(0, 0, pch=20, col=cols[1])
    if (legend) {
      legend("topright", legend = as.character(1:nb.h), bty="o",
             col = cols, lty = 1, lwd=1, bg="#FFFFFFCC", cex=0.7,
             title = "Number of harmonics")}
  })



# smooth.qual
setGeneric(name= "smooth.qual",   def=function(Coo, id=1,
           smooth.range = c(10, 50, 200, 500, 1000),
           palette = col.summer){standardGeneric("smooth.qual")})
setMethod(f="smooth.qual", signature="Coo", definition=
  function(Coo,
           id=1,
           smooth.range = c(10, 50, 200, 500, 1000),
           palette = col.summer){
    cont <- Coo@coo[[id]]
    coo.plot(coo.smooth(cont, 0), main = Coo@names[id])
    cols <- palette(length(smooth.range))
    for (i in seq(along = smooth.range)) {
      lines(coo.smooth(cont, smooth.range[i]), col = cols[i])
    }
    legend("topright", legend = as.character(smooth.range), 
           col = cols, lty = 1, bty = "o", bg="#FFFFFFCC", cex=0.7, title = "Smooth iterations")
  })

