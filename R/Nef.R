######################################################################
#######################   Nef methods   ##############################
######################################################################
### Nef class #############################################

# Domestic ###########################################################
Nef <- function(coeff, ...){new(Class="Nef", coeff=coeff,  names=rownames(coeff), nb.h=ncol(coeff)/4, ...)}

## Nef getters
# setMethod(
  # f = "[",
  # signature = "Nef",
  # definition = function(x, i, j, value){
    # switch(EXPR=i,
           # "coeff" = {return(x@coeff[j]) },
           # "names" = {return(x@names[j]) },
           # "fac"   = {return(x@fac[j])   },
           # "nb.h"  = {return(x@nb.h)     },
           # stop("This attribute does not exist or cannot be accessed"))})

## Nef setters
# setReplaceMethod(
  # f = "[",
  # signature = "Nef",
  # definition = function(x, i, j, value){
    # switch(EXPR=i,
           # "coeff"  = {x@coeff <- value
                       # x@nb.h  <- ncol(x@coeff)/4},
           # "names"  = {x@names[j] <- value
                       # rownames(x@coeff)[j] <- value},
           # "fac"    = {x@fac   <- value },
           # stop("This attribute does not exist or cannot be modified"))
    # validObject(x)
    # return(x)})

# show ###############################################################
setMethod(f="show", signature="Nef", definition=function(object){
  cat("A Fourier Elliptical coefficient matrix (see ?Nef) \n")
  cat(rep("*", 30),"\n", sep="")
  
  cat("\nGeneral\n")
  cat(rep("-", 20),"\n", sep="")
  cat(" -", nrow(object@coeff), ifelse(nrow(object@coeff)<2, "outline\n", "outlines\n"))
  cat(" -", object@nb.h, "harmonics\n")
  #if (length(object@ldk)!=0) cat(" -", length(object@ldk), "landmarks defined\n") else cat(" - No landmark defined\n")
  if (ncol(object@fac)!=0) cat(" -", ncol(object@fac), "grouping factors defined\n") else cat(" - No groups defined\n")
  
  cat("\nHarmonic coefficients (only the first individuals and harmonics are shown): @coeff\n")
  cat(rep("-", 20),"\n", sep="")
  nb.ind <- if(nrow(object@coeff) > 6) {1:6} else {1:nrow(object@coeff)}
  print(signif(object@coeff[nb.ind, coeff.sel(ceiling(object@nb.h/6), nb.h=object@nb.h)], 3))
  cat("...\n")
  
  cat("\nNames: @names\n")
  cat(rep("-", 20),"\n", sep="")
  if(length(object@names)>20) {
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
    }
  }
})

# boxplot ############################################################
setMethod(f="boxplot",
          signature="Nef",
          definition=
            function(x,
                     retain = 8,
                     drop   = 1,
                     palette   = col.gallus,
                     title  = "Variation of harmonic coefficients",
                     legend = TRUE){
              # we deduce and prepare
              nb.h <- x@nb.h
              x     <- x@coeff
              cs    <- coeff.sel(retain=retain, drop=drop, nb.h=nb.h)
              range <- (drop+1):retain
              # we save the old par and prepare the plot
              op <- par(no.readonly = TRUE)
              on.exit(par(op))
              cols <- palette(4)
              plot(NA, ylim=range(x[, cs]), xlim=range(range)+c(-0.6,0.6),
                   xlab="Harmonic number", ylab="Coefficient value", main=title,
                   axes=FALSE, xaxs="i")
              abline(v=range+0.4, col="grey80")
              abline(h=0, col="grey80")
              for (i in 1:4) {  
                boxplot(x[,(i-1)*nb.h+range],
                        range=0, boxwex=0.2, at=range-0.6 + (i*0.2),
                        col=cols[i], names=FALSE, border=cols[i], axes=FALSE, add=TRUE)}
              axis(1, at=range-0.1, labels=range)  
              axis(2)
              if (legend) {
                legend("topright", legend = LETTERS[1:4], bty="o",
                       fill = cols, bg="#FFFFFFCC", cex=0.7, title = "Harmonic coefficients")
              }
              box()
            })

# hist ###############################################################
setMethod(f="hist",
          signature="Nef",
          definition=
            function(x,
                     retain = 4,
                     drop   = 0,
                     palette   = col.gallus,
                     title  = "Variation of harmonic coefficients"){
              # we deduce and prepare
              cols <- palette(4)
              x     <- x@coeff
              cs    <- coeff.sel(retain=retain, drop=drop, nb.h=ncol(x)/4)
              range <- (drop+1):retain
              # we save the old par and prepare the plot
              op <- par(no.readonly = TRUE)
              on.exit(par(op))
              layout(matrix(1:length(cs), ncol=4, byrow=TRUE))
              par(oma=c(2, 2, 5, 2), mar=rep(2, 4))
              h.names <- paste(rep(LETTERS[1:4], each=retain-drop), range, sep="")
              cols    <- rep(cols, each=retain-drop)
              for (i in seq(along=cs)) { # thx Dufour, Chessel and Lobry
                h  <- x[, cs[i]] 
                h0 <- seq(min(h), max(h), len=50)
                y0 <- dnorm(h0, mean(h), sd(h))
                hist(h, main=h.names[i], col=cols[i], proba=TRUE, xlab="", ylab="")
                abline(v=mean(h), lwd=1)
                lines(h0, y0, col = "black", lwd = 2)}
              title(main=title, cex.main=2, font=2, outer=TRUE)
            })

# ade4 ###############################################################
setGeneric(name= "pca",   def=function(Nef,
                     subset=NULL,
                     row.w,
                     col.w,
                     center=TRUE,
                     scale=FALSE,
                     scannf=FALSE,
                     nf=3){standardGeneric("pca")})
setMethod(f="pca",
          signature="Nef",
          definition=
            function(Nef,
                     subset=NULL,
                     row.w,
                     col.w,
                     center=TRUE,
                     scale=FALSE,
                     scannf=FALSE,
                     nf=3){
              
              if (!is.null(subset)) {
                coe <- Nef@coeff[subset, ]
              } else { coe <- Nef@coeff }
              
              if (missing(row.w)) row.w <- rep(1, nrow(coe))/nrow(coe)
              if (missing(col.w)) col.w <- rep(1, ncol(coe))
               
              dude <- dudi.pca(df=as.data.frame(coe),
                               row.w = row.w,
                               col.w = col.w,
                               center= center,
                               scale = scale,
                               scannf = scannf,
                               nf = nf)
              dude$mean.shp <- apply(coe, 2, mean)
              if (!is.null(subset)) {
                dude$fac <- Nef@fac[subset, ] 
              } else { dude$fac <- Nef@fac }
              dude$coe <- coe
              return(dude)})

# harm.contrib #######################################################
setGeneric(name= "harm.contrib",   def=function(Nef,
                     id      = 1,
                     range.h = 1:floor(Nef@nb.h / 4),
                     amp.h   = c(0, 0.5, 1, 2, 5, 10),
                     palette = col.sari,
                     title   = "Harmonic contribution"){standardGeneric("harm.contrib")})
setMethod(f="harm.contrib",
          signature="Nef",
          definition=
            function(Nef,
                     id      = 1,
                     range.h = 1:floor(Nef@nb.h / 4),
                     amp.h   = c(0, 0.5, 1, 2, 5, 10),
                     palette = col.sari,
                     title   = "Harmonic contribution"){
              nb.h <- Nef@nb.h
              mult <- rep(1, Nef@nb.h)
              an <- Nef@coeff[id, 1:nb.h + 0*nb.h]
              bn <- Nef@coeff[id, 1:nb.h + 1*nb.h]
              cn <- Nef@coeff[id, 1:nb.h + 2*nb.h]
              dn <- Nef@coeff[id, 1:nb.h + 3*nb.h]
              res <- list()
              p <- 1 # bloody dirty
              for (j in seq(along=range.h)){
                for (i in seq(along=amp.h)){
                  mult.loc    <- mult
                  mult.loc[range.h[j]] <- amp.h[i]
                  res[[p]] <-
                    l2m(efourier.shape(an*mult.loc, bn*mult.loc,
                                       cn*mult.loc, dn*mult.loc,
                                       nb.h=nb.h, plot=FALSE))
                  p <- p+1
                }}
              cols <- rep(palette(length(amp.h)), length(range.h))
              coo.list.panel(res, dim=c(length(amp.h), length(range.h)),
                             byrow=FALSE, col=cols, mar=c(5.1, 5.1, 4.1, 2.1))
              axis(1, at=(1:length(range.h))-0.5, labels=range.h, line=2, lwd=1, lwd.ticks=0.5)
              mtext("Harmonic number", side=1, line=4)
              axis(2, at=(1:length(amp.h))-0.5,   labels=rev(amp.h), line=1, lwd=1, lwd.ticks=0.5)
              mtext("Amplification factor", side=2, line=3)
              title(main=title)
            })

# advancing the frontiers ############################################

# ellipses ###########################################################
# ellipse.par 
setGeneric(name= "ellipse.par",   def=function(Nef,
                     range  = 1:nrow(Nef@coeff),
                     nb.pts = 120){standardGeneric("ellipse.par")})
setMethod(f="ellipse.par",
          signature="Nef",
          definition=
            function(Nef,
                     range  = 1:nrow(Nef@coeff),
                     nb.pts = 120){
              nb.h  <- Nef@nb.h
              # we prepare the result matrices
              res.p <- res.e <- res.b  <- res.a <- matrix(nr=length(range), nc=nb.h)
              rownames(res.p) <- rownames(res.e) <- rownames(res.b) <- rownames(res.a) <- rownames(Nef@coeff)[range]
              colnames(res.a) <- paste("a",   1:nb.h, sep="") 
              colnames(res.b) <- paste("b",   1:nb.h, sep="")
              colnames(res.e) <- paste("e",   1:nb.h, sep="")
              colnames(res.p) <- paste("phi", 1:nb.h, sep="")
              # we calculate ellipse parameters
              for (i in seq(along=range)){
                coe <- Nef@coeff[range[i], ]
                for (j in 1:nb.h) {
                  a <- coe[j + 0*nb.h]
                  b <- coe[j + 1*nb.h]
                  c <- coe[j + 2*nb.h]
                  d <- coe[j + 3*nb.h]
                  # a,b and e
                  ell <- efourier.shape(a, b, c, d, nb.pts=nb.pts)
                  ep  <- ellpar(ell)
                  res.a[i, j] <- ep$a
                  res.b[i, j] <- ep$b
                  res.e[i, j] <- ep$e
                  # phi - more tricky
                  theta      <- (0.5 * atan(2 * (a*b + c*d)/(a^2 + c^2 - b^2 - d^2))) %% pi           
                  phaseshift <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
                  M2         <- matrix(c(a, c, b, d), 2, 2) %*% phaseshift
                  v          <- apply(M2^2, 2, sum)
                  if (v[1] < v[2]) {theta <- theta + pi/2}
                  theta       <- (theta + pi/2)%%pi - pi/2
                  a2          <- a*cos(theta) + b*sin(theta)
                  c2          <- c*cos(theta) + d*sin(theta)
                  res.p[i, j] <- atan(c2/a2) 
                }
              }
              return(list(a=res.a, b=res.b, e=res.e, phi=res.p))})

# Multivariate analysis of variance ##################################
setGeneric(name= "manova.nef",def=function(
  Nef,
  fac,
  harmonics.retained,
  retain=8, drop=0){standardGeneric("manova.nef")})
setMethod(f="manova.nef", signature="Nef", definition=
  function (Nef,
            fac,
            harmonics.retained,
            retain=8, drop=0){
    if (missing(fac)) stop("fac must be provided")
    if (!is.null(Nef@fac[,fac])) {
      fac <- Nef@fac[, fac]}
    if (missing(retain)) {
      retain <- floor((nrow(Nef@coeff)/4)-2)
    } else if(harmonics.retained>(4*floor((nrow(Nef@coeff)/4)-2))) {
      cat("The number of harmonics to retained is too high\n")
      harmonics.retained <- floor((nrow(Nef@coeff)/4)-2)}
    if(length(fac)!=nrow(Nef@coeff)) {
      stop("The length of the factor provided (", length(fac),") provided is different than the
      number of outlines (", nrow(Nef@coeff), ") in the Nef object provided")}
  harm.sel <- coeff.sel(retain, drop, ncol(Nef@coeff)/4)
    print(summary(manova(Nef@coeff[,harm.sel]~fac), test="Hotelling"))
    cat(retain, "harmonics are retained\n")})

###### Thin Plate Splines
# Multivariate analysis of variance ##################################
setGeneric(name= "tps",def=function(
  Nef,
  id.fr      = 1,
  id.to      = 5,
  nb.pts     = 64,
  amp        = 1,
  method     = c("grid", "arr", "iso")[1],
  grid.size  = 20,
  grid.col   = "grey50",
  
  arr.nb     = 100,
  arr.levels = 100 , 
  arr.len    = 0.1,
  arr.ang    = 30,
  arr.lwd    = 1,
  arr.col    = "grey50",
  
  
  iso.levels = 10, 
  iso.nb     = 5000,
  
  palette    = col.summer,
  cont       = TRUE,
  cont.col   = col.gallus(2),
  cont.lwd   = c(2, 2)){standardGeneric("tps")})
  
setMethod(f="tps", signature="Nef", definition=
  function (Nef,
                id.fr      = 1,
                id.to      = 5,
                nb.pts     = 64,
                amp        = 1,
                method     = c("grid", "arr", "iso")[1],
                grid.size  = 20,
                grid.col   = "grey50",
                
                arr.nb     = 100,
                arr.levels = 100 , 
                arr.len    = 0.1,
                arr.ang    = 30,
                arr.lwd    = 1,
                arr.col    = "grey50",
                
                
                iso.levels = 10, 
                iso.nb     = 5000,
                
                palette    = col.summer,
                cont       = TRUE,
                cont.col   = col.gallus(2),
                cont.lwd   = c(2, 2)){
  
  tps     <- function(matr, matt, n, ret=0, plot=1, col="black"){
    # slightly modified tps{sp}
    # Author: Julien Claude
    xm    <- min(matt[,1])
    ym    <- min(matt[,2])
    xM    <- max(matt[,1])
    yM    <- max(matt[,2])
    rX    <- xM-xm
    rY    <- yM-ym
    a     <- seq(xm-1/5*rX, xM+1/5*rX, length=n)
    b     <- seq(ym-1/5*rX, yM+1/5*rX,by=(xM-xm)*7/(5*(n-1)))
    m     <- round(0.5+(n-1)*(2/5*rX+ yM-ym) / (2/5*rX+ xM-xm))
    M     <- as.matrix(expand.grid(a,b))
    ngrid <- tps2d(M,matr,matt)
    if (plot) {
      plot(ngrid, cex=0.2, asp=1, axes=F, ann=FALSE, mar=rep(0,4))
      for (i in 1:m) lines(ngrid[(1:n) + (i-1)*n,], col=col)
      for (i in 1:n) lines(ngrid[(1:m) * n-i+1,],   col=col)
    }
    if (ret) ngrid
  }
  
  tps2d <- function(M, matr, matt) {
    # Interpolate coordinates for TPS deformation
    p  <- dim(matr)[1]
    q  <- dim(M)[1]
    n1 <- p+3
    P  <- matrix(NA, p, p)
    for (i in 1:p) {
      for (j in 1:p) {
        r2     <- sum((matr[i,]-matr[j,])^2)
        P[i,j] <- r2*log(r2)
      }
    }
    P[which(is.na(P))] <- 0
    Q  <- cbind(1, matr)
    L  <- rbind(cbind(P,Q), cbind(t(Q),matrix(0,3,3)))
    m2 <- rbind(matt, matrix(0, 3, 2))
    coefx <- solve(L)%*%m2[,1]
    coefy <- solve(L)%*%m2[,2]
    fx <- function(matr, M, coef) {
      Xn <- numeric(q)
      for (i in 1:q) {
        Z     <- apply((matr-matrix(M[i,],p,2,byrow=T))^2,1,sum)
        Xn[i] <- coef[p+1]+coef[p+2]*M[i,1]+coef[p+3]*M[i,2]+
          sum(coef[1:p]*(Z*log(Z)))
      }
      Xn
    }
    matg     <- matrix(NA, q, 2)
    matg[,1] <- fx(matr, M, coefx)
    matg[,2] <- fx(matr, M, coefy)
    matg
  }
  

  # we retrieve the configuration with PCA-based amplfications
  coe  <- Nef@coeff
  dudi <- pca(Nef)
  mean.shp <- apply(coe, 2, mean)
  fr   <- pca2shp(dudi$li[id.fr, 1:2], rot=dudi$c1[, 1:2], mean.shp=mean.shp, amp=amp, trans=FALSE, nb.pts=nb.pts)[[1]]
  to   <- pca2shp(dudi$li[id.to, 1:2], rot=dudi$c1[, 1:2], mean.shp=mean.shp, amp=amp, trans=FALSE, nb.pts=nb.pts)[[1]]
  
  # simple grid
  if (method=="grid") {
    tps(fr, to, grid.size, col=grid.col) 
  }
  
  # vector field ####
  if(method=="arr") {
    sR    <- spsample(Polygon(rbind(fr, fr[1,])), arr.nb, type="regular")@coords
    sT     <- tps2d(sR, fr, to)
    # grille simple, on affiche d'abord les deux courbes
    wdw      <- apply(rbind(fr, to), 2, range)
    plot(NA, xlim=wdw[, 1], ylim=wdw[, 2], asp=1, axes=FALSE, ann=FALSE, mar=rep(0,4))
    if (missing(arr.levels)) {arr.levels = arr.nb}
    if (!missing(palette)) {
      ed <- function(p1, p2) {return(sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2))}
      q.lev <- numeric(nrow(sR))
      for (i in seq(along=q.lev)) {q.lev[i] <- ed(sR[i,], sT[i,])}
      q.lev   <- cut(q.lev, breaks=arr.levels, labels=FALSE)
      arr.col <- palette(arr.levels)[q.lev]
    } else {
      arr.col <- rep(arr.col, nrow(sR))}
    # on affiche les flèches
    for (i in 1:nrow(sR)) {
      x0 <- sR[i, 1]
      y0 <- sR[i, 2]
      y1 <- sT[i, 2]
      x1 <- sT[i, 1]   
      arrows(x0, y0, x1, y1, 
             length=arr.len, angle=arr.ang, lwd=arr.lwd, col=arr.col[i])}
  }
    
    #iso
  if (method=="iso") {
    sR    <- spsample(Polygon(rbind(fr, fr[1,])), iso.nb, type="regular")@coords
    sT     <- tps2d(sR, fr, to)
    def    <- sqrt(apply((sT-sR)^2,1,sum))
    x1     <- length(unique(sR[,1]))
    y1     <- length(unique(sR[,2]))
    im     <- matrix(NA,x1,y1)
    xind   <- (1:x1)[as.factor(rank(sR[,1]))]
    yind   <- (1:y1)[as.factor(rank(sR[,2]))]
    n      <- length(xind)
    for (i in 1:n) im[xind[i], yind[i]] <- def[i]
    colors <- palette(iso.levels)
    x <- sort(unique(sR[,1]))
    y <- sort(unique(sR[,2]))
    image(x, y, im, col=colors, asp=1, axes=F, frame=F, ann=F, ylim=range(y)*1.1, xlim=range(x)*1.1)
    contour(x, y, im, nlevels=iso.levels, add=TRUE)
  }
    
    if (cont) {
    # common, gestion contours #####
    lines(coo.close(fr), lwd=cont.lwd[1], col=cont.col[1])
    lines(coo.close(to), lwd=cont.lwd[2], col=cont.col[2])
    points(fr[1, 1], fr[1, 2], col=cont.col[1])
    points(to[1, 1], to[1, 2], col=cont.col[2])
  }
}    )

    
    
	