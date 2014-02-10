#--------------------------------------------------------------------#
# Coe methods                                                        #
#--------------------------------------------------------------------#

# Coe declaration and builder ########################################
setClass(
  Class = "Coe",
  representation(coe  = "matrix",
                 names  = "character",
                 fac    = "data.frame",
                 nb.h   = "numeric",
                 method = "character"),
  validity=function(object){
    if(!is.matrix(object@coe)) cat("A matrix must be provided\n")
  }
)

# Coe builder ########################################################
Coe <- function(coe, names, method=c("eFourier", "rFourier", "tFourier"),...){
  if (missing(method)) {
    stop("method must be specified")    
    }
  # we allow partial matching of method argument
  methods <- c("eFourier", "rFourier", "tFourier")
  p <- pmatch(method, methods)
  if (is.na(p)) {
    stop("A valid method name must be passed to 'Coe' builder")
    } else {method <- methods[p]}
  if (missing(names)) {
    names <- if (is.null(rownames(coe))) {
                  names <- paste0("shp", 1:nrow(coe))
                  rownames(coe) <- names
                  names
              } else {rownames(coe)}
    }
  nb.h <- if (method == "eFourier") {ncol(coe)/4} else {ncol(coe)/2}
  # we want pretty names matrices
  if (is.null(colnames(coe))) {
    if (method == "eFourier") {
      colnames(coe) <- paste0(rep(LETTERS[1:4], each=nb.h), rep(1:nb.h, times=4))
    } else {
      colnames(coe) <- paste0(rep(LETTERS[1:2], each=nb.h), rep(1:nb.h, times=2))
    }
  }  
  new(Class="Coe",
      coe=coe,
      names=names,
      nb.h=nb.h,
      method=method,
      ...)}

# Coe getters ########################################################
# setMethod(
  # f = "[",
  # signature = "Coe",
  # definition = function(x, i, j, value){
    # switch(EXPR=i,
           # "coe" = {return(x@coe[j]) },
           # "names" = {return(x@names[j]) },
           # "fac"   = {return(x@fac[j])   },
           # "nb.h"  = {return(x@nb.h)     },
           # stop("This attribute does not exist or cannot be accessed"))})

# Coe setters ########################################################
# setReplaceMethod(
  # f = "[",
  # signature = "Coe",
  # definition = function(x, i, j, value){
    # switch(EXPR=i,
           # "coe"  = {x@coe <- value
                       # x@nb.h  <- ncol(x@coe)/4},
           # "names"  = {x@names[j] <- value
                       # rownames(x@coe)[j] <- value},
           # "fac"    = {x@fac   <- value },
           # stop("This attribute does not exist or cannot be modified"))
    # validObject(x)
    # return(x)})

# show(Coe) ##########################################################
setMethod(f="show", signature="Coe", definition=function(object){

  p <- pmatch(object@method, c("eFourier", "rFourier", "tFourier"))
  met <- switch(p, "elliptical Fourier", "radii variation", "tangent angle")
  
  cat("A matrix of harmonic coefficients obtained with", met, "analysis (see ?Coe) \n")
  cat(rep("*", 30),"\n", sep="")
  
  cat("\nGeneral\n")
  cat(rep("-", 20),"\n", sep="")
  cat(" -", nrow(object@coe), ifelse(nrow(object@coe)<2, "outline\n", "outlines\n"))
  cat(" -", object@nb.h, "harmonics\n")
  #if (length(object@ldk)!=0) cat(" -", length(object@ldk), "landmarks defined\n") else cat(" - No landmark defined\n")
  lf <- ncol(object@fac)
  if (lf!=0) cat(" - ", ncol(object@fac), " grouping factor", ifelse(lf>1, "s", ""), " defined\n", sep="") else cat(" - No groups defined\n")
  
  cat("\nSome harmonic coefficients: @coe\n")
  cat(rep("-", 20),"\n", sep="")
  ln <- nrow(object@coe)
  ind <- if (ln > 6) { c(1:3, (ln-2):ln) } else { 1:ln }
  print(signif(object@coe[ind,
                          coeff.sel(retain=ifelse(object@nb.h>3, 3, object@nb.h),
                                    drop=0,
                                    nb.h=object@nb.h,
                                    cph=ifelse(p==1, 4, 2))], 3))
  cat("...\n")
  
  cat("\nNames: @names\n")
  cat(rep("-", 20),"\n", sep="")
  ln <- length(object@names)
  if(ln>10) {
    cat(object@names[1:5], " ...\n", object@names[(ln-5):ln], sep="", fill=60)
  } else {
    cat(object@names, fill=60)}
  cat("\n")
  
  if (ncol(object@fac)){
    cat("\nFactors: @fac\n")
    cat(rep("-", 20),"\n", sep="")
    for (i in 1:ncol(object@fac)){
      cat(names(object@fac[i]), ":", levels(object@fac[, i]), "\n")
    }
  }
})

# boxplot(Coe) #######################################################
setMethod(f="boxplot",
          signature="Coe",
          definition=
            function(x,
                     retain = 8,
                     drop   = 0,
                     vcenter0 = FALSE,
                     palette   = col.gallus,
                     title  = "Variation of harmonic coefficients",
                     legend = TRUE){
              # we deduce and prepare
              nb.h  <- x@nb.h
              p <- pmatch(x@method, c("eFourier", "rFourier", "tFourier"))
              if (is.na(p)) {
                stop("Invalid method in Coe object")
                } else {
                cph <- ifelse(p==1, 4, 2)}
              if (p == 1 & missing(drop)) {drop = 1}
              x     <- x@coe      
              cs    <- coeff.sel(retain=retain, drop=drop, nb.h=nb.h, cph=cph)
              range <- (drop+1):retain
              # we save the old par and prepare the plot
              op <- par(no.readonly = TRUE)
              on.exit(par(op))
              cols <- palette(cph)
              ylim <- range(x[, cs])
              if (vcenter0) {
                mv <- max(abs(range(x[, cs])))
                ylim <- c(-mv, mv) 
                }
              plot(NA, ylim=ylim, xlim=range(range)+c(-0.6,0.6),
                   xlab="Harmonic rank", ylab="Coefficient value", main=title,
                   axes=FALSE, xaxs="i")
              abline(v=range+0.4, col="grey80")
              abline(h=0, col="grey80")
              for (i in 1:cph) {  
                boxplot(x[,(i-1)*nb.h+range],
                        range=0, boxwex=0.2, at=range-0.6 + (i*0.2),
                        col=cols[i], names=FALSE, border=cols[i], axes=FALSE, add=TRUE)
                }
              axis(1, at=range-0.1, labels=range)  
              axis(2)
              if (legend) {
                legend("topright", legend = LETTERS[1:cph], bty="o",
                       fill = cols, bg="#FFFFFFBB",
                       cex=0.7, inset=0.005,
                       title = " Harmonic coefficients ")
              }
              box()
            })

# hist(Coe) ##########################################################
setMethod(f="hist",
          signature="Coe",
          definition=
            function(x,
                     retain = 4,
                     drop   = 0,
                     palette   = col.gallus,
                     title  = "Variation of harmonic coefficients"){
  # we deduce and prepare
  p <- pmatch(x@method, c("eFourier", "rFourier", "tFourier"))
  if (is.na(p)) {
    stop("Invalid method in Coe object")
  } else {
    cph <- ifelse(p==1, 4, 2)}
  cols  <- palette(cph)
  #if (p == 1 & missing(drop)) {drop = 1}
  cs    <- coeff.sel(retain=retain, drop=drop, nb.h=x@nb.h, cph=cph)
  range <- (drop+1):retain
  x     <- x@coe
  # we save the old par and prepare the plot
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(matrix(1:length(cs), ncol=cph, byrow=TRUE))
  par(oma=c(2, 2, 5, 2), mar=rep(2, 4))
  h.names <- paste(rep(LETTERS[1:cph], each=retain-drop), range, sep="")
  cols    <- rep(cols, each=retain-drop)
  for (i in seq(along=cs)) { # thx Dufour, Chessel and Lobry
    h  <- x[, cs[i]] 
    h0 <- seq(min(h), max(h), len=50)
    y0 <- dnorm(h0, mean(h), sd(h))
    hist(h, main=h.names[i], col=cols[i], proba=TRUE, xlab="", ylab="", las=1)
    abline(v=mean(h), lwd=1)
    lines(h0, y0, col = "black", lwd = 2)}
  title(main=title, cex.main=2, font=2, outer=TRUE)
  })

# pca(Coe) ##########################################################
setGeneric(name= "pca",   def=function(Coe,
                                       subset.fac=NULL,
                                       subset.lev=NULL,
                                       row.w,
                                       col.w,
                                       center=TRUE,
                                       scale=FALSE,
                                       scannf=FALSE,
                                       nf=3){standardGeneric("pca")})

setMethod(f="pca",
          signature="Coe",
          definition=
            function(Coe,
                     subset.fac=NULL,
                     subset.lev=NULL,
                     row.w,
                     col.w,
                     center=TRUE,
                     scale=FALSE,
                     scannf=FALSE,
                     nf=5){
              
              if (!is.null(subset.fac) & !is.null(subset.lev)) {
                coe <- Coe@coe[Coe@fac[,subset.fac]==subset.lev, ]
              } else { coe <- Coe@coe }
              
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
              if (!is.null(subset.fac) & !is.null(subset.lev)) {
                dude$fac <- Coe@fac[Coe@fac[,subset.fac]==subset.lev,]
              } else { dude$fac <- Coe@fac }
              dude$coe <- coe
              dude$method <- Coe@method
              return(dude)})

# harm.contrib(Coe) ##################################################
setGeneric(name= "hcontrib",   def=function(Coe,
                     id      = 1,
                     harm.range = 1:floor(Coe@nb.h / 4),
                     amp.h   = c(0, 0.5, 1, 2, 5, 10),
                     palette = col.sari,
                     title   = "Harmonic contribution"){standardGeneric("hcontrib")})

setMethod(f="hcontrib",
          signature="Coe",
          definition=
            function(Coe,
                     id      = 1,
                     harm.range = 1:floor(Coe@nb.h / 4),
                     amp.h   = c(0, 0.5, 1, 2, 5, 10),
                     palette = col.sari,
                     title   = "Harmonic contribution"){
              nb.h <- Coe@nb.h
              mult <- rep(1, Coe@nb.h)
              
              # we handle the method
              p <- pmatch(tolower(Coe@method), c("efourier", "rfourier", "tfourier"))
              if (is.na(p)) { warning("Unvalid method. efourier is used.")
                } else {
                method.i <- switch(p, efourier.i, rfourier.i, tfourier.i)}
              # and prepare a xf accordingly
              if (p==1) {
                an <- Coe@coe[id, 1:nb.h + 0*nb.h]
                bn <- Coe@coe[id, 1:nb.h + 1*nb.h]
                cn <- Coe@coe[id, 1:nb.h + 2*nb.h]
                dn <- Coe@coe[id, 1:nb.h + 3*nb.h]
                xf <- list(an=an, bn=bn, cn=cn, dn=dn)
              } else {
                an <- Coe@coe[id, 1:nb.h + 0*nb.h]
                bn <- Coe@coe[id, 1:nb.h + 1*nb.h]
                xf <- list(an=an, bn=bn)
              }
              # the core below
              res <- list()
              p <- 1 # dirty
              for (j in seq(along=harm.range)){
                for (i in seq(along=amp.h)){
                  mult.loc    <- mult
                  mult.loc[harm.range[j]] <- amp.h[i]
                  xfi <- lapply(xf, function(x) x*mult.loc)
                  res[[p]] <-
                    l2m(method.i(xfi))
                  p <- p+1}}
              
              cols <- rep(palette(length(amp.h)), length(harm.range))
              coo.list.panel(res, dim=c(length(amp.h), length(harm.range)),
                             byrow=FALSE, cols=cols, mar=c(5.1, 5.1, 4.1, 2.1))
              axis(1, at=(1:length(harm.range))-0.5,
                   labels=harm.range, line=2, lwd=1, lwd.ticks=0.5)
              mtext("Harmonic rank", side=1, line=4)
              axis(2, at=(1:length(amp.h))-0.5,
                   labels=rev(amp.h), line=1, lwd=1, lwd.ticks=0.5)
              mtext("Amplification factor", side=2, line=3)
              title(main=title)
            })

# ellipse.par ########################################################
setGeneric(name= "ellipse.par",   def=function(Coe,
                     range  = 1:nrow(Coe@coe),
                     nb.pts = 120){standardGeneric("ellipse.par")})
setMethod(f="ellipse.par",
          signature="Coe",
          definition=
            function(Coe,
                     range  = 1:nrow(Coe@coe),
                     nb.pts = 120){
  if (Coe@method != "eFourier") stop("ellipse.par can only be calculated when Coe are obtained with eFourier")
  nb.h  <- Coe@nb.h
  # we prepare the result matrices
  res.p <- res.e <- res.b  <- res.a <- matrix(nrow=length(range), ncol=nb.h)
  rownames(res.p) <- rownames(res.e) <- rownames(res.b) <- rownames(res.a) <- rownames(Coe@coe)[range]
  colnames(res.a) <- paste("a",   1:nb.h, sep="") 
  colnames(res.b) <- paste("b",   1:nb.h, sep="")
  colnames(res.e) <- paste("e",   1:nb.h, sep="")
  colnames(res.p) <- paste("phi", 1:nb.h, sep="")
  if (length(range) > 10) {
    pb <- txtProgressBar(1, length(range))
    t <- TRUE } else {t <- FALSE}
  # we calculate ellipse parameters
  for (i in seq(along=range)){
    coe <- Coe@coe[range[i], ]
    for (j in 1:nb.h) {
      a <- coe[j + 0*nb.h]
      b <- coe[j + 1*nb.h]
      c <- coe[j + 2*nb.h]
      d <- coe[j + 3*nb.h]
      # a,b and e
      ell <- efourier.shape(a, b, c, d, nb.pts=nb.pts, plot=FALSE)
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
    if (t) setTxtProgressBar(pb, j)
  }
  return(list(a=res.a, b=res.b, e=res.e, phi=res.p))})

# Multivariate analysis of variance ##################################
setGeneric(name= "manova.Coe",def=function(
  Coe,
  fac,
  retain, drop=0){standardGeneric("manova.Coe")})

setMethod(f="manova.Coe", signature="Coe", definition=
  function (Coe,
            fac,
            retain, drop=0){
    p <- pmatch(Coe@method, c("eFourier", "rFourier", "tFourier"))
    if (is.na(p)) {
      stop("Invalid method in Coe object")
    } else {
      cph <- ifelse(p==1, 4, 2)}
    fr   <- floor((nrow(Coe@coe)-2)/cph) #full rank
    nb.h <- Coe@nb.h
    # dirty below. really needs clarification
    if (missing(fac)) stop("fac must be provided")
    if (!is.null(Coe@fac[,fac])) {
      fac <- Coe@fac[, fac]}
    if (!missing(retain)) {
      if ((retain - drop) > fr) {
        retain <- fr
        if (retain > nb.h) retain <- nb.h
        cat("The number or retained harmonics was too high. Analysis done with", retain, "harmonics\n")}
    } else {
      retain <- fr
      if (retain > nb.h) retain <- nb.h
      cat("The number or retained harmonics was not specified. Analysis done with", retain, "harmonics\n")
    }
    harm.sel <- coeff.sel(retain=retain, drop=drop, nb.h=Coe@nb.h, cph=cph)
    if(length(fac)!=nrow(Coe@coe)) {
      stop("The length of the factor provided (", length(fac),") provided is different than the
      number of outlines (", nrow(Coe@coe), ") in the Coe object provided")}
    mod <- summary(manova(Coe@coe[,harm.sel]~fac), test="Hotelling")
    return(mod)})

setGeneric(name= "meanShapes",def=function(
  Coe,
  fac,
  nb.pts=300){standardGeneric("meanShapes")})

setMethod(f="meanShapes", signature="Coe", definition=
  function (Coe,
            fac,
            nb.pts=300){
    # we handle the method
    p <- pmatch(tolower(Coe@method), c("efourier", "rfourier", "tfourier"))
    if (is.na(p)) { warning("Unvalid method. efourier is used.")
    } else {
      method.i <- switch(p, efourier.i, rfourier.i, tfourier.i)
      cph      <- ifelse(p==1, 4, 2)}
    fr   <- floor((nrow(Coe@coe)-2)/cph) #full rank
    nb.h <- Coe@nb.h
    # dirty below. really needs clarification
    if (ncol(Coe@fac)==1) { fac = 1 }
    if (missing(fac)) {
      cs  <- apply(Coe@coe, 2, mean)
      xf  <- coeff.split(cs=cs, nb.h=Coe@nb.h, cph=cph)
      res <- l2m(method.i(xf, nb.pts=nb.pts)) 
      warning("No fac provided, the mean shape is returned.")
      return(res)}
    
    if (!is.null(Coe@fac[,fac])) {
      fac <- Coe@fac[, fac]}
    f  <- factor(fac)
    fl <-  levels(f)
    
    res <- list()
    for (i in seq(along=fl)){
      cs <- apply(Coe@coe[f==fl[i], ], 2, mean)
      xf <- coeff.split(cs=cs, nb.h=Coe@nb.h, cph=cph)
      res[[i]] <- l2m(method.i(xf, nb.pts=nb.pts))
    }
    names(res) <- fl
    return(res)})

# a preliminary version of clust
setGeneric(name= "clust",def=function(
  Coe,
  fac,
  method = "euclidean",
  type="fan",
  palette=col.summer) standardGeneric("clust"))

setMethod(f="clust", signature="Coe", definition=
            function (Coe,
                      fac,
                      method = "euclidean",
                      type="fan",
                      palette=col.summer){
              if (missing(fac)) {
                if (ncol(Coe@fac)==1) {
                  fac <- 1
                  facs <- Coe@fac[, fac]
                  cols <- palette(nlevels(facs))[facs]
                } else {
                  cols <- rep("black", nrow(Coe@coe))}
              } else {
                facs <- Coe@fac[, fac]
                cols <- palette(nlevels(facs))[facs]
              }
              
              Coe.hc <- hclust(dist(Coe@coe, method=method))
              op <- par(no.readonly = TRUE)
              par(oma=rep(0, 4), mar=rep(0,4))
              plot(as.phylo.hclust(Coe.hc), tip.color=cols, type=type)
              par(op)
              return(Coe.hc)})

# end of Coe.R 

