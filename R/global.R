######################################################################
#######################   Coo methods   ##############################
######################################################################

# Classes stuff ######################################################
# Class declaration
# Coo definition
setClass(
  Class = "Coo",
  representation(coo       = "list",
                 names     = "character",
                 fac       = "data.frame",
                 ldk       = "list",
                 scale     = "numeric",
                 coo.nb    = "numeric",
                 coo.len   = "numeric",
                 details   = "list"),
  validity=function(object){
    if(any(lapply(object@coo, class) != "matrix")) cat("A list of coordinates as 2-col matrices must be provided\n")
    if(any(lapply(object@coo, ncol)  != 2))        cat("A list of coordinates as 2-col matrices must be provided\n")
  }
)

# Nef definition
setClass(
  Class = "Nef",
  representation(coeff  = "matrix",
                 names  = "character",
                 fac    = "data.frame",
                 nb.h   = "numeric"),
  validity=function(object){
    if(!is.matrix(object@coeff)) cat("A matrix must be provided\n")
  }
)

# Generic declaration
### Coo
setGeneric(name= "eFourier", 	def=function(Coo,
           nb.h      = 32,
           smooth.it = 0,
           normalize = TRUE,
           start     = FALSE){standardGeneric("eFourier")})
### Nef

setGeneric(name= "pca", 		  def=function(Nef, ...){standardGeneric("pca")})


# Coordinates utils ##################################################
coo.center     <- function(coo){
  if (is.matrix(coo)) {  
    return(apply(coo, 2, function(x) x-mean(x)))
  } else if (is.list(coo)){
    return(lapply(coo, function(x) x-mean(x)))
  } else stop("A list or a coordinate matrix must be provided to coo.center")}

coo.scale <- 
  function (coo, scale = 1) 
  {
    if (is.matrix(coo)) {
      r <- apply(coo, 2, function(x) diff(range(x)))
      r <- ifelse(r[1] > r[2], r[1], r[2])
      lm <- matrix(c(scale/r, 0, 0, scale/r), ncol=2)
      return(coo %*% lm)
    }
    else if (is.list(coo)) {
      r <- apply(coo, 2, function(x) diff(range(x)))
      r <- ifelse(r[1] > r[2], r[1], r[2])
      lm <- matrix(c(scale/r, 0, 0, scale/r), ncol=2)
      return(m2l(coo %*% lm))
    }
    else stop("A list or a coordinate matrix must be provided to coo.scale")
  }


coo.rotate     <- function(coo, theta){
  rmat <- matrix(c(cos(theta), sin(theta),
                   -sin(theta), cos(theta)), nrow=2)
  if (is.matrix(coo)) {  
    return(coo %*% rmat)
  } else if (is.list(coo)){
    coo <- cbind(coo$x, coo$y)
    coo <- coo %*% rmat
    return(list(x=coo[, 1], y=coo[, 2]))
  } else stop("A list or a coordinate matrix must be provided to coo.rotate")}

coo.align      <- function(coo){
  if (is.matrix(coo)) {  
    return(coo %*% svd(var(coo))$u)
  } else if (is.list(coo)){
    coo <- cbind(coo$x, coo$y)
    coo <- coo %*% svd(var(coo))$u
    return(list(x=coo[, 1], y=coo[, 2]))
  } else stop("A list of a coordinate matrix must be provided to coo.align")}

coo.trans      <- function(coo, x, y){
  if (is.matrix(coo)) {  
    return(cbind(coo[, 1] + x, coo[, 2] + y))
  } else if (is.list(coo)){
    return(list(x=coo$x + x, y=coo$y + y))
  } else stop("A list of a coordinate matrix must be provided to coo.trans")}

coo.slide      <- function(coo, id1){
  if (id1 == 0) {return(coo)}
  if (is.matrix(coo)) {
    n <- nrow(coo)
    id.slided <- c(id1:n, 1:(id1-1))
    return(coo[id.slided, ])
  } else if (is.list(coo)){
    n <- length(coo$x)
    id.slided <- c(id1:n, 1:(id1-1))
    return(list(x=coo$x[id.slided], y=coo$y[id.slided]))
  } else stop("A list of a coordinate matrix must be provided to coo.slide")}

coo.sample     <- function (coo, n) {
  if (is.matrix(coo)){
    sampled <- round(seq(1, nrow(coo), len = n + 1)[-(n + 1)])
    return(coo[sampled, ])
  } else if (is.list(coo)) {
    sampled <- round(seq(1, length(coo$x), len = n + 1)[-(n + 1)])
    return(list(x=coo$x[sampled], y=coo$y[sampled]))  
  } else stop("A list of a coordinate matrix must be provided to coo.sample")}

coo.sample.rr  <- function(coo, n){
 if (is.matrix(coo)) {
   Rx <- coo[, 1]
   Ry <- coo[, 2] }
  if (is.list(coo)){
   Rx <- coo$x
   Ry <- coo$y }
  le  <-length(Rx)
  M   <-matrix(c(Rx, Ry), le, 2)
  M1  <-matrix(c(Rx-mean(Rx), Ry-mean(Ry)), le, 2)
  V1  <-complex(real=M1[,1], imaginary=M1[,2])
  M2  <-matrix(c(Arg(V1), Mod(V1)), le, 2)
  V2  <-NA
  for (i in 0:(n-1)){
    V2[i+1]<-which.max((cos(M2[,1]-2*i*pi/n)))}
  V2<-sort(V2)
  list("pixindices"=V2,"radii"=M2[V2,2],"phase"=M2[V2,1],
      "coord"=M1[V2,], "orig.coord"=matrix(c(M1[,1]+mean(Rx), M1[,2]+mean(Ry)), ncol=2)[V2,])}

coo.smooth     <- function(coo, n=0){
  if (is.matrix(coo)) {
    p   <- nrow(coo)
    a   <- 0
    while (a < n) {
      a <- a + 1
      coo.i <- rbind(coo[-1, ], coo[1, ])
      coo.s <- rbind(coo[p, ],  coo[-p, ])
      coo   <- coo/2 + coo.i/4 + coo.s/4}
    return(coo)
  } else if (is.list(coo)){
    coo <- cbind(coo$x, coo$y)
    p   <- nrow(coo)
    a   <- 0
    while (a <= n) {
      a <- a + 1
      coo.i <- rbind(coo[-1, ], coo[1, ])
      coo.s <- rbind(coo[p, ],  coo[-p, ])
      coo   <- coo/2 + coo.i/4 + coo.s/4}
    return(list(x = coo[, 1], y = coo[, 2]))
  } else stop("A list of a coordinate matrix must be provided to coo.smooth")}

coo.close      <- function(coo){
  if (is.matrix(coo)) {  
    return(rbind(coo, coo[1, ]))
  } else if (is.list(coo)){
    return(list(x=c(coo$x, coo$x[1]), y=c(coo$y, coo$y[1])))
  } else stop("A list or a coordinate matrix must be provided to coo.center")}

coo.plot       <- function(coo=NA, col="#70809033",
                           border="#708090EE", xlim=c(-1, 1), ylim=c(-1, 1),
                           points=TRUE, first.point=TRUE,
                           points.col=border, pch=20, cex=0.25, main, ...){
  if (is.list(coo)) coo <- l2m(coo)
  if (missing(coo)) {
    plot(NA, asp=1, xlim=xlim, ylim=ylim, ann=FALSE, ...)
  } else {
	if (missing(xlim) & missing(ylim)) {
    plot(coo, type="n", asp=1, xlab="", ylab="", frame=FALSE, ...)
	} else {
	plot(coo, type="n", asp=1, xlim=xlim, ylim=ylim, xlab="", ylab="", frame=FALSE, ...)
	}
    polygon(coo, col=col, border=border)
    if (first.point) {points(coo[1, 1], coo[1, 2], col = border, ...)}
}
  if ((!missing(coo) & points)) {
    points(coo, pch=pch, cex=cex, col=points.col)}
  if (!missing(main)) title(main=main)}  

coo.draw <-
  function (coo = NA, col = "#70809033", border = "#708090EE",
            points = TRUE, first.point=TRUE,
            points.col = border, pch = 20, cex = 0.25,  ...)
  {
    if (is.list(coo)) {coo <- l2m(coo)}
    polygon(coo, col = col, border = border)
    if (first.point) {points(coo[1, 1], coo[1, 2], col = border, ...)}
    if ((missing(points) & nrow(coo) > 100) | points) {
      points(coo, pch = pch, cex = cex, col = points.col)
    }
  }

coo.template   <- function(coo, size=1) {
  # only for matrices
  coo      <- coo * min(size/apply(coo, 2, function(x) diff(range(x))))
  expected <- apply(coo, 2, function(x) diff(range(x)))/2
  observed <- apply(coo, 2, range)[2, ]
  shift    <-  expected - observed
  return(coo.trans(coo, shift[1], shift[2]))}

coo.list.panel <- function(coo.list, dim, byrow=TRUE,
                           fromtop=TRUE, mar=rep(0, 4), cols, borders, density = NULL, angle = 45){
  if (class(coo.list[[1]])=="list") coo.list <- lapply(coo.list, l2m)
  # if dim is missing, we define a square
  n <- length(coo.list)
  if(missing(dim)) {
    nc    <- ceiling(sqrt(n))
    nr    <- ceiling(n/nc)
    dim   <- c(nr, nc)}
  k       <- dim[1]*dim[2]
  if (k < n) stop("dim[1]*dim[2] must be >= the length of coo.list")
  pos <- matrix(1:k, dim[1], dim[2], byrow=byrow)
  if (fromtop) { pos <- pos[dim[1]:1,] }
  # we prepare the panel
  op <- par("mar")
  par(mar=mar)
  plot(NA, asp=1,
       xlim=c(0, dim[2]),
       ylim=c(0, dim[1]),
       xaxs="i", yaxs="i", frame=FALSE, ann=FALSE, axes=FALSE)
  # we template and plot shapes
  coo.tp  <- lapply(coo.list, coo.template, size=0.9)
  if (missing(cols))    cols     <- rep("grey80", n)
  if (missing(borders)) borders  <- rep("grey20", n)
  if (missing(density))    density     <- rep(NULL, n)
  if (missing(angle)) angle  <- rep(45, n)
  res <- data.frame(pos.x=numeric(), pos.y=numeric())
  for (i in 1:n){
    trans <- which(pos==i, arr.ind=TRUE) - 0.5
    res[i, ] <- c(trans[2], trans[1])
    polygon(coo.tp[[i]][, 1] + trans[2],
            coo.tp[[i]][, 2] + trans[1],
            col=cols[i], border=borders[i], density = density[i], angle = angle[i])
  }
  par(mar=op)
  invisible(res)}

coo.oscillo    <- function(coo, method=c("d0", "di")[1], plot=TRUE,
                           rug=TRUE, legend=TRUE, cols=col.gallus(2),
                           ref=FALSE, ref.nb=8, ...){
  if (is.list(coo)) coo <- l2m(coo)
  nr <- nrow(coo)
  if (method=="d0"){
    dx <- coo[, 1] - coo[1, 1]
    dy <- coo[, 2] - coo[1, 2]
  } else {
    if (method=="di"){
      dx <- coo[, 1] - coo[, 1][c(nr, (1:nr - 1))]
      dy <- coo[, 2] - coo[, 2][c(nr, (1:nr - 1))]
    } else {stop("inappropriate method in oscillo")}}
  if (ref & plot) {
    def.par <- par(no.readonly = TRUE)
    on.exit(par(def.par))
    layout(matrix(1:2, ncol=2), widths=c(1, 2))
    refs <- round(seq(1, nrow(coo), length=ref.nb+1)[-(ref.nb+1)])
    coo.plot(coo)
    text(coo[refs, 1], coo[refs, 2], labels=as.character(1:ref.nb), cex=0.7)
    }
  if (plot) {
    ry <- max(abs(range(c(dx, dy))))
    plot(NA, xlim=c(1, nr), xlab="Points sampled along the outline",
         ylim=c(-ry, ry)*1.1,   ylab="Deviation",
         xaxs="i", las=1, col=cols[1])
    lines(dx, col=cols[1], ...)
    lines(dy, col=cols[2], ...)
    abline(h=0, col="grey80", lty=2)
    if (method=="d0" & rug) {
      dx.i <- coo[, 1] - coo[, 1][c(nr, (1:nr - 1))]
      dy.i <- coo[, 2] - coo[, 2][c(nr, (1:nr - 1))]
      axis(3, at=(1:nr)[dx.i>0], line=-0.5, labels=FALSE, col=cols[1])
      axis(3, at=(1:nr)[dy.i>0], line=-1, labels=FALSE, col=cols[2])
    }
    if (method=="d0" & legend) {
      legend("bottomright", legend = c("xi - x0", "yi - y0"),
             col = cols, bg="#FFFFFFCC", cex=0.7, lty = 1, lwd=1)}
    if (method=="di" & legend) {
      legend("bottomright", legend = c("dx", "dy"),
             col = cols, bg="#FFFFFFCC", cex=0.7, lty = 1, lwd=1)}
    if (ref){
      #text((1:nr)[refs], dx[refs], labels=as.character(1:ref.nb), cex=0.7, col=cols[1])
      #text((1:nr)[refs], dy[refs], labels=as.character(1:ref.nb), cex=0.7, col=cols[2])
      text((1:nr)[refs], 0, labels=as.character(1:ref.nb), cex=0.7)
      }
    }
  return(list(x=dx, y=dy))}

coo.oscillo1    <- function(coo, method=c("d0", "di")[1], plot=TRUE,
                            rug=TRUE, legend=TRUE, cols=col.gallus(1), 
                            xlab="Points sampled along the outline", ylab="Deviation", ...){
  if (!is.numeric(coo)) stop("A numeric must be provided to coo.oscillo1")
  nr <- length(coo)
  if (method=="d0"){
    dx <- coo - coo[1]
  } else {
    if (method=="di"){
      dx <- coo - coo[c(nr, (1:nr - 1))]
    } else {stop("inappropriate method in oscillo")}}
  if (plot) {
    plot(dx, type="l", ylim=range(dx)*1.1,
         xlab=xlab, ylab=ylab, las=1, col=cols, ...)
    abline(h=0, col="grey80", lty=2)
    if (method=="d0" & rug) {
      dx.i <- coo - coo[c(nr, (1:nr - 1))]
      axis(3, at=(1:nr)[dx.i>0], line=-0.5, labels=FALSE, col=cols[1])
    }
    if (method=="d0" & legend) {
      legend("bottomright", legend = "xi - x0",
             col = cols, bg="#FFFFFFCC", cex=0.7, lty = 1, lwd=1)}
    if (method=="di" & legend) {
      legend("bottomright", legend = "dx",
             col = cols, bg="#FFFFFFCC", cex=0.7, lty = 1, lwd=1)}}
  return(list(x=dx))}

coo.amplify <- function(coo, amp=rep(0.5, 4), nb.h, draw=FALSE, ...){
  if (is.list(coo)) {coo <- l2m(coo)}
  if (length(amp) == 1) {amp <- rep(amp, 4)}
  if (missing(nb.h)) {nb.h <- dim(coo)[1]/2}
  coo.ef <- efourier(coo, nb.h=nb.h, smooth.it=0)
  coo.ef.amp <- ef.amplify(coo.ef, amp=amp)
  coo.amp <- efourier.i(coo.ef.amp, nb.pts=nrow(coo))
  if (draw) {
    coo.draw(coo.amp, ...) }
  return(coo.amp)}

l2m            <- function(l) {return(cbind(l$x, l$y))}
m2l            <- function(m) {return(list(x=m[,1], y=m[,2]))}
edm            <- function(m1, m2){return(sqrt((m1[, 1] - m2[, 1])^2 + (m1[, 2] - m2[, 2])^2))}
l2a            <- function(l){  return(array(unlist(l), dim=c(nrow(l[[1]]), ncol(l[[1]]), length(l))))}
# Import section #####################################################
import.txt <- function(txt.list){
  cat("Extracting", length(txt.list), ".jpg outlines...\n")
  if (length(txt.list) > 10) {
    pb <- txtProgressBar(1, length(txt.list))
    t <- TRUE } else {t <- FALSE}
  res <- list()
  for (i in seq(along = txt.list)) {
    coo <- read.table(txt.list[i], header=FALSE)
    res[[i]] <- as.matrix(coo)
    if (t) setTxtProgressBar(pb, i)
  }
  names(res) <- substr(txt.list, start=1, stop=nchar(txt.list)-4)
  return(res)}

import.jpg <- function(jpg.list) {
  cat("Extracting", length(jpg.list), ".jpg outlines...\n")
  if (length(jpg.list) > 10) {
    pb <- txtProgressBar(1, length(jpg.list))
    t <- TRUE } else {t <- FALSE}
  res <- list()
  for (i in seq(along=jpg.list)) {
    img <- import.img.prepare(jpg.list[i])
    res[[i]] <- import.img.Conte(img)
    if (t) setTxtProgressBar(pb, i)
  }
  names(res) <- substr(jpg.list, start=1, stop=nchar(jpg.list)-4) 
  return(res)}

import.multi1.jpg <- function(path){
  img <- import.img.prepare(path)
  res <- list()
  on.exit(return(res))
  i <- 1
  res[[i]] <- import.img.Conte(img, auto=FALSE)
  cat("1 outline extracted. Press 'ESC' to finish.\n")
  while (TRUE) {
    i <- i+1
    res[[i]] <- import.img.Conte(img, auto=FALSE, plot=FALSE)
    cat(paste(i, "outlines extracted.\n"))
  }  
return(res)}

import.img.prepare <- function(path){
  img <- read.jpeg(path)
  if (class(img)[2] == "array") {img <- rgb2grey(img)}
  img[img >  0.5] <- 1
  img[img <= 0.5] <- 0
  img <- imagematrix(img)
  # if(any(dim(img)>500)) {cat("\t(large image)")}
  if(any(c(img[1, ], img[nrow(img), ], img[, 1], img[, ncol(img)]) != 1)){
    # cat("\t(outline spans image border)")
    img <- rbind(rep(1, ncol(img)), img, rep(1, ncol(img)))
    img <- cbind(rep(1, nrow(img)), img, rep(1, nrow(img)))
    img <- imagematrix(img)}              
  return(img)}

import.img.Conte <- 
  function (img, x, auto=TRUE, plot=TRUE) 
  {
    if (class(img)[1] != "imagematrix") {
      stop("An 'imagematrix' object is expected")
    }
    if (missing(x)) {
      if (auto) {
        x <- round(dim(img)/2)
      } else { x <- c(1, 1)}}
    while (img[x[1], x[2]] != 0) {
      if (plot) {
        plot(img, main = "Click a point within the shape")
        rect(0, 0, ncol(img), nrow(img), border = "red")}
      click <- lapply(locator(1), round)
      x <- c(nrow(img) - click$y, click$x)
      if (any(x > dim(img))) {
        x <- round(dim(img)/2)
      }
    }
    while (abs(img[x[1], x[2]] - img[x[1], (x[2] - 1)]) < 0.1) {
      x[2] <- x[2] - 1
    }
    a <- 1
    M <- matrix(c(0, -1, -1, -1, 0, 1, 1, 1, 1, 1, 0, -1, -1, 
                  -1, 0, 1), 2, 8, byrow = TRUE)
    M <- cbind(M[, 8], M, M[, 1])
    X <- 0
    Y <- 0
    x1 <- x[1]
    x2 <- x[2]
    SS <- NA
    S <- 6
    while ((any(c(X[a], Y[a]) != c(x1, x2)) | length(X) < 3)) {
      if (abs(img[x[1] + M[1, S + 1], x[2] + M[2, S + 1]] - 
        img[x[1], x[2]]) < 0.1) {
        a <- a + 1
        X[a] <- x[1]
        Y[a] <- x[2]
        x <- x + M[, S + 1]
        SS[a] <- S + 1
        S <- (S + 7)%%8
      }
      else if (abs(img[x[1] + M[1, S + 2], x[2] + M[2, S + 
        2]] - img[x[1], x[2]]) < 0.1) {
        a <- a + 1
        X[a] <- x[1]
        Y[a] <- x[2]
        x <- x + M[, S + 2]
        SS[a] <- S + 2
        S <- (S + 7)%%8
      }
      else if (abs(img[x[1] + M[1, S + 3], x[2] + M[2, S + 
        3]] - img[x[1], x[2]]) < 0.1) {
        a <- a + 1
        X[a] <- x[1]
        Y[a] <- x[2]
        x <- x + M[, S + 3]
        SS[a] <- S + 3
        S <- (S + 7)%%8
      }
      else {
        S <- (S + 1)%%8
      }
    }
    return(cbind((Y[-1]), ((dim(img)[1] - X))[-1]))
  }

# elliptical Fourier analysis ########################################
efourier  <- function (coo, nb.h = 32, smooth.it = 0) {
  if (is.matrix(coo)) coo <- m2l(coo)
  coo <- coo.smooth(coo, smooth.it)
  coo <- coo.sample(coo, nb.h * 2)
  p <- length(coo$x)
  Dx <- coo$x - coo$x[c(p, (1:p - 1))]
  Dy <- coo$y - coo$y[c(p, (1:p - 1))]
  Dt <- sqrt(Dx^2 + Dy^2)
  t1 <- cumsum(Dt)
  t1m1 <- c(0, t1[-p])
  T <- sum(Dt)
  an <- bn <- cn <- dn <- numeric(nb.h)
  for (i in 1:nb.h) {
    Ti <- (T/(2 * pi^2 * i^2))
    r <- 2 * i * pi
    an[i] <- Ti * sum((Dx/Dt) * (cos(r * t1/T) - cos(r * t1m1/T)))
    bn[i] <- Ti * sum((Dx/Dt) * (sin(r * t1/T) - sin(r * t1m1/T)))
    cn[i] <- Ti * sum((Dy/Dt) * (cos(r * t1/T) - cos(r * t1m1/T)))
    dn[i] <- Ti * sum((Dy/Dt) * (sin(r * t1/T) - sin(r * t1m1/T)))
    }
  ao <- 2 * sum(coo$x * Dt/T)
  co <- 2 * sum(coo$y * Dt/T)
  return(list(an = an, bn = bn, cn = cn, dn = dn, ao = ao, co = co))}



efourier.i <- function(ef, nb.h, nb.pts = 300) {
  #if (any(names(ef) != c("an", "bn", "cn", "dn"))) {
  #  stop("a list containing 'an', 'bn', 'cn' and 'dn' harmonic coefficients must be provided")}
  if (is.null(ef$ao)) ef$ao <- 0
  if (is.null(ef$co)) ef$co <- 0
  an <- ef$an
  bn <- ef$bn
  cn <- ef$cn
  dn <- ef$dn
  ao <- ef$ao
  co <- ef$co
  if (missing(nb.h)) nb.h <- length(an)
  theta <- seq(0, 2 * pi, length = nb.pts + 1)[-(nb.pts + 1)]
  hx <- matrix(NA, nb.h, nb.pts)
  hy <- matrix(NA, nb.h, nb.pts)
  for (i in 1:nb.h) {
    hx[i, ] <- an[i] * cos(i * theta) + bn[i] * sin(i * theta)
    hy[i, ] <- cn[i] * cos(i * theta) + dn[i] * sin(i * theta)}
  x <- (ao/2) + apply(hx, 2, sum)
  y <- (co/2) + apply(hy, 2, sum)
  list(x = x, y = y)}

efourier.norm <- function(ef, start = FALSE) {
  A1 <- ef$an[1]
  B1 <- ef$bn[1]
  C1 <- ef$cn[1]
  D1 <- ef$dn[1]
  nb.h <- length(ef$an)
  theta      <- 0.5 * atan(2 * (A1 * B1 + C1 * D1)/(A1^2 + C1^2 - B1^2 - D1^2)) %% pi
  phaseshift <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  M2 <- matrix(c(A1, C1, B1, D1), 2, 2) %*% phaseshift
  v <- apply(M2^2, 2, sum)
  if (v[1] < v[2]) {theta <- theta + pi/2}
  theta <- (theta + pi/2)%%pi - pi/2
  Aa <- A1*cos(theta) + B1*sin(theta)
  Cc <- C1*cos(theta) + D1*sin(theta)
  scale <- sqrt(Aa^2 + Cc^2)
  psi   <- atan(Cc/Aa)%%pi
  if (Aa<0){psi<-psi+pi}
  size  <- 1/scale
  rotation <- matrix(c(cos(psi), -sin(psi), sin(psi), cos(psi)), 2, 2)
  A <- B <- C <- D <- numeric(nb.h)
  if (start) {theta <- 0}
  lnef <- NULL
  for (i in 1:nb.h) {
    mat <- size * rotation %*%
      matrix(c(ef$an[i], ef$cn[i], ef$bn[i], ef$dn[i]), 2, 2) %*%
      matrix(c(cos(i*theta), sin(i*theta), -sin(i*theta), cos(i*theta)), 2, 2)
      A[i] <- mat[1, 1]
      B[i] <- mat[1, 2]
      C[i] <- mat[2, 1]
      D[i] <- mat[2, 2]
      lnef <- c(lnef, c(A[i], B[i], C[i], D[i]))}
  list(A = A, B = B, C = C, D = D, size = scale, theta = theta, 
      psi = psi, ao = ef$ao, co = ef$co, lnef = lnef)}

efourier.shape <- function(an, bn, cn, dn, nb.h, nb.pts=80, alpha=2, plot=TRUE){
  if (missing(nb.h) &  missing(an)) nb.h <- 1
  if (missing(nb.h) & !missing(an)) nb.h <- length(an)
  if (missing(an)) an <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  if (missing(bn)) bn <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  if (missing(cn)) cn <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  if (missing(dn)) dn <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  ef  <- list(an=an, bn=bn, cn=cn, dn=dn, ao=0, co=0)
  shp <- efourier.i(ef, nb.h=nb.h, nb.pts=nb.pts)      
  if (plot) coo.plot(shp)
  return(shp)}

ef.amplify <- function(ef, amp=rep(0.5, 4)){
  ef$an <- ef$an*amp[1]
  ef$bn <- ef$bn*amp[2]
  ef$cn <- ef$cn*amp[3]
  ef$dn <- ef$dn*amp[4]
  return(ef)
}

# Radius variation Fourier ###########################################
rfourier<-function(coo, nb.h=32, smooth.it=0){
  if (is.list(coo)) coo <- l2m(coo)
  coo <- coo.smooth(coo, smooth.it)
  coo <- coo.sample.rr(coo, nb.h * 2)$coord
  p <- nrow(coo)
  bn <- an <- numeric(nb.h)
  Z  <- complex(real=coo[, 1], imaginary=coo[, 2])
  r     <- Mod(Z)
  angle <-Arg(Z)
  ao <- 2* sum(r)/p
  for (i in 1:nb.h){
    an[i]<-(2/p)*sum(r * cos(i*angle))
    bn[i]<-(2/p)*sum(r * sin(i*angle))}
  list(an=an, bn=bn, ao=ao, r= r)}

rfourier.i<-function(rf, nb.pts=300, nb.h=length(rf$an)) {
  an <- rf$an
  bn <- rf$bn
  ao <- rf$ao
  theta <- seq(0, 2*pi, length=nb.pts+1)[-(nb.pts+1)]
  harm  <- matrix (NA, nb.h, nb.pts)
  for (i in 1:nb.h){
    harm[i, ] <-an[i]*cos(i*theta)+ bn[i]*sin(i*theta)}
  r<-(ao/2) + apply(harm, 2, sum)
  Z<-complex(modulus=r, argument=theta)
  list(x = Re(Z), y = Im(Z), angle = theta, r = r)}

rfourier.shape <- function(an, bn, nb.h, nb.pts=80, alpha=2, plot=TRUE){
  if (missing(nb.h) &  missing(an)) nb.h <- 1
  if (missing(nb.h) & !missing(an)) nb.h <- length(an)
  if (missing(an)) an <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  if (missing(bn)) bn <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  rf  <- list(an=an, bn=bn, ao=0, co=0)
  shp <- rfourier.i(rf, nb.h=nb.h, nb.pts=nb.pts)      
  if (plot) coo.plot(shp)
  return(shp)}

# Tangent angle Fourier ##############################################
tfourier<-function(coo, nb.h=32, smooth.it=0){
  if (is.list(coo)) coo <- l2m(coo)
  coo <- coo.smooth(coo, smooth.it)
  coo <- coo.sample.rr(coo, nb.h * 2)$coord
  p   <- nrow(coo)
  bn  <- an <- numeric(nb.h)
  tangvect  <- coo-rbind(coo[p,], coo[-p,])
  perim     <- sum(sqrt(apply((tangvect)^2, 1, sum)))
  v0        <- (coo[1,]-coo[p,])
  tet1 <- Arg(complex(real=tangvect[, 1], imaginary = tangvect [, 2]))
  tet0 <- tet1[1]
  t1   <- (seq(0, 2*pi, length=(p+1)))[1:p]
  phi  <- (tet1-tet0-t1)%%(2*pi)
  ao<- 2* sum(phi)/p
  for (i in 1:nb.h){
    an[i]<- (2/p) * sum(phi * cos (i*t1))
    bn[i]<- (2/p) * sum(phi * sin (i*t1))}
  list(ao=ao, an=an, bn=bn, phi=phi, t=t1, perimeter=perim, thetao=tet0)}

tfourier.i<-function(tf, nb.pts, nb.h=length(tf$an), thetao=0){
  ao <- tf$ao
  an <- tf$an
  bn <- tf$bn
  theta <- seq(0, 2*pi, length=nb.pts+1)[-(nb.pts+1)]
  harm  <- matrix (NA, nb.h, nb.pts)
  for (i in 1:nb.h){
    harm[i,] <- an[i]*cos(i*theta) + bn[i]*sin(i*theta)}
  phi  <- (ao/2) + apply(harm, 2, sum)
  vect <- matrix(NA, 2, nb.pts)
  Z    <- complex(modulus=(2*pi)/nb.pts, argument=phi+theta+thetao)
  Z1   <- cumsum(Z)
  list(x=Re(Z1), y=Im(Z1), angle=theta, phi=phi)}

tfourier.shape <- function(an, bn, ao=0, nb.h, nb.pts=80, alpha=2, plot=TRUE){
  if (missing(nb.h) &  missing(an)) nb.h <- 1
  if (missing(nb.h) & !missing(an)) nb.h <- length(an)
  if (missing(an)) an <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  if (missing(bn)) bn <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  tf  <- list(an=an, bn=bn, ao=ao)
  shp <- tfourier.i(tf, nb.h=nb.h, nb.pts=nb.pts)      
  if (plot) coo.plot(shp)
  return(shp)}

# Coefficient matrices utils #########################################
coeff.sel <- function(retain=8, drop=0, nb.h=32, cph=4){
  cs <- numeric()
  for (i in 1:cph) {
    cs <- c(cs, (1+drop):retain + nb.h*(i-1))}
  return(cs)}

coeff.split <- function(cs, nb.h=8, cph=4){
  if (missing(nb.h)) {nb.h <- length(cs)/cph }
  cp <- list()
  for (i in 1:cph) {
    cp[[i]] <- cs[1:nb.h + (i-1)*nb.h]
    }
  names(cp) <- paste(letters[1:cph], "n", sep="")
  return(cp)}

# Miscellaneous ######################################################
# Color palettes
col.summer <- colorRampPalette(c("#4876FF", "#FFFF00", "#FF3030"))
col.gallus <- colorRampPalette(c("#000080", "#FFFFFF", "#EE0000"))
col.sari   <- colorRampPalette(c("#551A8B", "#FF7F00"))
col.india  <- colorRampPalette(c("#FF9933", "#138808"))
col.nb     <- colorRampPalette(c("#FFFFFF", "#000000"))
col.wcol   <- function(col.hex) colorRampPalette(c("#FFFFFF", col.hex))
col.bcol   <- function(col.hex) colorRampPalette(c("#000000", col.hex))

pca2shp <- function (pos, rot, mean.shp, scale=1, amp=1, trans=TRUE, nb.pts=64) {
  # we check a bit
  if (!is.matrix(pos))        pos <- as.matrix(pos)
  if (ncol(pos) != ncol(rot)) stop("rot an pos must have the same ncol")
  if(length(mean.shp) != nrow(rot)) stop("mean.shp length must equals the col number of rot")
  # stupid function
  mprod <- function(m, s){
    res <- m
    for (i in 1:ncol(m)) {
      res[, i] <- m[, i]*s[i]}
    return(res)}
  res <- list()
  nb.h <- length(mean.shp)/4
  n  <- nrow(pos)
  for (i in 1:n) {
    ax.contrib <- mprod(rot, pos[i, ])*amp
    coe        <- mean.shp + apply(ax.contrib, 1, sum)
    coo        <- l2m(efourier.shape(coe[1:nb.h + 0*nb.h], coe[1:nb.h + 1*nb.h],
                                     coe[1:nb.h + 2*nb.h], coe[1:nb.h + 3*nb.h],
                                     nb.h = nb.h, nb.pts=nb.pts, plot=FALSE))
    coo <- coo.scale(coo, scale=scale)
    if (trans) { res[[i]] <- coo.trans(coo, pos[i, 1], pos[i, 2]) } else { res[[i]] <- coo }
  }
  return(res)}

dev.plot       <- function(mat, matv, cols, x=1:ncol(mat), bw=0.1) {
  r <- nrow(mat)
  if (missing(cols) | length(cols)!=r) {cols <- rep("black", r)}
  if (!missing(matv)){
    if (any(dim(x)!=dim(x))) {stop("mat and matv must be of the same dimension")}
    for (i in 1:nrow(mat)){
      segments(x,    mat[i, ] - matv[i, ], x, mat[i, ]    + matv[i, ], col=cols[i], lwd=0.5)
      segments(x-bw, mat[i, ] - matv[i, ], x+bw, mat[i, ] - matv[i, ], col=cols[i], lwd=0.5)
      segments(x-bw, mat[i, ] + matv[i, ], x+bw, mat[i, ] + matv[i, ], col=cols[i], lwd=0.5)
    }
  }
  for (i in 1:nrow(mat)) {
    lines(x, mat[i, ], col=cols[i], type="o", pch=20, cex=0.6)}}


dudi.plot <- function(dudi, fac = NULL, xax = 1, yax = 2, grid = TRUE,
                      points     = TRUE,  pch.points=1,  col.points="black", cex.points=0.8,
                      labels     = FALSE, label=rownames(dudi$li), boxes=TRUE, clabel=1,
                      neighbors  = FALSE, col.nei="grey90",  lwd.nei=0.5,
                      star       = TRUE,  col.star="grey60", cstar=1,
                      ellipses   = TRUE,  col.ellipse="grey30", cellipse=1, axesell=TRUE,
                      chull      = FALSE, col.chull="grey30", optchull = c(0.5, 1),
                      arrows     = FALSE, edge.arrow=FALSE, box.arrow=TRUE, maxnb.arrow=10, dratio.arrow=0.2, 
                      shapes     = TRUE,  pos.shp=c("li", "usr", "circle", "full", "range")[4],
                      nr.shp=6, nc.shp=5, amp.shp=1, scale.shp=0.8,
                      circle.nb.shp=12, circle.r.shp,
                      col.shp="#70809011", border.shp="#708090",      
                      rug        = TRUE, rug.ticksize=0.01, rug.col="#708090",
                      eigen      = FALSE, eigen.ratio=0.2,
                      palette    = col.sari,
                      title      = substitute(dudi),
                      center.orig= FALSE,
                      zoom.plot  = 1){
  
  # we prepare and check a bit
  if (!missing(fac)) {
  if (ncol(dudi$fac)==0) { fac <- factor(rep("", nrow(dudi$li))) } else {fac <- dudi$fac[, fac]}
  if ((nlevels(fac) > 1)) {
    if (missing(col.star))    col.star    <- paste(palette(nlevels(fac)), "33", sep="")
    if (missing(col.ellipse)) col.ellipse <- palette(nlevels(fac))
    if (missing(col.chull))   col.chull   <- palette(nlevels(fac))} }
  
  # we initialize the factorial map
  if (center.orig) {
    li.2      <- apply(dudi$li[,c(xax, yax)], 2, function(x) x^2)
    li.len    <- apply(li.2, 1, function(x) sqrt(sum(x)))
    w <- max(li.len)*(1/zoom.plot)
    s.label(dudi$li, xax=xax, yax=yax, xlim=c(-w, w), clabel=0, cpoint=0, sub=title, grid=grid)
  } else {     
    s.label(dudi$li, xax=xax, yax=yax, clabel=0, cpoint=0, sub=title, grid=grid)}
  
  # size of the grid
  xaxp <- par("xaxp")
  ax <- (xaxp[2] - xaxp[1])/xaxp[3]
  yaxp <- par("yaxp")
  ay <- (yaxp[2] - yaxp[1])/yaxp[3]
  d <- min(ax, ay)
  
  # we redefine the right margins
  op <- par("mar")
  par(mar=rep(0.1, 4))
  
  # rug
  if (rug) {
    rug(dudi$li[, xax], side=1, ticksize=rug.ticksize, col=rug.col, lwd=0.4)
    rug(dudi$li[, yax], side=2, ticksize=rug.ticksize, col=rug.col, lwd=0.4)
    box()}
  
  # neighboring graph
  if (neighbors) {
    fun <- function(x, coo, col, lwd) {
      segments(coo$x[x[1]], coo$y[x[1]],
               coo$x[x[2]], coo$y[x[2]], 
               col = col, lwd = lwd)}
    neig <- nb2neig(tri2nb(dudi$li[, c(xax, yax)]))
    coo  <- list(x=dudi$li[, xax], y=dudi$li[, yax])
    apply(unclass(neig), 1, fun,
          coo = coo, col=col.nei, lwd=lwd.nei)}
  
  # star*
  if (star & !is.null(fac)) {
    s.class(dudi$li, xax=xax, yax=yax, fac=fac, clabel=0, cpoint=0, add.plot=TRUE,
            cstar=cstar, col=col.star, cellipse=0)}
  
  # ellipses*
  if (ellipses & !is.null(fac)) {
    s.class(dudi$li,  xax=xax, yax=yax, fac=fac,
            clabel=0, cpoint=0, add.plot=TRUE,
            cstar=0, col=col.ellipse, cellipse=cellipse, axesell=axesell)}
  
  # chull*
  if (chull & !is.null(fac)) {
    s.chull(dudi$li,xax=xax, yax=yax, fac=fac, col=col.chull, 
            optchull=optchull, add.plot=TRUE)}
  
  # arrow
  if (arrows) {
    arr.2      <- apply(dudi$co[,c(xax, yax)], 2, function(x) x^2)
    arr.len    <- apply(arr.2, 1, function(x) sqrt(sum(x)))
    arr.sorted <- order(arr.len, decreasing=TRUE)[1:maxnb.arrow]
    arr.disp   <- arr.sorted[ arr.len[arr.sorted] > d*dratio.arrow ]
    s.arrow(dudi$co[arr.disp,],
            label = rownames(dudi$co[arr.disp, c(xax, yax)]), edge=edge.arrow, add.plot=TRUE, boxes=box.arrow)}
  # shapes
  if (shapes & !is.null(dudi$mean.shp)) {
    if (pos.shp=="li") {
      pos <- dudi$li[, c(xax, yax)]
    } else if (pos.shp=="usr") {
    } else if (pos.shp=="circle") {
      if (missing(circle.r.shp)) {
        li.2      <- apply(dudi$li[,c(xax, yax)], 2, function(x) x^2)
        li.len    <- apply(li.2, 1, function(x) sqrt(sum(x)))
        circle.r.shp <- mean(li.len)}
      t <- seq(0, 2*pi, len=circle.nb.shp+1)[-(circle.nb.shp+1)]
      pos <- cbind(circle.r.shp*cos(t), circle.r.shp*sin(t))
    } else if (pos.shp=="range") {
      pos <- expand.grid(seq(min(dudi$li[, xax]), max(dudi$li[, xax]), len=nr.shp),
                         seq(min(dudi$li[, yax]), max(dudi$li[, yax]), len=nc.shp))
    } else if  (pos.shp=="full") {
      w <- par("usr")
      pos <- expand.grid(seq(w[1]+d, w[2]-d, len=nr.shp),
                         seq(w[3]+d, w[4]-d, len=nc.shp))
    } else {stop("shp.pos must be passed with values li, usr, circle, full or range")}    
    x <- pca2shp(pos, rot=dudi$c1[, c(xax, yax)],
                 mean.shp=dudi$mean.shp, scale=d*scale.shp, amp=amp.shp)
    b.quiet <- lapply(x, polygon, col=col.shp, border=border.shp, pch=20)
  }
  
  # labels and points
  if (points) {
  repeach <- function(x, each){
    if (length(x) != length(each)) return(rep(x[1], sum(each)))
    res <- vector(mode = class(x[1]))
    for (i in seq(along=x)) {
      res <- append(res, rep(x[i], each[i]))}
    return(res)}
  if (!is.null(fac)) {
    nb <- table(fac)
    pch.points <- repeach(pch.points, nb)
    if (missing(col.points)) {col.points <- repeach(palette(nlevels(fac)), nb)}
    cex.points <- repeach(cex.points, nb)} 
  points(dudi$li[, c(xax, yax)], pch=pch.points, col=col.points, cex=cex.points)
  }
  
  #  s.label(dudi$li, xax=xax, yax=yax, clabel=0, cpoint=cpoint, pch=pch, add.plot=TRUE)}
  # labels
  if (labels) {
    s.label(dudi$li, xax=xax, yax=yax, clabel=clabel, cpoint=0, boxes=boxes, add.plot=TRUE)}
  
  # ellipses* (only for the labels) #probably not the most orthodox option
  if (ellipses & !is.null(fac)) {
    s.class(dudi$li,  xax=xax, yax=yax, fac=fac,
            clabel=clabel, cpoint=0, add.plot=TRUE,
            cstar=0, col=NA, cellipse=0, axesell=FALSE)}
  
  # eigen
  if (eigen) {
    par("mar"=op)
    add.scatter.eig(dudi$eig, nf=dudi$nf, xax=xax, yax=yax, eigen.ratio, posi="bottomright")}
  
  # we restore the margins
  par("mar"=op)
}

ellpar <- function(coo){
  if (is.list(coo)) coo <- cbind(coo$x, coo$y) ### ML
  coo <- coo %*% svd(var(coo))$u
  a <- max(coo[, 1])
  b <- max(coo[, 2])
  e <- sqrt((a^2 - b^2)/a^2)
  list(a=a, b=b, e=e)}

