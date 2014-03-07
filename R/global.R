#--------------------------------------------------------------------#
# Global functions                                                   #
#--------------------------------------------------------------------#

# coo utils ##########################################################

coo.center     <- function(coo){
  if (is.matrix(coo)) {  
    return(apply(coo, 2, function(x) x-mean(x)))
  } else if (is.list(coo)){
    return(lapply(coo, function(x) x-mean(x)))
  } else stop("A list or a coordinate matrix must be provided to coo.center")}

coo.scale <- 
  function (coo, scale) 
  {
    if (missing(scale)) {
      scale <- ifelse(is.list(coo), coo.centsize(l2m(coo)), coo.centsize(coo))}
    if (is.list(coo)) {
      return(lapply(coo, function(x) x/scale))}
    if (is.matrix(coo)) {
      return(coo/scale)}
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

coo.sample.int <- function(coo, n){
  if (n < 3) { stop("n must be >= 3") }
  if (is.list(coo)) { coo <- l2m(coo) }
  if (!is.closed(coo)) { coo <- coo.close(coo) }
  orig <- cumsum(coo.perim.pts(coo))
  targ <- seq(0, coo.perim(coo), length=n+1)[-(n+1)]
  res <- matrix(c(coo[1, ], rep(NA, n*2 - 2)), byrow=TRUE,
                nrow=n, ncol=2, dimnames=c(list(paste0("pt", 1:n), c("x", "y"))))
  for (i in 2:n) {
    k <- max(which(orig <= targ[i]))
    r <- (targ[i] - orig[k]) / (orig[k+1]- orig[k])
    res[i, ] <- edi(coo[k, ], coo[k+1, ], r)}
  return(res)}

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
    if (is.closed(coo)) { return(coo) } else { return(rbind(coo, coo[1, ])) }
  } else if (is.list(coo)){
    if (is.closed(coo)) { return(coo) } else {
    return(list(x=c(coo$x, coo$x[1]), y=c(coo$y, coo$y[1])))}
  } else stop("A list or a coordinate matrix must be provided to coo.center")}

coo.unclose      <- function(coo){
  if(!is.closed(coo)) return(coo)
  if (is.matrix(coo)) {  
    return(coo[-nrow(coo), ])
  } else if (is.list(coo)){
    return(list(x=coo$x[-length(coo$x)], y=coo$y[-length(coo$y)]))
  } else stop("A list or a coordinate matrix must be provided to coo.center")}

is.closed <- function(coo){
  if (is.list(coo)) {
    n <- length(coo$x)
    return((coo$x[1] == coo$x[n]) & (coo$y[1] == coo$y[n]))}
  if (is.matrix(coo)) {
    n <- nrow(coo)
    return((coo[1, 1] == coo[n, 1]) & (coo[1, 2] == coo[n, 2]))}}   

coo.ldk <- function(coo, nb.ldk) {
  if (is.list(coo)) coo <- l2m(coo)
  coo.plot(coo)
  ldk <- numeric(nb.ldk)
  cat("[")
  for (i in 1:nb.ldk){
    p <- l2m(locator(1))
    l <- apply(coo, 1, function(y) sqrt(sum((p-y)^2)))
    ldk[i] <- which.min(l)
    points(coo[ldk[i], 1], coo[ldk[i], 2], pch=20, col="red")
    cat("*")
    }
  cat("]\n")
  return(ldk)}

coo.centpos <- function(coo){
  return(apply(coo, 2, mean))}

coo.centsize <- function(coo){
  cent <- coo.centpos(coo)
  return(mean(apply(coo, 1, function(x) sqrt(sum((x-cent)^2)))))}

coo.baseline <- function(coo, ldk1=1, ldk2=2, t1=c(-0.5, 0), t2=c(0.5, 0)){
  if (is.list(coo)) {coo <- l2m(coo)}
  t1x <- t1[1]
  t1y <- t1[2]
  t2x <- t2[1]
  t2y <- t2[2]
  r1x <- coo[ldk1, 1]
  r1y <- coo[ldk1, 2]
  r2x <- coo[ldk2, 1]
  r2y <- coo[ldk1, 2]
  # translation based on the first landmark
  ref <- coo.trans(coo, t1x - coo[ldk1, 1] , t1y - coo[ldk1, 2])
  # we calculate dx and dy for the two vectors
  rx <- ref[ldk2, 1] - t1x
  ry <- ref[ldk2, 2] - t1y 
  tx <- t2x - t1x
  ty <- t2y - t1y
  vi <- vecs.param(rx, ry, tx, ty)
  # we rotate accordingly with a center defined as the first landmark (trans, rot, untrans)
  ref <- coo.trans(ref, -t1x, -t1y)
  ref <- ref / vi$r.norms
  ref <- coo.rotate(ref, -vi$d.angle)
  ref <- coo.trans(ref, t1x, t1y)
  return(ref)}

coo.rotate.center <- function(coo, theta, center=c(0, 0)){
  coo <- coo.trans(coo, -center[1], -center[2])
  coo <- coo.rotate(coo, theta)
  coo <- coo.trans(coo, center[1], center[2])
  return(coo)}

coo.perim <- function(coo){
  if (is.list(coo)) { coo <- l2m(coo) }
  n <- nrow(coo)
  perim <- sum(sqrt(apply((coo-coo.slide(coo, n))^2, 1, sum)))
  return(perim)}

coo.perim.pts <-  function (coo){
    if (is.list(coo)) {
      coo <- l2m(coo)
    }
    n <- nrow(coo)
    perim <- sqrt(apply((coo - coo.slide(coo, n))^2, 1, sum))
    return(perim)}

coo.force2close <- function(coo){
  if (is.list(coo)) { coo <- l2m(coo)}
  if (is.closed(coo)) {return(coo)}
  n  <- nrow(coo)
  d  <- coo[1, ] - coo[n, ]
  dm <- cbind(seq(0, d[1], length=n), seq(0, d[2], length=n))
  return(coo+dm)}


# coo plotting functions #############################################

coo.plot       <- function(coo=NA, col="#F5F5F5",
                           border="#1A1A1A", lwd=1, lty = 1, xlim=c(-1, 1), ylim=c(-1, 1),
                           points=FALSE, first.point=TRUE, centroid=TRUE,
                           points.col=border, pch=1, cex=0.8, main, ...){
  #   rewrite below without the entire par list
  #   op <- par(no.readonly = TRUE)
  #   on.exit(par(op))
  #   if (missing(main)) {
  #     par(mar=c(2.2, 2.2, 0.6, 0.6))
  #   } else {
  #     par(mar=c(2.2, 2.2, 2.2, 1.2))}
  if (is.list(coo)) coo <- l2m(coo)
  if (missing(coo)) {
    plot(NA, asp=1, xlim=xlim, ylim=ylim, ann=FALSE, frame=FALSE,
         cex.axis=0.7, cex=0.7, ...)
  } else {
	if (missing(xlim) & missing(ylim)) {
    plot(coo, type="n", asp=1, ann=FALSE, frame=FALSE,
         cex.axis=0.7, cex=0.7, las=1, ...)
	} else {
	plot(coo, type="n", asp=1, ann=FALSE, frame=FALSE,
       xlim=xlim, ylim=ylim,
	     cex.axis=0.7, cex=0.7, las=1, ...)
	}
    polygon(coo, col=col, border=border, lwd=lwd, lty=lty)
    if (first.point) {points(coo[1, 1], coo[1, 2], col = border, pch=20, ...)}
}
  if ((!missing(coo) & missing(points))) {
    if (nrow(coo)<=100) points(coo, pch=pch, cex=cex, col=points.col)}
  if ((!missing(coo) & points)) {
    points(coo, pch=pch, cex=cex, col=points.col)}
  if (!missing(coo) & centroid) {
    cent <- coo.centpos(coo)
    points(cent[1], cent[2], pch=3, cex=0.5)
    }
  if (!missing(main)) title(main=main)}  

coo.draw <-
  function (coo = NA, col = "#70809033", border = "#708090EE", lwd=1, lty=1,
            points = FALSE, first.point=TRUE, centroid = FALSE,
            points.col = border, pch = 20, cex = 0.25,  ...)
  {
    if (is.list(coo)) {coo <- l2m(coo)}
    polygon(coo, col = col, border = border, lwd=lwd, lty=lty)
    if (first.point) {points(coo[1, 1], coo[1, 2], col = border, pch=20, ...)}
    if (points) {
      points(coo, pch = pch, cex = cex, col = points.col)
    }
    if (centroid) {
      cent <- coo.centpos(coo)
      points(cent[1], cent[2], pch=3, cex=0.5)}
  }

coo.template   <- function(coo, size=1) {
  # only for matrices
  coo      <- coo * min(size/apply(coo, 2, function(x) diff(range(x))))
  expected <- apply(coo, 2, function(x) diff(range(x)))/2
  observed <- apply(coo, 2, range)[2, ]
  shift    <-  expected - observed
  return(coo.trans(coo, shift[1], shift[2]))}

coo.list.panel <- function(coo.list, dim, byrow=TRUE,
                           fromtop=TRUE, mar=rep(0, 4),
                           cols, borders, density = NULL, angle = 45){
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
  if (fromtop & dim[1]>1) { pos <- pos[dim[1]:1,] }
  # we prepare the panel
  op <- par("mar")
  par(mar=mar)
  plot(NA, asp=1,
       xlim=c(0, dim[2]),
       ylim=c(0, dim[1]),
       xaxs="i", yaxs="i", frame=FALSE, ann=FALSE, axes=FALSE)
  # we template and plot shapes
  coo.tp  <- lapply(coo.list, coo.template, size=0.9)
  if (missing(cols))    { cols      <- rep("grey80", n) }
  if (missing(borders)) { borders   <- rep("grey20", n) }
  if (missing(density)) { density   <- rep(NULL, n) }
  if (missing(angle))   { angle     <- rep(45, n) }
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
             col = cols, bg="#FFFFFFCC", cex=0.7, lty = 1, lwd=1, inset=0.005)}
    if (method=="di" & legend) {
      legend("bottomright", legend = c("dx", "dy"),
             col = cols, bg="#FFFFFFCC", cex=0.7, lty = 1, lwd=1, inset=0.005)}
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
             col = cols, bg="#FFFFFFCC", cex=0.7, lty = 1, lwd=1, inset=0.005)}
    if (method=="di" & legend) {
      legend("bottomright", legend = "dx",
             col = cols, bg="#FFFFFFCC", cex=0.7, lty = 1, lwd=1, inset=0.005)}}
  return(list(x=dx))}

coo.ef.amplify <- function(coo, amp=rep(0.5, 4), nb.h=5, draw=FALSE, ...){
  if (is.list(coo)) {coo <- l2m(coo)}
  if (length(amp) == 1) {amp <- rep(amp, 4)}
  if (missing(nb.h)) {nb.h <- floor(dim(coo)[1]/2)-1}
  coo.ef <- efourier(coo, nb.h=nb.h, smooth.it=0)
  coo.ef.amp <- ef.amplify(coo.ef, amp=amp)
  coo.amp <- efourier.i(coo.ef.amp, nb.pts=nrow(coo))
  if (draw) {
    coo.draw(coo.amp, ...) }
  return(coo.amp)}

#  A utilities & plotting ############################################
A.mshape<-function(A){apply(A, c(1,2), mean)}

A.plot <- function(A, col, palette=col.summer, xlim, ylim,
                  border.col=NA, border.lty=2, pch=3, cex=1){
  if (!is.array(A)) {stop("An array must be provided.")}
  if (!is.array(A)) {stop("An 3 dimension array must be provided.")}
  if (missing(col)) {
    cols <- palette(dim(A)[1])
  } else {
    cols <- rep(col, dim(A)[1])
  }
  if (missing(xlim) & missing(ylim)) {
    wdw  <- apply(apply(A, 1:2, range), 3, range)
    xlim <- wdw[, 1]
    ylim <- wdw[, 2]}
  plot(NA, xlim=xlim, ylim=ylim, asp=1, ann=FALSE, frame=FALSE, cex.axis=0.7)
  A.segments(A, border.col=border.col, border.lty=border.lty)
  be.quiet <- apply(A, 3, points, pch=pch, cex=cex, col=cols)}

A.points <- function(A, col, palette=col.summer, pch=3, cex=1){
  if (!is.array(A)) {stop("An array must be provided.")}
  if (!is.array(A)) {stop("An 3 dimension array must be provided.")}
  if (missing(col)) {
    cols <- palette(dim(A)[1])
  } else {
    cols <- rep(col, dim(A)[1])
  }
  be.quiet <- apply(A, 3, points, pch=pch, cex=cex, col=cols)}


A.segments <- function(A, col=NA, border.col="grey60", border.lty=1){
  if (!is.array(A)) {stop("An array must be provided.")}
  if (!is.array(A)) {stop("An 3 dimension array must be provided.")}
  be.quiet <- apply(A, 3, coo.draw, col=col, border=border.col, lty=border.lty)}

# Other utilities and class converters ###############################

l2m            <- function(l) {return(cbind(l$x, l$y))}

m2l            <- function(m) {return(list(x=m[,1], y=m[,2]))}

ed             <- function(pt1, pt2){return(sqrt((pt1[1]-pt2[1])^1+(pt1[2]-pt2[2])^2))}

edi <- function(pt1, pt2, r=0.5){
  return(r*(pt2-pt1) + pt1) }

edm            <- function(m1, m2){return(sqrt((m1[, 1] - m2[, 1])^2 + (m1[, 2] - m2[, 2])^2))}

edm.nearest <- function(m1, m2, full=FALSE){
  if (!is.matrix(m1) | !is.matrix(m2)) stop("Matrices must be provided")
  if (ncol(m1)!=2    | ncol(m2)!=2)    stop("2-cols matrices must be provided")
  nr <- nrow(m1)
  pos <- d  <- numeric(nr)
  for (i in 1:nr){
    m1.i   <- m1[i, ]
    di     <- apply(m2, 1, function(x) sqrt(sum((x - m1.i)^2)))
    d[i]   <- min(di)
    pos[i] <- which.min(di)}
  if (full) return(list(d=d, pos=pos)) else return(d) }



l2a            <- function(l){return(array(unlist(l), dim=c(nrow(l[[1]]), ncol(l[[1]]), length(l))))}
a2l <- function(a){
  if (!is.array(a)) stop("An array of dimension 3 must be provided")
  k <- dim(a)[3]
  l <- list()
  for (i in 1:k) {l[[i]] <- a[,,i]}
  return(l)}

dev.segments <-function(coo, cols, lwd=1){
  nr <- nrow(coo)
  coo <- rbind(coo, coo[1, ])
  for (i in 1:nr) {
    segments(coo[i, 1], coo[i, 2], coo[i+1, 1], coo[i+1, 2],
             col=cols[i], lwd=lwd)}}

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

# Import section #####################################################
import.txt <- function(txt.list, ...){
  cat("Extracting", length(txt.list), ".jpg outlines...\n")
  if (length(txt.list) > 10) {
    pb <- txtProgressBar(1, length(txt.list))
    t <- TRUE } else {t <- FALSE}
  res <- list()
  for (i in seq(along = txt.list)) {
    coo <- read.table(txt.list[i], ...)
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
  img <- readJPEG(path)
  #if (class(img)[2] == "array") {img <- rgb2grey(img)} #to2fixed...ReadImages
  img[img >  0.5] <- 1
  img[img <= 0.5] <- 0
  #img <- imagematrix(img) #to2fixed...ReadImages
  # if(any(dim(img)>500)) {cat("\t(large image)")}
  if(any(c(img[1, ], img[nrow(img), ], img[, 1], img[, ncol(img)]) != 1)){
    # cat("\t(outline spans image border)")
    img <- rbind(rep(1, ncol(img)), img, rep(1, ncol(img)))
    img <- cbind(rep(1, nrow(img)), img, rep(1, nrow(img)))
    #img <- imagematrix(img) #to2fixed...ReadImages
  }              
  return(img)}

import.img.Conte <- 
  function (img, x, auto=TRUE, plot=TRUE) 
  {
    #if (class(img)[1] != "imagematrix") {
    #  stop("An 'imagematrix' object is expected")}2befixe...ReadImages
    img <- t(img[nrow(img):1,])
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

# Babel ##############################################################

pix2chc <- function(coo) {
  if (is.list(coo)) {
    coo <- l2m(coo)}
  if (is.matrix(coo) & ncol(coo)!=2) {
    stop("A 2 col matrix must be provided")}
  coo.d <- apply(coo, 2, diff)
  if (!all(coo.d %in% -1:1)) {
    stop("Matrix must contain only entire pixels indices")}
  if (any(apply(coo.d, 1, function(x) all(x==rep(0, 2))))) {
    stop("At least two succesive coordinates don't code for a displacement")}
  m   <- as.matrix(expand.grid(-1:1, -1:1))[-5,]
  g   <- c(5, 6, 7, 4, 0, 3, 2, 1)
  chc <- g[apply(coo.d, 1, function(x) which(x[1]==m[, 1] & x[2]==m[, 2]))] #dirty
  return(chc)}

chc2pix <- function(chc){
  if (!all(chc %in% 0:7)) {
    stop("chc string must only contain integers between 0 and 7")}
  m <- matrix(c(1, 0, 1, 1, 0, 1, -1, 1,
                -1, 0, -1, -1, 0, -1, 1, -1), ncol=2, byrow=TRUE)
  pix <- apply(m[chc+1,], 2, cumsum)
  return(pix)}

Coo2chc <- function(Coo, file="chc.chc"){ 
  res <- list()
  pb <- txtProgressBar(1, Coo@coo.nb)
  for (i in 1:Coo@coo.nb){
    res[[i]] <- c(Coo@names[i], rep(1, 3), Polygon(list(Coo@coo[[i]]))@area,
                  pix2chc(Coo@coo[[i]]), -1, "\n")
    setTxtProgressBar(pb, i)}
  cat(unlist(res), file=file)
  cat(".chc file succesfully written here:", file)}

chc2Coo <- function(chc.path){
  chc <- readLines(chc.path)
  coo.list <- list()
  coo.names <- character()
  for (i in seq(along=chc)) {
    chc.i <- unlist(strsplit(chc[i], " "))
    rm    <- match(chc.i, c("", " ", "-1"), , nomatch=0)
    if (any(rm)) {chc.i <- chc.i[rm==0]}
    coo.names[i] <- chc.i[1]
    st    <- as.numeric(chc.i[2:3])
    pix.i <- chc2pix(as.numeric(chc.i[-(1:5)]))
    coo.list[[i]] <- coo.trans(pix.i, st[1], st[2])
  }
  names(coo.list) <- coo.names
  return(Coo(coo.list))}

nef2Coe <- function(nef.path) {
  # change nef to coe one day
  nef     <- readLines(nef.path)
  HARMO.l <- grep(pattern="HARMO", nef)
  nb.h    <- as.numeric(substring(nef[HARMO.l], 8))
  nef     <- nef[-(1:HARMO.l)]
  nb.coo  <- length(nef)/(nb.h+1)
  coo.i   <- 1:nb.coo
  coo.beg <- (coo.i-1)*(nb.h + 1)+1
  coo.end <- coo.beg + nb.h
  res     <- matrix(NA, nrow=nb.coo, ncol=nb.h*4, dimnames=
    list(nef[coo.beg],
         paste(rep(LETTERS[1:4], each=nb.h), 1:nb.h, sep="")))
  for (i in seq(along=coo.i)) {
    nef.i    <- nef[(coo.beg[i]+1) : coo.end[i]]
    x        <- as.numeric(unlist(strsplit(nef.i, " ")))
    res[i, ] <- x[!is.na(x)]}
  return(Coe(res))}

Coe2nef <- function(Coe, file="nef.nef"){
  nb.h      <- Coe@nb.h
  coo.names <- Coe@names
  coo.i   <- 1:length(coo.names)
  coo.beg <- (coo.i-1)*(nb.h + 1)+1
  coo.end <- coo.beg + nb.h
  #nef <- c("#CONST ", constant.coeff, " \n", "#HARMO ", nb.h, " \n")
  nef <- character()
  for (i in seq(along=coo.names)) {
    nef <- append(nef, c(coo.names[i], "\n"))
    for (j in 1:nb.h){
      coeff.i <- round(as.numeric(Coe@coeff[i, (0:3)*nb.h+j]), 8)
      nef     <- append(nef, c(coeff.i, "\n"))}
  }
  cat(nef, file=file)
  cat(".nef file succesfully written here:", file)
}

# xFourier core functions ############################################

# efourier 
efourier  <- function (coo, nb.h = 32, smooth.it = 0, silent = FALSE) {
  if (is.matrix(coo)) coo <- m2l(coo)
  if (is.closed(coo)) coo <- coo.unclose(coo)
  if (missing(nb.h))  {
    nb.h <- length(coo$x)/2 - 1 # should not be 1
    warning(paste(" * 'nb.h' not provided and set to", nb.h))}
  if(nb.h * 2 > length(coo$x)) {
    nb.h = floor(length(coo$x)/2)-1 # should not be -1
    if (!silent){
    warning(" * The number of harmonics to calculate should be lower than half the number of points. 
    'The number of harmonics used 'nb.h' has been set to: ", nb.h)}}
  if (nb.h == -1) {
    nb.h = floor(length(coo$x)/2)-1 # should not be -1
    if (!silent){
    cat(" * The number of harmonics used has been set to: ", nb.h)}}
  if (smooth.it!=0) { coo <- coo.smooth(coo, smooth.it)}
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
  for (i in 1:nb.h) {
    mat <- size * rotation %*%
      matrix(c(ef$an[i], ef$cn[i], ef$bn[i], ef$dn[i]), 2, 2) %*%
      matrix(c(cos(i*theta), sin(i*theta), -sin(i*theta), cos(i*theta)), 2, 2)
      A[i] <- mat[1, 1]
      B[i] <- mat[1, 2]
      C[i] <- mat[2, 1]
      D[i] <- mat[2, 2]
      lnef <- c(A[i], B[i], C[i], D[i])}
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

# rfourier 
rfourier <- function(coo, nb.h, smooth.it=0, norm=FALSE, silent=FALSE){
  if (missing(nb.h))  {stop("nb.h must be provided")}
  if (is.list(coo))   {coo <- l2m(coo)}
  if (is.closed(coo)) {coo <- coo.unclose(coo)}
  if(nb.h * 2 > nrow(coo)) {
    nb.h = floor(nrow(coo)/2)-1 # should not be -1
    if (!silent){
      warning("The number of harmonics to calculate should be lower than half the number of points. 
    The number of harmonics used has been set to: ", nb.h)}}
  if (nb.h == -1) {
    nb.h = floor(nrow(coo)/2)-1 # should not be -1
    if (!silent){
      cat("The number of harmonics used has been set to: ", nb.h)}}
  if (smooth.it!=0) { coo <- coo.smooth(coo, smooth.it)}
  if (norm) {
    coo   <- coo.scale(coo.center(coo))
    rsize <- mean(apply(coo, 1, function(x) sqrt(sum(x^2))))
    coo   <- coo.scale(coo, 1/rsize)}
  
  # from Claude
  p     <- nrow(coo)
  an    <- bn <- numeric(nb.h)
  Z     <- complex(real=coo[, 1], imaginary=coo[, 2])
  r     <- Mod(Z)
  angle <- Arg(Z)
  ao    <- 2*sum(r)/p
  for (i in 1:nb.h){
    an[i]<-(2/p)*sum(r * cos(i*angle))
    bn[i]<-(2/p)*sum(r * sin(i*angle))}
  list(an=an, bn=bn, ao=ao, r=r)}


rfourier.i <- function(rf, nb.h, nb.pts=300) {
  if (!all(c("an", "bn") %in% names(rf))) {
    stop("a list containing 'an' and 'bn' harmonic coefficients must be provided")}
  ao <- ifelse(is.null(rf$ao), 1, rf$ao)
  an <- rf$an
  bn <- rf$bn
  if (missing(nb.h)) {nb.h <- length(an)}
  if (nb.h > length(an)) {
    nb.h <- length(an)
    warning("nb.h cannot be higher than length(rf$an) and has been set to: ", nb.h)}
  theta <- seq(0, 2*pi, length=nb.pts)
  harm  <- matrix(NA, nrow=nb.h, ncol=nb.pts)
  for (i in 1:nb.h){
    harm[i, ]<- an[i]*cos(i*theta) + bn[i]*sin(i*theta)}
  r <- (ao/2) + apply(harm, 2, sum)
  Z <- complex(modulus=r, argument=theta)
  list(x=Re(Z), y=Im(Z), angle=theta, r=r)}

rfourier.shape <- function(an, bn, nb.h, nb.pts=80, alpha=2, plot=TRUE){
  if (missing(nb.h) &  missing(an)) nb.h <- 1
  if (missing(nb.h) & !missing(an)) nb.h <- length(an)
  if (missing(an)) an <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  if (missing(bn)) bn <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  rf  <- list(an=an, bn=bn, ao=0)
  shp <- rfourier.i(rf, nb.h=nb.h, nb.pts=nb.pts)      
  if (plot) coo.plot(shp)
  return(shp)}

# tfourier 
tfourier <- function(coo, nb.h, smooth.it=0, norm=FALSE, silent=TRUE){
  if (missing(nb.h))  {stop("nb.h must be provided")}
  if (is.list(coo))   {coo <- l2m(coo)}
  if (is.closed(coo)) {coo <- coo.unclose(coo)}
  if(nb.h * 2 > nrow(coo)) {
    nb.h = floor(nrow(coo)/2)-1 # should not be -1
    if (!silent){
      warning("The number of harmonics to calculate should be lower than half the number of points. 
    The number of harmonics used has been set to: ", nb.h)}}
  if (nb.h == -1) {
    nb.h = floor(nrow(coo)/2)-1 # should not be -1
    if (!silent){
      cat("The number of harmonics used has been set to: ", nb.h)}}
  if (smooth.it!=0) { coo <- coo.smooth(coo, smooth.it)}
  if (norm) {
    coo <- coo.scale(coo.center(coo))
    coo <- coo.trans(coo, -coo[1, 1], -coo[1, 2])
  }
 
  p <- nrow(coo)
  an <- bn <- numeric(nb.h)
  tangvect <- coo - rbind(coo[p,], coo[-p,])
  perim <- sum(sqrt(apply((tangvect)^2, 1, sum)))
  v0    <- coo[1,]-coo[p,]
  tet1   <- Arg(complex(real=tangvect[,1], imaginary = tangvect[,2]))
  tet0   <- tet1[1]
  t1     <- seq(0, 2*pi, length= (p+1))[1:p]
  phi    <- (tet1-tet0-t1)%%(2*pi)
  ao     <- 2*sum(phi)/p
  for (i in 1:nb.h){
    an[i]<- (2/p) * sum( phi * cos (i*t1))
    bn[i]<- (2/p) * sum( phi * sin (i*t1))}
  list(ao=ao, an=an, bn=bn, phi=phi, t=t1, perimeter=perim,
       thetao=tet0, x1=coo[1, 1], y1=coo[1, 2])}


tfourier.i<-function(tf, nb.h, nb.pts=300, force2close=FALSE, rescale=TRUE, perim=2*pi, thetao=0){
  if (!all(c("an", "bn") %in% names(tf))) {
    stop("a list containing 'an' and 'bn' harmonic coefficients must be provided")}
  ao <- ifelse(is.null(tf$ao), 0, tf$ao)
  if (missing(thetao)) { thetao <- ifelse(is.null(tf$thetao), 0, tf$thetao) }
  an <- tf$an
  bn <- tf$bn
  if (missing(nb.h)) {nb.h <- length(an)}
  if (nb.h > length(an)) {
    nb.h <- length(an)
    warning("nb.h cannot be higher than length(rf$an) and has been set to: ", nb.h)}
  #if (missing(nb.pts)) {nb.pts=nb.h*2}
  theta <- seq(0, 2*pi, length=nb.pts)
  harm  <- matrix(NA, nrow=nb.h, ncol=nb.pts)
  for (i in 1:nb.h){
    harm[i,] <- an[i]*cos(i*theta) + bn[i]*sin(i*theta)}
  phi  <- (ao/2) + apply(harm, 2, sum)
  vect <- matrix(NA, 2, nb.pts)
  Z    <- complex(modulus=(2*pi)/nb.pts, argument=phi+theta+thetao)
  Z1   <- cumsum(Z)
  coo <- cbind(Re(Z1), Im(Z1))
  if (force2close) { coo <- coo.force2close(coo)}
  if (rescale)     {
    if(missing(perim)) {
      perim <- ifelse(is.null(tf$perim), 2*pi, tf$perim)}
    coo <- coo.scale(coo, coo.perim(coo)/perim) }
  if (!all(is.null(tf$x1) & is.null(tf$x1))) {
    coo <- coo.trans(coo, tf$x1, tf$y1)}
  return(list(x=coo[, 1], y=coo[, 2], angle=theta, phi=phi))
}



tfourier.shape <- function(an, bn, ao=0, nb.h, nb.pts=80, alpha=2, plot=TRUE){
  if (missing(nb.h) &  missing(an)) nb.h <- 1
  if (missing(nb.h) & !missing(an)) nb.h <- length(an)
  if (missing(an)) an <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  if (missing(bn)) bn <- runif(nb.h, -pi, pi) / (1:nb.h)^alpha
  tf  <- list(an=an, bn=bn, ao=ao)
  shp <- tfourier.i(tf, nb.h=nb.h, nb.pts=nb.pts)      
  if (plot) coo.plot(shp)
  return(shp)}

# PCA and Cie ########################################################

pca2shp <- function (pos, rot, mean.shp,
                     method=c("efourier", "rfourier", "tfourier"),
                     scale=1, amp=1, trans=TRUE, nb.pts=64, rotate.shp) {
  # we check a bit
  if (!is.matrix(pos))        pos <- as.matrix(pos)
  if (ncol(pos) != ncol(rot)) stop("rot an pos must have the same ncol")
  if(length(mean.shp) != nrow(rot)) stop("mean.shp length must equals the col number of rot")
  # we handle method argument
  if (missing(method)) {
    warning("Method not provided. efourier is used.")
    p <- 1 # we also need to switch below
    method.i   <- efourier.i
  } else {
    p <- pmatch(tolower(method), c("efourier", "rfourier", "tfourier"))
    if (is.na(p)) { warning("Unvalid method. efourier is used.")
    } else {
      method.i   <- switch(p, efourier.i,   rfourier.i,   tfourier.i)}}
  
  # stupid function
  mprod <- function(m, s){
    res <- m
    for (i in 1:ncol(m)) {
      res[, i] <- m[, i]*s[i]}
    return(res)}
  nb.h <- length(mean.shp)/ifelse(p==1, 4, 2)
  n  <- nrow(pos)
  # we prepare the array
  res <- array(NA, dim=c(nb.pts, 2, n),
               dimnames=list(paste0("pt", 1:nb.pts),
                             c("x", "y"),
                             paste0("shp", 1:n)))
  for (i in 1:n) {
    ax.contrib <- mprod(rot, pos[i, ])*amp
    coe        <- mean.shp + apply(ax.contrib, 1, sum)
    if (p==1) {
      xf <- list(an=coe[1:nb.h + 0*nb.h],
                 bn=coe[1:nb.h + 1*nb.h],
                 cn=coe[1:nb.h + 2*nb.h],
                 dn=coe[1:nb.h + 3*nb.h])
    } else {
      xf <- list(an=coe[1:nb.h + 0*nb.h],
                 bn=coe[1:nb.h + 1*nb.h])
      }
    coo        <- l2m(method.i(xf, nb.h = nb.h, nb.pts=nb.pts))
    coo <- coo.template(coo, size=scale)
    # if required we rotate shapes
    if (!missing(rotate.shp)) { coo <- coo.rotate(coo, rotate.shp) }
    # by default we return translated shapes, ready to draw as a layer
    if (trans)                { coo <- coo.trans(coo, pos[i, 1], pos[i, 2]) }
    res[,,i] <- coo.force2close(coo)
  }
  invisible(res)}

# Morphological space
morpho.space <- function(dudi, xax = 1, yax = 2, xlim, ylim, nb.pts=300,
                         pos.shp=c("li", "circle", "range")[3], 
                         nr.shp=6, nc.shp=5, amp.shp=1, scale.shp=1, rotate.shp=0,
                         circle.nb.shp=12, circle.r.shp,
                         plot=TRUE, layer=TRUE, col.shp="#70809011", border.shp="#708090",
                         pch.pts=20, col.pts="grey40", first.point=FALSE){
  # we first check argument passed to pos.shp, to define pos
  if (is.data.frame(pos.shp)) pos.shp <- as.matrix(pos.shp) # e.g. when passed with expand.grid
  if (is.matrix(pos.shp)) {
    if (ncol(pos.shp)!=2) {stop("When passed with a matrix, pos.shp requires a two columns matrix")}
    pos <- pos.shp
  } else if (pos.shp=="li") { # we retrieve coordinates from the dudi.object
    pos <- dudi$li[, c(xax, yax)]
  } else if (pos.shp=="circle") {
    if (missing(circle.r.shp)) { # if missing we define it as the mean distance from the origin
      li.2      <- apply(dudi$li[, c(xax, yax)], 2, function(x) x^2)
      li.len    <- apply(li.2, 1, function(x) sqrt(sum(x)))
      circle.r.shp <- mean(li.len)}
    t <- seq(0, 2*pi, len=circle.nb.shp+1)[-(circle.nb.shp+1)]
    pos <- cbind(circle.r.shp*cos(t), circle.r.shp*sin(t))
  } else if (pos.shp=="range") { # by default, we aim at covering the range covered by points on the two PC axes
    pos <- expand.grid(seq(min(dudi$li[, xax]), max(dudi$li[, xax]), len=nr.shp),
                       seq(min(dudi$li[, yax]), max(dudi$li[, yax]), len=nc.shp))
    pos <- as.matrix(pos)
  } else {
    stop("shp.pos must be passed with values li, circle, range or a matrix of coordinates")}
  # We define a scale.shp that should fit, and then modify it with scale.shp
  if (missing(scale.shp)) {
    scale.shp <- min(apply(dudi$li[,c(xax, yax)], 2, function(x) diff(range(x)))/(c(nr.shp, nc.shp)-1))}
  # Here we (finally) calculate the shapes
  shapes <- pca2shp(pos, rot=dudi$c1[, c(xax, yax)],
                    mean.shp=dudi$mean.shp, method=dudi$method,
                    scale=scale.shp, amp=amp.shp, rotate.shp=rotate.shp, nb.pts=nb.pts)
  if (plot) {
    if (missing(xlim) & missing(ylim)) {
      w <- apply(shapes, 2, range)
    } else {
      w <- cbind(xlim, ylim)}
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mar=c(3, 3, 1, 1))
    plot(dudi$li[, c(xax, yax)], xlim=w[, 1], ylim=w[,2], asp=1, las=1,
         col=col.pts, pch=pch.pts, cex=1, cex.axis=0.7, ann=FALSE)
    abline(h=0, v=0, lty=2, col="grey80")
    box() }
  if (layer) {
    apply(shapes, 3, coo.draw, points=FALSE, border=border.shp, col=col.shp, first.point=first.point)
    #A.segments(shapes, col=col.shp, border.col=border.shp)
  }
  invisible(shapes)
}

# dudi.plot
dudi.plot <- function(dudi, fac = NULL, xax = 1, yax = 2, grid = TRUE,
                      points     = TRUE,  pch.points=1,  col.points="black", cex.points=0.8,
                      labels     = FALSE, label=rownames(dudi$li), boxes=TRUE, clabel=0.6,
                      neighbors  = FALSE, col.nei="grey90",  lwd.nei=0.5,
                      star       = TRUE,  col.star="grey60", cstar=1,
                      ellipses   = TRUE,  col.ellipse="grey30", cellipse=1, axesell=TRUE,
                      chull      = FALSE, col.chull="grey30", optchull = c(0.5, 1),
                      arrows     = FALSE, edge.arrow=FALSE, box.arrow=TRUE, maxnb.arrow=10, dratio.arrow=0.2, 
                      shapes     = TRUE,  pos.shp=c("li", "circle", "range", "full")[3],
                      nr.shp=6, nc.shp=5, amp.shp=1, scale.shp=0.666, nb.pts.shp=300, first.point.shp=FALSE, rotate.shp=0,
                      circle.nb.shp=12, circle.r.shp,
                      col.shp="#70809011", border.shp="#708090",      
                      rug        = TRUE, rug.ticksize=0.01, rug.col="#708090",
                      eigen      = FALSE, eigen.ratio=0.2,
                      palette    = col.sari,
                      title      = substitute(dudi), legend=FALSE,
                      center.orig= FALSE,
                      zoom.plot  = 1){
  
  # we prepare and check a bit
  if (!missing(fac)) {
    if (!is.factor(fac)) {
      if (ncol(dudi$fac)==0) { fac <- factor(rep("", nrow(dudi$li))) } else {fac <- dudi$fac[, fac]}}
    if ((nlevels(fac) > 1)) {
      if (missing(col.star))    col.star    <- paste(palette(nlevels(fac)), "33", sep="")
      if (missing(col.ellipse)) col.ellipse <- palette(nlevels(fac))
      if (missing(col.chull))   col.chull   <- palette(nlevels(fac))} }
  
  # we initialize the factorial map
  if (center.orig) {
    li.2      <- apply(dudi$li[,c(xax, yax)], 2, function(x) x^2)
    li.len    <- apply(li.2, 1, function(x) sqrt(sum(x)))
    xw <- max(li.len)*(1/zoom.plot)
    yw <- min(li.len)*(1/zoom.plot)
    s.label(dudi$li, xax=xax, yax=yax, xlim=c(-xw, xw), ylim=c(-yw, yw), clabel=0, cpoint=0, sub=title, grid=grid)
  } else {     
    s.label(dudi$li, xax=xax, yax=yax, clabel=0, cpoint=0, sub=title, grid=grid)}
  
  # size of the grid
  xaxp <- par("xaxp")
  ax <- (xaxp[2] - xaxp[1])/xaxp[3]
  yaxp <- par("yaxp")
  ay <- (yaxp[2] - yaxp[1])/yaxp[3]
  d <- min(ax, ay)
  
  # we redefine shorter margins
  op <- par("mar")
  par(mar=rep(0.1, 4))
  
  # rug
  if (rug) {
    rug(dudi$li[, xax], side=1, ticksize=rug.ticksize, col=rug.col, lwd=0.4)
    rug(dudi$li[, yax], side=2, ticksize=rug.ticksize, col=rug.col, lwd=0.4)
    box()}
  
  # neighbors network
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
  
  # arrow - bloody dirty below
  if (arrows) {
    arr.2      <- apply(dudi$co[,c(xax, yax)], 2, function(x) x^2)
    arr.len    <- apply(arr.2, 1, function(x) sqrt(sum(x)))
    if (maxnb.arrow > length(dudi$cw)) { maxnb.arrow <- length(dudi$cw) }
    arr.sorted <- order(arr.len, decreasing=TRUE)[1:maxnb.arrow]
    arr.disp <- if (missing(dratio.arrow)) {
      arr.len[arr.sorted] > 0
      } else {
      arr.len[arr.sorted] > d*dratio.arrow }
    if (sum(arr.disp)>0) {
      #arr.disp   <- arr.sorted[ arr.len[arr.sorted] > d*dratio.arrow ]
      arr.co <- dudi$co[names(which(arr.disp)), c(xax, yax)]
      #s.arrow(dudi$co[names(which(arr.disp)), c(xax, yax)],
      #        label = rownames(dudi$co[arr.disp, c(xax, yax)]),
      #        edge=edge.arrow, add.plot=TRUE, boxes=box.arrow)
      s.arrow(arr.co, 1, 2, 
              label = rownames(arr.co),
              edge=edge.arrow, add.plot=TRUE,boxes=box.arrow, clabel=clabel)
      
      
      }}
  
  # shapes
  if (!is.null(dudi$method)) { # for use on dudi.pca objects
  if ((dudi$method != "tFourier")) {
    if (shapes) {
    if (!is.matrix(pos.shp)) {
      if (pos.shp=="full") {
        w <- par("usr")
        pos.shp <- as.matrix(expand.grid(seq(w[1]+d/2, w[2]-d/2, len=nr.shp),
                                         seq(w[3]+d/2, w[4]-d/2, len=nc.shp)))}}
    shapes <- morpho.space(dudi, xax=xax, yax=yax, plot=FALSE, layer=TRUE,
                           nb.pts=nb.pts.shp, pos.shp=pos.shp,
                           nr.shp = nr.shp, nc.shp = nc.shp, amp.shp = 1, 
                           scale.shp = d*scale.shp, rotate.shp = rotate.shp,
                           circle.nb.shp = circle.nb.shp, circle.r.shp = circle.r.shp,                         
                           col.shp="#70809011", border.shp="#708090", first.point=first.point.shp, pch.pts=NA)}}
  }
  # labels and points
  if (points) {
  repeach <- function(x, each){ # forgot that but probbaly dirty
    if (length(x) != length(each)) return(rep(x[1], sum(each)))
    res <- vector(mode = class(x[1]))
    for (i in seq(along=x)) {
      res <- append(res, rep(x[i], each[i]))}
    return(res)}
  if (!is.null(fac)) {
    nb <- table(fac)
    if (missing(pch.points)) {
      pch.points <- repeach(pch.points, nb)}
    if (missing(col.points)) {
      #col.points <- repeach(palette(nlevels(fac)), nb)
      col.points <- palette(nlevels(fac))[fac]
      }
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
  
  #legend
  if ((legend) | (missing(legend) & !is.null(fac))) {
    legend("topright", col=palette(nlevels(fac)), lwd=2,
           legend=levels(fac), bty="n")}
  # eigen
  if (eigen) {
    par("mar"=op)
    add.scatter.eig(dudi$eig, nf=dudi$nf, xax=xax, yax=yax, eigen.ratio, posi="bottomright")}
  
  # we restore the margins
  par("mar"=op)
}

# PC contribution to shape
PC.contrib <- function(dudi, PC.r=1:dudi$nf, sd=2,
                       cols=rep(NA, 3), borders=c("#000080", "#000000", "#EE0000"),
                       lwd=1, nb.pts=300, plot=TRUE, legend=TRUE){
  # mean.coo <- efourier.i(coeff.split(dudi$mean.shp, nb.h=length(dudi$mean.shp)/4))
  if ((length(PC.r) > dudi$nf) | (max(PC.r) > dudi$nf)) {
    stop("The PC.r must correspond to PC axes present in the dudi object")}
  res <- list()
  for (i in seq(along=PC.r)) {
    pos.i <- sd*sd(dudi$li[, PC.r[i]])
    shp.i <- pca2shp(pos=matrix(c(-pos.i, 0, pos.i), nrow=3),
                     rot=as.matrix(dudi$c1[, PC.r[i]]), 
                     mean.shp=dudi$mean.shp, method=dudi$method,
                     trans=FALSE, nb.pts=nb.pts)
    shp.i <- a2l(shp.i) # we reconvert to list
    names(shp.i) <- paste0(rep(paste0("PC", PC.r[i]), 3), c("-", "m", "+"))
    res <- append(res, shp.i)}
  if (plot) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mar=c(1, 2, 1, 1), xpd=NA)
    n   <- length(PC.r)
    pos <- cbind(1:n, matrix((n+1):(4*n), nrow=n, ncol=3, byrow=TRUE))
    plot(NA, asp = 1, xlim = c(0, 4), ylim = c(0, n), 
         xaxs = "i", yaxs = "i", frame = FALSE, ann = FALSE, axes = FALSE)
    res.t <- lapply(res, coo.template, size=0.9)
    for (i in 1:n) {
      coo.draw(coo.trans(res.t[[(i-1)+1]], 0.5, n-((i-1)+0.5)),
               col=cols[1], border=borders[1], lwd=lwd, points = FALSE, first.point=FALSE)
      coo.draw(coo.trans(res.t[[(i-1)+2]], 0.5, n-((i-1)+0.5)),
               col=cols[2], border=borders[2], lwd=lwd, points = FALSE, first.point=FALSE)
      coo.draw(coo.trans(res.t[[(i-1)+3]], 0.5, n-((i-1)+0.5)),
               col=cols[3], border=borders[3], lwd=lwd, points = FALSE, first.point=FALSE)
    }
    for (i in 1:(n*3)) {
      pos.x <- rep(0:2 + 1.5, times=n)
      pos.y <- rep((n-1):0*1 + 0.5, each=3)
      coo.draw(coo.trans(res.t[[i]], pos.x[i], pos.y[i]),
               col=cols[((i-1) %% 3) +1], border=borders[((i-1) %% 3) +1],
               lwd=lwd, points=FALSE, first.point=FALSE)}
    if (legend) {
      text(1.5, n, labels=paste("-", sd, "s.d.", sep=""), adj=0.5)
      text(2.5, n, labels="Mean", adj=0.5)
      text(3.5, n, labels=paste("+", sd, "s.d.", sep=""), adj=0.5)
      text(0, (n:1) - 0.5, labels=paste("PC", PC.r), adj=1)
    }
  }
  invisible(res)}

# Thin Plate Spline ##################################################

tps2d <- function(grid0, fr, to){
  if (is.closed(fr)) fr <- coo.unclose(fr)
  if (is.closed(to)) to <- coo.unclose(to)
  p  <- nrow(fr)
  q  <- nrow(grid0)
  P  <- matrix(NA, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      r2     <- sum((fr[i,]-fr[j,])^2)
      P[i,j] <- r2*log(r2)}}
  P[is.na(P)] <- 0
  Q  <- cbind(1, fr)
  L  <- rbind(cbind(P, Q), cbind(t(Q), matrix(0,3,3)))
  m2 <- rbind(to, matrix(0, 3, 2))
  coefx <- solve(L)%*%m2[, 1]
  coefy <- solve(L)%*%m2[, 2]
  fx <- function(fr, grid0, coef) {
    Xn <- numeric(q)
    for (i in 1:q) {
      Z     <- apply((fr-matrix(grid0[i, ], p, 2, byrow=TRUE))^2, 1, sum)
      Xn[i] <- coef[p+1]+coef[p+2]*grid0[i,1]+coef[p+3]*grid0[i,2]+
        sum(coef[1:p]*(Z*log(Z)))}
    return(Xn)}
  grid1 <- cbind(fx(fr, grid0, coefx), fx(fr, grid0, coefy))
  return(grid1)}

tps.grid <- function(fr, to, amp=1, plot.full=TRUE, grid.outside = 0.2,
                     grid.size = 20, grid.col   = "grey40",
                     shp = TRUE, shp.col =  rep(NA, 2), shp.border=col.gallus(2),
                     shp.lwd = c(2, 2), shp.lty = c(1, 1)){
  # simple magnification
  if (!missing(amp)) to <- to + (to-fr)*amp
  # we prepare the grid
  x1     <- min(to[, 1])
  x2     <- max(to[, 1])
  y1     <- min(to[, 2])
  y2     <- max(to[, 2])
  rx     <- x2 - x1
  ry     <- y2 - y1
  dim.grid <- if (rx > ry) { c(grid.size, round(grid.size*ry / rx)) } else { c(round(grid.size*rx / ry), grid.size) }
  xgrid0 <- seq(x1-rx*grid.outside, x2+rx*grid.outside, length=dim.grid[1])
  ygrid0 <- seq(y1-ry*grid.outside, y2+ry*grid.outside, length=dim.grid[2])
  grid0 <- as.matrix(expand.grid(xgrid0, ygrid0))
  grid1 <- tps2d(grid0, fr, to)
  if (plot.full){
    wdw <- apply(rbind(grid0, grid1), 2, range)
  } else {
    wdw <- apply(rbind(fr, to), 2, range)}
  plot(NA, xlim=wdw[, 1], ylim=wdw[, 2], asp=1, ann=FALSE, axes=FALSE, mar=rep(0, 4))
  for (i in 1:dim.grid[2]) lines(grid1[(1:dim.grid[1]) + (i-1)*dim.grid[1],], col=grid.col)
  for (i in 1:dim.grid[1]) lines(grid1[(1:dim.grid[2]) * dim.grid[1]-i+1,],   col=grid.col)
  if (shp) {
    coo.draw(fr, border=shp.border[1], col=shp.col[1], lwd=shp.lwd[1], lty=shp.lty[1])
    coo.draw(to, border=shp.border[2], col=shp.col[2], lwd=shp.lwd[2], lty=shp.lty[2])}
}


tps.arr <- function(fr, to, amp=1, palette = col.summer,
                    arr.nb = 100, arr.levels = 100, arr.len = 0.1,
                    arr.ang = 30, arr.lwd = 1, arr.col = "grey50",
                    shp = TRUE, shp.col =  rep(NA, 2), shp.border=col.gallus(2),
                    shp.lwd = c(2, 2), shp.lty = c(1, 1)){
  if (!missing(amp)) to <- to + (to-fr)*amp
  grid0  <- spsample(Polygon(coo.close(fr)), arr.nb, type="regular")@coords
  grid1     <- tps2d(grid0, fr, to)
  # grille simple, on affiche d'abord les deux courbes
  wdw      <- apply(rbind(fr, to), 2, range)
  plot(NA, xlim=wdw[, 1]*1.05, ylim=wdw[, 2]*1.05, asp=1, axes=FALSE, ann=FALSE, mar=rep(0,4))
  if (missing(arr.levels)) {arr.levels = arr.nb}
  if (!missing(palette)) {
    q.lev   <- cut(edm(grid0, grid1), breaks=arr.levels, labels=FALSE)
    arr.cols <- palette(arr.levels)[q.lev]
  } else {
    arr.cols <- rep(arr.col, nrow(grid0))}
  arrows(grid0[, 1], grid0[, 2], grid1[, 1], grid1[, 2],
         length=arr.len, angle=arr.ang, lwd=arr.lwd, col=arr.cols)
  if (shp) {
    coo.draw(fr, border=shp.border[1], col=shp.col[1], lwd=shp.lwd[1], lty=shp.lty[1])
    coo.draw(to, border=shp.border[2], col=shp.col[2], lwd=shp.lwd[2], lty=shp.lty[2])}
}


tps.iso <- function(fr, to, amp=1, palette = col.summer,
                    iso.nb = 500, iso.levels = 12, cont=TRUE, cont.col="black",
                    shp = TRUE, shp.col =  rep(NA, 2), shp.border=col.gallus(2),
                    shp.lwd = c(2, 2), shp.lty = c(1, 1)){  
  if (!missing(amp)) to <- to + (to-fr)*amp
  grid0  <- spsample(Polygon(coo.close(fr)), iso.nb, type="regular")@coords
  grid1  <- tps2d(grid0, fr, to)
  def    <- edm(grid0, grid1)
  x1     <- length(unique(grid0[,1]))
  y1     <- length(unique(grid0[,2]))
  im     <- matrix(NA,x1,y1)
  xind   <- (1:x1)[as.factor(rank(grid0[,1]))]
  yind   <- (1:y1)[as.factor(rank(grid0[,2]))]
  n      <- length(xind)
  for (i in 1:n) im[xind[i], yind[i]] <- def[i]
  iso.cols <- palette(iso.levels)
  x <- sort(unique(grid0[,1]))
  y <- sort(unique(grid0[,2]))
  image(x, y, im, col=iso.cols, asp=1, xlim=range(x)*1.05, ylim=range(y)*1.05,
        axes=FALSE, frame=FALSE, ann=FALSE)
  if (cont) contour(x, y, im, nlevels=iso.levels, add=TRUE, drawlabels=FALSE, col=cont.col)
  if (shp) {
    coo.draw(fr, border=shp.border[1], col=shp.col[1], lwd=shp.lwd[1], lty=shp.lty[1])
    coo.draw(to, border=shp.border[2], col=shp.col[2], lwd=shp.lwd[2], lty=shp.lty[2])}}

# Miscellaneous ######################################################

ellpar <- function(coo){
  if (is.list(coo)) coo <- cbind(coo$x, coo$y) ### ML
  coo <- coo %*% svd(var(coo))$u
  a <- max(coo[, 1])
  b <- max(coo[, 2])
  e <- sqrt((a^2 - b^2)/a^2)
  list(a=a, b=b, e=e)}

dev.plot       <- function(mat, dev, cols, x=1:ncol(mat), 
                           lines=TRUE, poly=TRUE, segments=FALSE, bw=0.1,
                           plot=FALSE, main="Deviation plot", xlab="", ylab="Deviations") {
  # we prepare and check a bit
  r <- nrow(mat)
  if (!missing(dev)){
    if (any(dim(mat)!=dim(dev))) {
      stop("mat and dev must be of the same dimension")}}
  if (missing(cols))   {cols <- rep("#000000", r)}
  if (length(cols)!=r) {cols <- rep("#000000", r)}
  # we call a new plot if required
  if (plot) {
    if (missing(dev)) {
      ylim <- range(mat)
      } else { ylim <- c(min(mat+dev), max(mat+dev)) }
    plot(NA, xlim=range(x), ylim=ylim, main=main, xlab=xlab, ylab=ylab, las=1, xaxs="i", yaxs="i")
    axis(1, at=1:ncol(mat))}
  # if a deviation matrix is provided
  if (!missing(dev)){
    for (i in 1:r){
      # if required, we draw the background polygons
      if (poly) {
        polygon(x = c(x, rev(x)),
                y = c(mat[i, ] - dev[i, ], rev(c(mat[i, ] + dev[i, ]))),
                col=paste0(cols[i], "55"), border=NA)}
      # if required we draw the dev segments
      if (segments) {
      segments(x,    mat[i, ] - dev[i, ], x, mat[i, ]    + dev[i, ], col=cols[i], lwd=0.5)
      segments(x-bw, mat[i, ] - dev[i, ], x+bw, mat[i, ] - dev[i, ], col=cols[i], lwd=0.5)
      segments(x-bw, mat[i, ] + dev[i, ], x+bw, mat[i, ] + dev[i, ], col=cols[i], lwd=0.5)}}}
  # if a dev matrix is not provided, we simply draw lines
  if (lines) {
    for (i in 1:nrow(mat)) {
      if (lines) {
      lines(x, mat[i, ], col=cols[i], type="o", cex=0.25, pch=20)}}}}

# returns difference angle and norm ratios between two vectors given as 4 numeric.
vecs.param <- function(r1, i1, r2, i2){
  x <- c(r1, i1, r2, i2)
  if (!is.numeric(x)) {stop("4 numeric must be passed.")}
  if (length(x)!=4)   {stop("4 numeric must be passed.")}
  r.norms <- sqrt((r2^2 + i2^2)) / sqrt((r1^2 + i1^2))
  d1 <- sqrt(sum(r1^2 + i1^2))
  d2 <- sqrt(sum(r2^2 + i2^2))
  return(list(r.norms=d1/d2, d.angle=atan2(i2, r2) - atan2(i1, r1)))}

# Calculates harmonic power given a list from e/t/rfourier
harm.pow <- function(xf){
  if (is.list(xf)) {
    if (all(c("an", "bn", "cn", "dn") %in% names(xf))) {
      return((xf$an^2 + xf$bn^2 + xf$cn^2 + xf$dn^2)/2)
    } else {
      if (all(c("an", "bn") %in% names(xf))) {
        return((xf$an^2 + xf$bn^2)/2)}
    }
  } else {
    stop("a list containing 'an', 'bn' ('cn', 'dn') harmonic coefficients must be provided")}}

# Color palettes
col.summer <- colorRampPalette(c("#4876FF", "#FFFF00", "#FF3030"))
col.gallus <- colorRampPalette(c("#000080", "#FFFFFF", "#EE0000"))
col.blackgallus <- colorRampPalette(c("#000080", "#000000", "#EE0000"))
col.sari   <- colorRampPalette(c("#551A8B", "#FF7F00"))
col.india  <- colorRampPalette(c("#FF9933", "#138808"))
col.bw     <- colorRampPalette(c("#FFFFFF", "#000000"))
col.wcol   <- function(col.hex) colorRampPalette(c("#FFFFFF", col.hex))
col.bcol   <- function(col.hex) colorRampPalette(c("#000000", col.hex))
              