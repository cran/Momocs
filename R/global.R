setGeneric(name= "get.Nef", 	def=function(Coo, nb.h=32, smooth.it=0, fromrt=FALSE)
											{standardGeneric("get.Nef")})
setGeneric(name= "dev.qual", 	def=function(Coo, id=1:length(Coo@coo), nb.h = 32, smooth.it=0, range= seq(1, nb.h, len=4))
											{standardGeneric("dev.qual")})  
setGeneric(name= "dev.quant", 	def=function(Coo, id=1:length(Coo@coo), nb.h = 32, smooth.it=0, plot=TRUE)
											{standardGeneric("dev.quant")})
setGeneric(name= "harm.pow", 	def=function(Coo, id, nb.h = 32, smooth.it = 0, plot = TRUE, max.h = nb.h, first = TRUE)
											{standardGeneric("harm.pow")})
setGeneric(name= "manova.nef", 	def=function(Nef, fac, harmonics.retained, drop=FALSE)
											{standardGeneric("manova.nef")})
setGeneric(name= "pca", 		def=function(Nef, fac = NA,PCa = 1, PCb = 2, col = "black", pch = 1, lty=1, shp.nb=NA, shp.size, shp.col="#00000022",
											 shp.border="black", title = "Principal Component Analysis", legend = TRUE, lab = FALSE, lab.txt = rownames(Nef@coeff),
											 lab.cex = 1, lab.box = TRUE, ell = TRUE, r = 1, lwd = 1, zoom.x = 0.25, zoom.y = 0.3)
											{standardGeneric("pca")})
setGeneric(name= "pca3", 		def=function(Nef,  fac = NA, col = 1:nlevels(fac), pch = 1:nlevels(fac), lty = rep(1,nlevels(fac)),
											 lab = FALSE, lab.txt = rownames(Nef@coeff), lab.cex = 1, lab.box = TRUE, ell = 1, r = 1,
											 lwd = 1, zoom = 1.4, legend = FALSE){standardGeneric("pca3")})
setGeneric(name= "pca.tps", 	def=function(Nef, fac = NA, PCa = 1, PCb = 2, col = "black", pch = 1, ell = TRUE, zoom = 1.4, 
											 ncells = 20, title = "Deformations alongs PC axes")
											 {standardGeneric("pca.tps")})
setGeneric(name= "morph.sp",    def=function(Nef, PCa = 1, PCb = 2, nb.PCa = 5, nb.PCb = 6, fac = NA, morph.sp.extend = 1, zoom.extend = 1.2, asp,
											pch = 20, shp.col = NA, shp.lwd = 1, shp.size,
											col = "grey40", ell = FALSE, r = 1, lwd = 1, title = "Morphological space")
											{standardGeneric("morph.sp")})
setGeneric(name= "morph.PC", 	def=function(Nef, sd.nb=1, pca.ax=seq(1, 3))
											{standardGeneric("morph.PC")})
setGeneric(name= "traj", 		def=function(Nef, fr = c(0, 0), to = c(1, 1), nb.int = 50,
											 nb.pts=500, save = FALSE, prog = TRUE, pause=TRUE)
											{standardGeneric("traj")})
setGeneric(name= "tps.grid", 	def=function(Nef, fr, to, nb.pts = 50, amp = 1, grid.size = 50, grid.col = "grey40", cont = TRUE,
											 cont.col = c("dodgerblue3", "firebrick3"), cont.lwd = rep(3,2))
											 {standardGeneric("tps.grid")})
setGeneric(name= "tps.iso", 	def=function(Nef, fr, to, nb.pts = 200, amp = 1, iso.pts = 1000, col.pal = topo.colors, col.lev = 500, cont.to = TRUE, cont.fr = TRUE,
											 cont.lev = 10, cont.col = c("dodgerblue3", "firebrick3"),  cont.lwd = rep(3,2))
											{standardGeneric("tps.iso")})
setGeneric(name= "tps.vf", 		def=function(Nef, fr, to, nb.pts = 100, amp = 1, arr.nb = 300, arr.len = 0.05, arr.ang = 30, arr.col = "grey40",
											 arr.pal = FALSE, arr.palette = topo.colors, arr.pal.lev = 20, arr.lwd = 1, cont.col = c("dodgerblue3", "firebrick3"),
											 cont.lwd = rep(3, 2))
											 {standardGeneric("tps.vf")})

############################# Internal functions ###############################
# Import .txt or .jpg outlines
get.cont <- function (path) {
  if (missing(path)) {path <- choose.dir()}
  img.list <- list.files(path, full.names = TRUE)
  img.lite <- list.files(path, full.names = FALSE)
  ext <- substr(img.list, nchar(img.list) - 2, nchar(img.list))
  ext.types <- unique(ext)
  if(any(ext.types=="jpg") & any(ext.types=="txt")){
    stop("A folder with only .jpg or only .txt files must be provided")}
  if (any(ext == "jpg")) {
    retain <- which(ext=="jpg")
    img.list <- img.list[retain]
    img.lite <- img.lite[retain]
    coo <- list()
    #t0 <- Sys.time()
    for (i in seq(along = img.list)) {
      #t1 <- Sys.time()
      cat("\nOutline ", img.lite[i], "\t\t", sep="")
      img <- read.jpeg(img.list[i])
      if (class(img)[2] == "array") {img <- rgb2grey(img)}
      img[img >  0.5] <- 1
      img[img <= 0.5] <- 0
      img <- imagematrix(img)
      if(any(dim(img)>500)) {cat("\t(large image)")}
      if(any(c(img[1, ], img[nrow(img), ], img[, 1], img[, ncol(img)]) != 1)){
        cat("\t(outline spans image border)")
        img <- rbind(rep(1, ncol(img)), img, rep(1, ncol(img)))
        img <- cbind(rep(1, nrow(img)), img, rep(1, nrow(img)))
        img <- imagematrix(img)}              
      cont <- Conte(img)
      coo[i] <- list(cont)
      #t2 <- Sys.time()
      #cat("\tDone in", format(t2-t1, digits=3),"\t")
      #cat("~", format((t2-t0)*(length(img.list)/i), digits=3),
      #    " remaining", sep="")
	  }
  } else { # if .txt files and not .jpg
    retain <- which(ext=="txt")
    img.list <- img.list[retain]
    img.lite <- img.lite[retain]
    coo <- list()
    #t0 <- Sys.time()
    for (i in seq(along = img.list)) {
      cont <- read.table(img.list[i], header = TRUE)
      coo[i] <- list(list(x = cont[, 1], y = cont[, 2]))
    }
  }
  names(coo) <- substr(img.lite, 1, nchar(img.lite) - 4)
  cat("\n")
  return(Coo(coo))}

# Core function that calculates Elliptical Fourier Analysis
efourier <- function (coo, nb.h = 32, smooth.it = NULL) {
  if (!is.list(coo)) 
      stop("A list of coordinates (x,y) must be provided")
	if (is.numeric(smooth.it)) {coo <- cont.smooth(coo, smooth.it)}
  coo <- cont.sample(coo, nb.h * 2)
  p <- length(coo$x)
  Dx <- coo$x - coo$x[c(p, (1:p - 1))]
  Dy <- coo$y - coo$y[c(p, (1:p - 1))]
  Dt <- sqrt(Dx^2 + Dy^2)
  t1 <- cumsum(Dt)
  t1m1 <- c(0, t1[-p])
  T <- sum(Dt)
  an <- bn <- cn <- dn <- numeric(nb.h)
  for (i in 1:nb.h) {
      an[i] <- (T/(2 * pi^2 * i^2)) * sum((Dx/Dt) * 
          (cos(2 * i * pi * t1/T) - cos(2 * pi * i * t1m1/T)))
      bn[i] <- (T/(2 * pi^2 * i^2)) * sum((Dx/Dt) *
          (sin(2 * i * pi * t1/T) - sin(2 * pi * i * t1m1/T)))
      cn[i] <- (T/(2 * pi^2 * i^2)) * sum((Dy/Dt) * 
          (cos(2 * i * pi * t1/T) - cos(2 * pi * i * t1m1/T)))
      dn[i] <- (T/(2 * pi^2 * i^2)) * sum((Dy/Dt) *
          (sin(2 * i * pi * t1/T) - sin(2 * pi * i * t1m1/T)))}
  ao <- 2 * sum(coo$x * Dt/T)
  co <- 2 * sum(coo$y * Dt/T)
  return(list(an = an, bn = bn, cn = cn, dn = dn, ao = ao, co = co))}

# Normalized elliptical Fourier analysis that returns a Nef object
eFa <- function(coo, nb.h = 32, smooth.it = 0, fromrt = FALSE) {
  ef <- efourier(coo, nb.h = nb.h, smooth.it = smooth.it)
  A1 <- ef$an[1]
  B1 <- ef$bn[1]
  C1 <- ef$cn[1]
  D1 <- ef$dn[1]
  theta <- 0.5 * atan(2 * (A1 * B1 + C1 * D1)/(A1^2 + C1^2 - B1^2 - D1^2)) %% pi
  phaseshift <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  M2 <- matrix(c(A1, C1, B1, D1), 2, 2) %*% phaseshift
  v <- apply(M2^2, 2, sum)
  if (v[1] < v[2]) {theta <- theta + pi/2}
  theta <- (theta + pi/2)%%pi - pi/2
  Aa <- A1 * cos(theta) + B1 * sin(theta)
  Cc <- C1 * cos(theta) + D1 * sin(theta)
  scale <- sqrt(Aa^2 + Cc^2)
  psi   <- atan(Cc/Aa)%%pi
  if (Aa<0){psi<-psi+pi}
  size  <- 1/scale
  rotation <- matrix(c(cos(psi), -sin(psi), sin(psi), cos(psi)), 2, 2)
  A <- B <- C <- D <- numeric(nb.h)
  if (fromrt) {theta <- 0}
  lnef <- NULL
  for (i in 1:nb.h) {
      mat <- size * rotation %*% matrix(c(ef$an[i], ef$cn[i], 
          ef$bn[i], ef$dn[i]), 2, 2) %*% matrix(c(cos(i * theta), 
          sin(i * theta), -sin(i * theta), cos(i * theta)), 
          2, 2)
      A[i] <- mat[1, 1]
      B[i] <- mat[1, 2]
      C[i] <- mat[2, 1]
      D[i] <- mat[2, 2]
      lnef <- c(lnef, c(A[i], B[i], C[i], D[i]))}
  list(A = A, B = B, C = C, D = D, size = scale, theta = theta, 
      psi = psi, ao = ef$ao, co = ef$co, lnef = lnef)}  

# Sample points along an outline
cont.sample <- function (coo, n) {
  if (!is.list(coo)) 
      stop("A list of coordinates must be provided")
  x <- coo$x[seq(1, length(coo$x), len = n + 1)][-(n + 1)]
  y <- coo$y[seq(1, length(coo$y), len = n + 1)][-(n + 1)]
  new.coo <- list(x = x, y = y)
  return(new.coo)}

# Simplistic smoothing algorithm
cont.smooth <- function (M, n) {
    if (is.list(M)) {M <- matrix(c(M$x, M$y), ncol = 2)}
    p <- ncol(M)
    a <- 0
    while (a <= n) {
      a <- a + 1
      Ms <- rbind(M[p, ], M[-p, ])
      Mi <- rbind(M[-1, ], M[1, ])
      M <- M/2 + Ms/4 + Mi/4}
    return(list(x = M[, 1], y = M[, 2]))}

# select colums that correspond to harmonics we want to retain
col.sel <- function (h.fr = 1, h.to = 8, h.max = 32, drop = FALSE){
    fr <- h.max * 0:3 + h.fr
    to <- fr + h.to - 1
    cs <- c(fr[1]:to[1], fr[2]:to[2], fr[3]:to[3], fr[4]:to[4])
    if (drop) {cs <- cs[-c(1, h.max + 2, 2 * h.max + 2)]}
    return(cs)}

# Close an outline
closed.outline <- function (cont) {
  if (is.list(cont)) {
    cont.c <- list(x = c(cont$x, cont$x[1]), y = c(cont$y, cont$y[1]))
  } else {
    stop("A list of (x;y) coordinates must be provided")}
  return(cont.c)}

# Algorithm for outline extraction
Conte <- 
  function (img) 
  {
    if (class(img)[1] != "imagematrix"){ 
      stop("An 'imagematrix' object is expected")}
    x <- round(dim(img)/2)
    while(img[x[1], x[2]] != 0) {
      plot(img, main = "Click a point within the shape")
      cat("\t(click)")
      rect(0, 0, ncol(img), nrow(img), border="red")
      click <- lapply(locator(1), round)
      x     <- c(nrow(img)-click$y, click$x)
      if(any(x > dim(img))) {x <- round(dim(img)/2)}
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
      if (abs(img[x[1] + M[1, S + 1], x[2] + M[2, S + 1]] - img[x[1], 
                                                                x[2]]) < 0.1) {
        a <- a + 1
        X[a] <- x[1]
        Y[a] <- x[2]
        x <- x + M[, S + 1]
        SS[a] <- S + 1
        S <- (S + 7)%%8
      }
      else if (abs(img[x[1] + M[1, S + 2], x[2] + M[2, S + 2]] - 
        img[x[1], x[2]]) < 0.1) {
        a <- a + 1
        X[a] <- x[1]
        Y[a] <- x[2]
        x <- x + M[, S + 2]
        SS[a] <- S + 2
        S <- (S + 7)%%8
      }
      else if (abs(img[x[1] + M[1, S + 3], x[2] + M[2, S + 3]] - 
        img[x[1], x[2]]) < 0.1) {
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
    return(list(x = (Y[-1]), y = ((dim(img)[1] - X))[-1]))
  }

# Draw a Fourier Ellipse
draw.Fell <- function (an = pi, bn = -pi, cn = pi, dn = pi,
  					n = 200, cols = topo.colors, title = FALSE){
    theta <- seq(0, 2 * pi, length = n + 1)[-(n + 1)]
    ell <- list(x = an * cos(theta) + bn * sin(theta), y = cn * 
        cos(theta) + dn * sin(theta))
    plot(ell, asp = 1, col = cols(length(theta)), pch = 20, xlab = paste("x = ", 
        signif(an, 3), "*cos(2.n.pi/T)+", signif(bn, 3), "*sin(2.n.pi)", 
        sep = ""), ylab = paste("y = ", signif(cn, 3), "*cos(2.n.pi/T)+", 
        signif(dn, 3), "*sin(2.n.pi)", sep = ""))
    if (title) {
        title(paste("an=", round(an, 3), ", bn=", round(bn, 3), 
            ", cn=", round(cn, 3), ", dn=", round(dn, 3), sep = ""))}
    arrows(ell$x[1], ell$y[1], ell$x[2], ell$y[2])
    return(ell)}

# Inverse EFT
iefourier <- function (an, bn, cn, dn, k, n, ao = 0, co = 0) {
    theta <- seq(0, 2 * pi, length = n + 1)[-(n + 1)]
    harmx <- matrix(NA, k, n)
    harmy <- matrix(NA, k, n)
    for (i in 1:k) {
        harmx[i, ] <- an[i] * cos(i * theta) + bn[i] * sin(i * theta)
        harmy[i, ] <- cn[i] * cos(i * theta) + dn[i] * sin(i * theta)}
    x <- (ao/2) + apply(harmx, 2, sum)
    y <- (co/2) + apply(harmy, 2, sum)
    list(x = x, y = y)}

# Draw confidence ellipses given a nuage de points    
panel.lm <- function (x, y, r = 1, col = "black", lwd = 1, lty = 1) {
    dfn <- 2
    dfd <- length(x) - 1
    shape <- var(cbind(x, y), na.rm = TRUE)
    keep <- (!is.na(x) & !is.na(y))
    center <- c(mean(x[keep]), mean(y[keep]))
    radius <- sqrt(dfn * qf(0.68, dfn, dfd)) * r
    segments <- 75
    angles <- seq(0, 2 * pi, length = segments)
    unit.circle <- cbind(cos(angles), sin(angles))
    ellipse.pts <- t(center + radius * t(unit.circle %*% chol(shape)))
    ellx <- ellipse.pts[, 1]
    elly <- ellipse.pts[, 2]
    usr <- par()$usr
    minx <- usr[1]
    maxx <- usr[2]
    miny <- usr[3]
    maxy <- usr[4]
    ellx <- ifelse(ellx < minx, minx, ellx)
    ellx <- ifelse(ellx > maxx, maxx, ellx)
    elly <- ifelse(elly < miny, miny, elly)
    elly <- ifelse(elly > maxy, maxy, elly)
    lines(ellx, elly, col = col, lwd = lwd, lty = lty)}

# Plots a shape from the PCA
pca2shp <- function (pc1 = 0, pc2 = 0, data, nb.h = ncol(data)/4, nb.pts = 500, 
    amp = 1, col = "black", lwd = 2, plot = TRUE) {
    pca <- prcomp(data, center = TRUE, scale = FALSE)
    shp.mu <- as.numeric(apply(data, 2, mean))
    coe <- shp.mu + amp * pc1 * pca$rotation[, 1] + amp * pc2 * 
        pca$rotation[, 2]
    st <- 1 + nb.h * 0:3
    end <- nb.h + nb.h * 0:3
    cont <- iefourier(coe[st[1]:end[1]], coe[st[2]:end[2]], coe[st[3]:end[3]], 
        coe[st[4]:end[4]], nb.h, nb.pts)
    if (plot) {
        plot(closed.outline(cont), asp = 1, typ = "l", ax = FALSE, 
            ann = FALSE, lwd = lwd, col = col)}
    return(cont)}
	
# Core tps
tps <- function (matr, matt, n, plot = TRUE, col = "black") {
    xm <- min(matt[, 1])
    ym <- min(matt[, 2])
    xM <- max(matt[, 1])
    yM <- max(matt[, 2])
    rX <- xM - xm
    rY <- yM - ym
    a <- seq(xm - 1/5 * rX, xM + 1/5 * rX, length = n)
    b <- seq(ym - 1/5 * rX, yM + 1/5 * rX, by = (xM - xm) * 7/(5 * 
        (n - 1)))
    m <- round(0.5 + (n - 1) * (2/5 * rX + yM - ym)/(2/5 * rX + 
        xM - xm))
    M <- as.matrix(expand.grid(a, b))
    ngrid <- tps2d(M, matr, matt)
    if (plot) {
        plot(ngrid, cex = 0.2, asp = 1, axes = FALSE, ann = FALSE, 
            mar = rep(0, 4))
        for (i in 1:m) lines(ngrid[(1:n) + (i - 1) * n, ], col = col)
        for (i in 1:n) lines(ngrid[(1:m) * n - i + 1, ], col = col)
    }
    return(ngrid)}
	
# tps grids
tps2d <- function (M, matr, matt) {
    p <- dim(matr)[1]
    q <- dim(M)[1]
    n1 <- p + 3
    P <- matrix(NA, p, p)
    for (i in 1:p) {
        for (j in 1:p) {
            r2 <- sum((matr[i, ] - matr[j, ])^2)
            P[i, j] <- r2 * log(r2)
        }
    }
    P[which(is.na(P))] <- 0
    Q <- cbind(1, matr)
    L <- rbind(cbind(P, Q), cbind(t(Q), matrix(0, 3, 3)))
    m2 <- rbind(matt, matrix(0, 3, 2))
    coefx <- solve(L) %*% m2[, 1]
    coefy <- solve(L) %*% m2[, 2]
    fx <- function(matr, M, coef) {
        Xn <- numeric(q)
        for (i in 1:q) {
            Z <- apply((matr - matrix(M[i, ], p, 2, byrow = TRUE))^2, 
                1, sum)
            Xn[i] <- coef[p + 1] + coef[p + 2] * M[i, 1] + coef[p + 
                3] * M[i, 2] + sum(coef[1:p] * (Z * log(Z)))
        }
        Xn
    }
    matg <- matrix(NA, q, 2)
    matg[, 1] <- fx(matr, M, coefx)
    matg[, 2] <- fx(matr, M, coefy)
    matg}
	