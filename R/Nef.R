############################# Nef methods ##############################################

# NEF class and builder
setClass("Nef", representation(coeff="matrix"))
Nef <- function(mat) {new(Class="Nef", coeff=mat)}

# show(Nef)
setMethod(f="show", signature="Nef", definition=function(object){
cat(
rep("*", 23, sep=""), "*** A Nef class object ", rep("*", 23, sep=""), "\n",
" ->  It contains ", ncol(object@coeff), " (", ncol(object@coeff)/4, " harmonics) normalized harmonic coefficients
     for the  ", nrow(object@coeff), " shapes and a sample of it is shown below.\n",
" ->  Various statistical and graphical methods can be passed to it.\n",
" ->  The matrix of harmonic coefficients is stored in the @coeff slot.\n",
" ->  See ?Nef for further details.\n",
rep("*", 67, sep=""), "\n", sep="")
print(signif(object@coeff[c(1:3, (nrow(object@coeff)-2):nrow(object@coeff)),
                                col.sel(1, 3, ncol(object@coeff)/4)],3))})
				
# Perfoms a MANOVA on a Nef
setMethod(f="manova.nef", signature="Nef", definition=
  function (Nef,
            fac,
            harmonics.retained,
            drop=FALSE){
  if (missing(harmonics.retained)) {
    harmonics.retained <- floor((nrow(Nef@coeff)/4)-2)
  } else if(harmonics.retained>(4*floor((nrow(Nef@coeff)/4)-2))) {
    cat("The number of harmonics to retained is too high\n")
    harmonics.retained <- floor((nrow(Nef@coeff)/4)-2)}
  if(length(fac)!=nrow(Nef@coeff)) {
  stop("The length of the factor provided (", length(fac),") provided is different than the
      number of outlines (", nrow(Nef@coeff), ") in the Nef object provided")}
  harm.sel <- col.sel(1, harmonics.retained, ncol(Nef@coeff)/4, drop=drop)
  print(summary(manova(Nef@coeff[,harm.sel]~fac), test="Hotelling"))
  cat(harmonics.retained, "harmonics are retained\n")})
  
# Plots PCA
setMethod(f="pca", signature="Nef", definition=
  function(Nef,
           fac = NA,
           PCa = 1,
           PCb = 2,
           col = "black",
           pch = 1,
           lty=1,
           shp.nb=NA,
           shp.size,
           shp.col="#00000022",
           shp.border="black",
           title = "Principal Component Analysis",
           legend = TRUE,
           lab = FALSE,
           lab.txt = rownames(Nef@coeff),
           lab.cex = 1,
           lab.box = TRUE, 
           ell = TRUE,
           r = 1,
           lwd = 1,
           zoom.x = 0.25,
           zoom.y = 0.3){
  zoom.wdw <- function (x, y, zoom.x, zoom.y = zoom.x) {
    dx <- (max(x) - min(x)) * zoom.x
    dy <- (max(y) - min(y)) * zoom.y
    min.x <- min(x) - dx
    max.x <- max(x) + dx
    min.y <- min(y) - dy
    max.y <- max(y) + dy
    list(x = c(min.x, max.x), y = c(min.y, max.y))}
  def.par <- par(no.readonly = TRUE)
  pca <- prcomp(Nef@coeff)
  proj <- list(x = pca$x[, PCa], y = pca$x[, PCb])
  wdw <- zoom.wdw(proj$x, proj$y, zoom.x, zoom.y)
  plot(NA, xlim = wdw$x, ylim = wdw$y, main = title,
       xlab = paste("PC", PCa, sep = ""),
       ylab = paste("PC", PCb, sep = ""))
  abline(h = 0, v = 0, col = "grey90")
  box()
  if (lab) {
    for (i in seq(along = lab.txt)) {
      mx <- diff(wdw$x) * 0.005
      my <- diff(wdw$y) * 0.005
      x0 <- proj$x[i]
      y0 <- proj$y[i]
      k  <- lab.txt[i]
      x1 <- x0 + 2 * mx
      y1 <- y0 + 2 * my
      x2 <- x1 + strwidth(k)
      y2 <- y1 + strheight(k)
      rect(x1 - mx, y1 - my, x2 + mx, y2 + my, col = "white", 
            border = if(lab.box) NULL else NA)
      segments(x0, y0, x1 - mx, y1 - my)
      text(x1, y1, k, adj = rep(0, 2), cex=lab.cex)}}
    if (is.factor(fac)) {
        k <- levels(fac)
        if (missing(col) | length(col)!=nlevels(fac)) {col <- 1:nlevels(fac)}
        if (missing(pch) | length(pch)!=length(pch))  {pch <- 1:nlevels(fac)}
        if (missing(lty) | length(lty)!=length(lty))  {lty <- rep(lty, nlevels(fac))}
        for (i in seq(along = k)) {
            s <- k == k[i]
            points(proj$x[s], proj$y[s], pch = pch[i], col = col[i])
            if (ell) {
                panel.lm(proj$x[s], proj$y[s], col = col[i], lty=lty[i], 
                          r = r, lwd = lwd)}
    }
    if (legend) {
      legend(wdw$x[2] * 0.75, wdw$y[2] * 0.95, legend = k, 
      pch = pch, col = col)}
  } else (points(proj, pch = pch, col = col))
  if (is.numeric(shp.nb)) {
    cat("Click on the plot to display the", shp.nb, "shapes\n")
    if (missing(shp.size)) {
        w <- par("usr")
        cont <- pca2shp(0, 0, data = Nef@coeff, plot = 0)
        range.cont <- lapply(cont, function(x) diff(range(x)))
        shp.size <- ifelse(range.cont$x >= range.cont$y, diff(w[1:2]/(9 * 
            range.cont$x)), diff(w[3:4]/(6 * range.cont$y)))}
    for (i in 1:round(shp.nb)) {
        q <- locator(1)
        cont <- pca2shp(q$x, q$y, data = Nef@coeff, plot = 0)
        cont$x <- cont$x * shp.size + q$x
        cont$y <- cont$y * shp.size + q$y
        cont <- closed.outline(cont)
        polygon(cont, col = shp.col, border=shp.border)}}})

# Plots a triple PCA and the cum axis variance
setMethod(f="pca3", signature="Nef", definition=
    function(Nef,
             fac = NA,
             col = 1:nlevels(fac),
             pch = 1:nlevels(fac),
             lty = rep(1,nlevels(fac)),
             lab = FALSE,
             lab.txt = rownames(Nef@coeff),
             lab.cex = 1,
             lab.box = TRUE,
             ell = 1,
             r = 1,
             lwd = 1,
             zoom = 1.4,
             legend = FALSE){
  def.par <- par(no.readonly = TRUE)
  pca <- prcomp(Nef@coeff)
  ev <- pca$sdev^2
  ev <- ev/sum(ev)
  layout(matrix(1:4, 2, 2))
  pca(Nef, fac, PCa = 3, PCb = 2, title = "PC3-PC2", 
      col = col, pch = pch, lty=lty, lab = lab, lab.txt = lab.txt, lab.cex = lab.cex, 
      lab.box = lab.box, ell = ell, r = r, 
      lwd = lwd, zoom.x = zoom, leg = FALSE)
  pca(Nef, fac, PCa = 1, PCb = 2, title = "PC1-PC2", 
      col = col, pch = pch, lty=lty, lab = lab, lab.txt = lab.txt, lab.cex = lab.cex, 
      ell = ell, r = r, lwd = lwd, zoom.x = zoom, leg = FALSE)
  barplot(ev[1:12], col = c(rep("grey40", 3), rep("grey90", 
      9)))
  if (legend) {
      if (missing(fac)) 
          stop("Legend cannot be plotted if fac is not provided")
      w <- par("usr")
      legend(w[2] * 0.66, w[4] * 0.95, legend = levels(fac), 
          pch = pch, col = col, bty = "n")
  }
  pca(Nef, fac, PCa = 1, PCb = 3, title = "PC1-PC3", 
      col = col, pch = pch, lty=lty, lab = lab, lab.txt = lab.txt, lab.cex = lab.cex, 
      leg = FALSE, lab.box = lab.box, ell = ell,  
      r = r, lwd = lwd, zoom.x = zoom)
  par(def.par)})
 
# Plots a PCA plus the grid deformations
setMethod(f="pca.tps", signature="Nef", definition=
  function(Nef,
           fac = NA,
           PCa = 1,
           PCb = 2,
           col = "black",
           pch = 1,
           ell = TRUE,
           zoom = 1.4,
           ncells = 20, 
           title = "Deformations alongs PC axes") {
  if (is.factor(fac) & missing(col)) col <- 1:nlevels(fac)
  if (is.factor(fac) & missing(pch)) pch <- 1:nlevels(fac)
  nb.h <- ncol(Nef@coeff)/4
  pca(Nef, fac = fac, PCa = PCa, PCb = PCb, lab = FALSE, 
      zoom.x = zoom, ell = ell, pch = pch, col = col, title = title)
  def.par <- par(no.readonly = TRUE)
  pca <- prcomp(Nef@coeff)
  extrem <- apply(pca$x[, c(PCa, PCb)], 2, range)
  l2m <- function(m) matrix(c(m$x, m$y), nc = 2)
  W <- l2m(pca2shp(extrem[1, 1], 0, nb.h = nb.h, nb.pts = 500, 
        data = Nef@coeff, pl = 0))
    E <- l2m(pca2shp(extrem[2, 1], 0, nb.h = nb.h, nb.pts = 500, 
        data = Nef@coeff, pl = 0))
    S <- l2m(pca2shp(0, extrem[1, 2], nb.h = nb.h, nb.pts = 500, 
        data = Nef@coeff, pl = 0))
    N <- l2m(pca2shp(0, extrem[2, 2], nb.h = nb.h, nb.pts = 500, 
        data = Nef@coeff, pl = 0))
    cent <- l2m(pca2shp(0, 0, nb.h = nb.h, nb.pts = 500,
        data = Nef@coeff, pl = 0))
    plt <- par("plt")
    plt.x <- diff(plt[1:2])
    plt.y <- diff(plt[3:4])
    ins <- 0.02 #former inset adjustment                     
    W.x1 <- plt[1] + plt.x * ins
    W.x2 <- plt[1] + plt.x * (0.25 - ins)
    W.y1 <- mean(plt[3:4]) - plt.y * 0.25
    W.y2 <- mean(plt[3:4]) + plt.y * 0.25
    par(plt = c(W.x1, W.x2, W.y1, W.y2), new = TRUE)
    tps(cent, W, n = ncells, plot = 1)
    E.x1 <- plt[2] - plt.x * (0.25 + ins)
    E.x2 <- plt[2] - plt.x * ins
    par(plt = c(E.x1, E.x2, W.y1, W.y2), new = TRUE)
    tps(cent, E, n = ncells, plot = 1)
    S.x1 <- mean(plt[1:2]) - plt.x * 0.25
    S.x2 <- mean(plt[1:2]) + plt.x * 0.25
    S.y1 <- plt[3] + plt.x * ins
    S.y2 <- plt[3] + plt.y * (0.25 + ins)
    par(plt = c(S.x1, S.x2, S.y1, S.y2), new = TRUE)
    tps(cent, S, n = ncells, plot = 1)
    N.y1 <- plt[4] - plt.y * (0.25 + ins)
    N.y2 <- plt[4] - plt.x * ins
    par(plt = c(S.x1, S.x2, N.y1, N.y2), new = TRUE)
    tps(cent, N, n = ncells, plot = 1)
    par(def.par)})

# Plots the morphological space
setMethod(f="morph.sp", signature="Nef", definition=
  function(Nef,
           PCa = 1,
           PCb = 2,
           nb.PCa = 5,
           nb.PCb = 6,
		   fac = NA,
           morph.sp.extend = 1,
           zoom.extend = 1.2,
		   asp,
           pch = 20,
           shp.col = NA,
           shp.lwd = 1,
           shp.size,
           col = "grey40", 
           ell = FALSE,
           r = 1,
           lwd = 1,
           title = "Morphological space") {
    nb.h <- ncol(Nef@coeff)/4
    proj <- prcomp(Nef@coeff)$x[, c(PCa, PCb)]
    wdw <- apply(proj, 2, range) * zoom.extend * morph.sp.extend
    if (!missing(asp)) {
	plot(NA, xlim = wdw[, 1], ylim = wdw[, 2], asp = asp, las=1,
         xlab = paste("Principal Component", PCa),
         ylab = paste("Principal Component", PCb),
         main = title)
	} else {
	plot(NA, xlim = wdw[, 1], ylim = wdw[, 2], las=1,
         xlab = paste("Principal Component", PCa),
         ylab = paste("Principal Component", PCb),
         main = title)}
    abline(h = 0, v = 0, col = "grey90") ;  box()
    # Start factor
    if (is.factor(fac)) {
      if (missing(pch)) {pch <- (1:nlevels(fac))}
      if (missing(col)) {col <- (1:nlevels(fac))}
      k <- levels(fac)
      for (i in seq(along = k)) {
        s <- fac == k[i]
        points(proj[s, 1], proj[s, 2], pch = pch[i], col = col[i])
        if (ell) {panel.lm(proj[s, 1], proj[s, 2], col = col[i],
                           r = r, lwd = lwd)}}
    } else {points(proj[, 1], proj[, 2], pch = pch, col = col)}
    # End factor
    w <- apply(proj, 2, range) * morph.sp.extend
    if (missing(shp.size)) {
      cont <- pca2shp(0, 0, data = Nef@coeff, plot = 0)
      range.cont <- lapply(cont, function(x) diff(range(x)))
      shp.size <- ifelse(range.cont$x >= range.cont$y,
                         diff(w[1:2]/(12 * range.cont$x)),
                         diff(w[3:4]/(8 * range.cont$y)))}
    sp.pc1 <- seq(w[1, 1], w[2, 1], len = nb.PCa)
    sp.pc2 <- seq(w[1, 2], w[2, 2], len = nb.PCb)
    for (i in seq(along = sp.pc1)) {
      for (j in seq(along = sp.pc2)) {
        cont <- pca2shp(pc1 = sp.pc1[i], pc2 = sp.pc2[j],
                        nb.h = nb.h, nb.pts = 300,
                        data = Nef@coeff, plot = FALSE)
        shp <- list(x = sp.pc1[i] + cont$x * shp.size,
                    y = sp.pc2[j] + cont$y * shp.size)
        shp <- closed.outline(shp)
        if (!missing(shp.col)) {polygon(shp, col = shp.col, border = NA)}
        lines(shp, lwd = shp.lwd)}}})

# Plots the morphological space along PC axis
setMethod(f="morph.PC", signature="Nef", definition=
	function(Nef,
			 sd.nb=1,
			 pca.ax=seq(1, 3)){
  layout(matrix(1:(3*length(pca.ax)), nc=3, byrow=T))
  op <- par(no.readonly = TRUE)
  par(oma=rep(1,4), mar=rep(1,4))
  x.mu    <- apply(Nef@coeff, 2, mean)
  pca     <- prcomp(Nef@coeff)
  nb.h <- ncol(Nef@coeff)/4
  for (i in seq(along=pca.ax)) {
    x.sd    <- sd(pca$x[, pca.ax[i]])
    coe.min <- x.mu - sd.nb * x.sd*pca$rotation[, pca.ax[i]] 
    coe.max <- x.mu + sd.nb * x.sd*pca$rotation[, pca.ax[i]]
    coe.comp.plot <- function(coe, nb.h, title="") {
      st      <- 1 + nb.h * 0:3
      end     <- nb.h + nb.h * 0:3
      coe.shp <- iefourier(coe[st[1]:end[1]], coe[st[2]:end[2]],
                           coe[st[3]:end[3]], coe[st[4]:end[4]],
                           nb.h, 300)
      plot(NA, xlim=range(coe.shp$x)*1.2, ylim=range(coe.shp$y)*1.2, asp=1, main=title, frame=F, axes=F)
      polygon(closed.outline(coe.shp), col="grey90")}
    coe.comp.plot(coe.min, nb.h, title=paste("PC", pca.ax[i], " - ", sd.nb, "SD", sep=""))
    coe.comp.plot(x.mu,    nb.h, title="Mean shape")
    coe.comp.plot(coe.max, nb.h, title=paste("PC", pca.ax[i], " + ", sd.nb, "SD", sep=""))}
  par(op)})
  
# Calculates shape intermediates
setMethod(f="traj", signature="Nef", definition=
  function(Nef,
           fr = c(0, 0),
           to = c(1, 1),
           nb.int = 50,
		   nb.pts = 500,
           save = FALSE,
           prog = TRUE,
		   pause=TRUE){
  X11()
  pca(Nef)
  if (missing(fr)) {
    cat("Click on the plot to select the starting point\n")
    loc <- locator(1)
    fr <- c(loc$x, loc$y)}
  if (missing(to)) {
    cat("Click on the plot to select the ending point\n")
    loc <- locator(1)
    to <- c(loc$x, loc$y)}
  if (save) {
    dir <- paste(format(Sys.time(), "%m%d%y-%H%M"), "-traj", sep = "")
    dir.create(dir, showWarnings = FALSE)}
  x <- seq(fr[1], to[1], len = nb.int)
  y <- seq(fr[2], to[2], len = nb.int)
  par(mar=rep(0,4))
  w <- c(-1.3, 1.3)
  for (i in 1:nb.int) {
    cont <- pca2shp(pc1 = x[i], pc2 = y[i], data = Nef@coeff, nb.pts= nb.pts)
    plot(closed.outline(cont), xlim = w, ylim = w, asp = 1, type = "l", 
          xaxs = "i", yaxs = "i", ax = FALSE, ann = FALSE, lwd = 2)
    if (prog) {
      int <- (2.6/nb.int)
      rect(-1.3, -1.2, -1.3 + int * i, -1.15, border = NA, col = "red")
      rect(-1.3, -1.2, 1.3, -1.15, border = "black", col = NA)}
    if (pause) readline(prompt = "Press <Enter>") else Sys.sleep(0.01)
    if (save){ 
      savePlot(filename = paste(dir, "/", i, ".jpg", sep = ""), type = c("tiff"))}}
  readline(prompt = "Press <Enter> to continue...")
  dev.off()})

# Calculates shape intermediates
setMethod(f="tps.grid", signature="Nef", definition=
function (Nef,
          fr,
          to,
          nb.pts = 50,
          amp = 1,
          grid.size = 50,
          grid.col = "grey40",
          cont = TRUE,
          cont.col = c("dodgerblue3", "firebrick3"),
          cont.lwd = rep(3,2)) {
  # X11()
  def.par <- par(no.readonly = TRUE)
  if (missing(fr)) {
    morph.sp(Nef)
    cat("Click on the plot to select the starting point\n")
    loc <- locator(1)
    fr <- c(loc$x, loc$y)}
  if (missing(to)) {
    cat("Click on the plot to select the ending point\n")
    loc <- locator(1)
    to <- c(loc$x, loc$y)}
  nb.h <- ncol(Nef@coeff)/4
  shp.fr <- pca2shp(fr[1], fr[2], Nef@coeff, nb.h = nb.h, nb.pts = nb.pts, 
                    amp = amp, plot = FALSE)
  shp.to <- pca2shp(to[1], to[2], Nef@coeff, nb.h = nb.h, nb.pts = nb.pts, 
                    amp = amp, plot = FALSE)
  fr.m <- matrix(c(shp.fr$x, shp.fr$y), nc = 2)
  to.m <- matrix(c(shp.to$x, shp.to$y), nc = 2)
  tps(fr.m, to.m, grid.size, col = grid.col)
  if (cont) {
    fr.l <- list(x = fr.m[, 1], y = fr.m[, 2])
    to.l <- list(x = to.m[, 1], y = to.m[, 2])
    lines(closed.outline(fr.l), lwd = cont.lwd[1], col = cont.col[1])
    lines(closed.outline(to.l), lwd = cont.lwd[2], col = cont.col[2])}
  # dev.off()
  par(def.par)})

# Calculates shape intermediates
setMethod(f="tps.iso", signature="Nef", definition=
  function (Nef,
            fr,
            to,
            nb.pts = 200,
            amp = 1,
            iso.pts = 1000, 
            col.pal = topo.colors,
            col.lev = 500,
            cont.to = TRUE,
            cont.fr = TRUE, 
            cont.lev = 10,
            cont.col = c("dodgerblue3", "firebrick3"), 
            cont.lwd = rep(3,2)){
  # X11()
  def.par <- par(no.readonly = TRUE)
  if (missing(fr)) {
    morph.sp(Nef)
    cat("Click on the plot to select the starting point\n")
    loc <- locator(1)
    fr <- c(loc$x, loc$y)}
  if (missing(to)) {
    cat("Click on the plot to select the ending point\n")
    loc <- locator(1)
    to <- c(loc$x, loc$y)}
  nb.h <- ncol(Nef@coeff)/4
  shp.fr <- pca2shp(fr[1], fr[2], Nef@coeff, nb.h = nb.h, nb.pts = nb.pts, 
                    amp = amp, plot = FALSE)
  shp.to <- pca2shp(to[1], to[2], Nef@coeff, nb.h = nb.h, nb.pts = nb.pts, 
                  amp = amp, plot = FALSE)
  m.fr <- matrix(c(shp.fr$x, shp.fr$y), nc = 2)
  m.to <- matrix(c(shp.to$x, shp.to$y), nc = 2)
  sFE  <- spsample(Polygon(rbind(m.fr, m.fr[1, ])), iso.pts, 
      type = "regular")
  sR  <- sFE@coords
  sT  <- tps2d(sR, m.fr, m.to)
  def <- sqrt(apply((sT - sR)^2, 1, sum))
  x1  <- length(unique(sR[, 1]))
  y1  <- length(unique(sR[, 2]))
  im  <- matrix(NA, x1, y1)
  xind <- (1:x1)[as.factor(rank(sR[, 1]))]
  yind <- (1:y1)[as.factor(rank(sR[, 2]))]
  n <- length(xind)
  for (i in 1:n) {im[xind[i], yind[i]] <- def[i]}
  colors <- col.pal(col.lev)
  x <- sort(unique(sR[, 1]))
  y <- sort(unique(sR[, 2]))
  image(x, y, im, col = colors, asp = 1,
        axes = FALSE, frame = FALSE, ann = FALSE)
  contour(x, y, im, nlevels = cont.lev, add = TRUE)
  if (cont.fr) {
      lines(closed.outline(shp.fr), col = cont.col[1], lwd = cont.lwd[1])}
  if (cont.to) 
      {lines(closed.outline(shp.to), col = cont.col[2], lwd = cont.lwd[2])}
  # dev.off()
  par(def.par)})

# Calculates shape intermediates
setMethod(f="tps.vf", signature="Nef", definition=
  function (Nef,
            fr,
            to,
            nb.pts = 100,
            amp = 1,
            arr.nb = 300,
            arr.len = 0.05,
            arr.ang = 30,
            arr.col = "grey40",
            arr.pal = FALSE,
            arr.palette = topo.colors,
            arr.pal.lev = 20,
            arr.lwd = 1, 
            cont.col = c("dodgerblue3", "firebrick3"),
            cont.lwd = rep(3, 2)){
  e.d <- function (pt1, pt2) {sqrt(sum((pt1 - pt2)^2))}
  # X11()
  def.par <- par(no.readonly = TRUE)
  if (missing(fr)) {
    morph.sp(Nef)
    cat("Click on the plot to select the starting point\n")
    loc <- locator(1)
    fr <- c(loc$x, loc$y)}
  if (missing(to)) {
    cat("Click on the plot to select the ending point\n")
    loc <- locator(1)
    to <- c(loc$x, loc$y)}
  nb.h <- ncol(Nef@coeff)/4
  shp.fr <- pca2shp(fr[1], fr[2], Nef@coeff, nb.h = nb.h, nb.pts = nb.pts, 
                    amp = amp, plot = 0)
  shp.to <- pca2shp(to[1], to[2], Nef@coeff, nb.h = nb.h, nb.pts = nb.pts, 
                    amp = amp, plot = 0)
  m.fr <- matrix(c(shp.fr$x, shp.fr$y), nc = 2)
  m.to <- matrix(c(shp.to$x, shp.to$y), nc = 2)
  sFE <- spsample(Polygon(rbind(m.fr, m.fr[1, ])), arr.nb, type = "regular")
  sR <- sFE@coords
  sT <- tps2d(sR, m.fr, m.to)
  w <- c(range(c(shp.fr$x, shp.to$x)), range(c(shp.fr$y, shp.to$y))) * 1.1
  plot(NA, xlim = w[1:2], ylim = w[3:4], asp = 1,
       axes = FALSE, ann = FALSE, mar = rep(0, 4))
  if (arr.pal) {
    q.lev <- numeric(nrow(sR))
    for (i in seq(along = q.lev)) q.lev[i] <- e.d(sR[i, ], sT[i, ])
      q.lev <- cut(q.lev, breaks = arr.pal.lev, labels = FALSE)
      arr.col <- arr.palette(arr.pal.lev)[q.lev]
    } else {arr.col <- rep(arr.col, nrow(sR))}
  for (i in 1:dim(sR)[1]) {
    x0 <- sR[i, 1]
    y0 <- sR[i, 2]
    y1 <- sT[i, 2]
    x1 <- sT[i, 1]
    arrows(x0, y0, x1, y1, length = arr.len, angle = arr.ang, 
           lwd = arr.lwd, col = arr.col[i])}
  lines(rbind(m.fr, m.fr[1, ]), lwd = cont.lwd[1], col = cont.col[1])
  lines(rbind(m.to, m.to[1, ]), lwd = cont.lwd[2], col = cont.col[2])
  # dev.off()
  par(def.par)})