############################# Coo methods #######################################
# Coo class and builder
setClass("Coo", representation(coo="list"))
Coo <- function(list) {new(Class="Coo", coo=list)}

# show(Coo)
setMethod(f="show", signature="Coo", definition=function(object){
  out.len <- numeric()
  for (i in 1:length(object@coo)) {
    out.len <- append(out.len, length(object@coo[[i]]$x))} 
cat(
rep("*", 23, sep=""), "*** A Coo class object ", rep("*", 23, sep=""), "\n",
" ->  It contains ", length(out.len), " outlines with on average ",
  round(mean(out.len)), "+/-", round(sd(out.len)), " coordinates.\n",
" ->  Calibration methods can be passed to it.\n",
" ->  The list of (x; y) coordinates is stored in the @coo slot.\n",
" ->  See ?Coo for further details.\n",
rep("*", 67, sep=""), sep="")})  

# plot one or several outlines from a Coo object
setMethod(f="plot", signature="Coo", definition=
  function(x,y,range=1:length(x@coo), prompt=TRUE, col="grey40",...) {
    for (i in seq(along=range)){
        plot(x=x@coo[[range[i]]]$x, y=x@coo[[range[i]]]$y, col="red", type="l", lwd=3, an=FALSE, asp=1, axes=FALSE, frame=FALSE)
		title(main=names(x@coo)[range[i]], sub=paste(length(x@coo[[range[i]]]$x), "coordinates"))
        polygon(x@coo[[range[i]]], col=col)
        if(prompt) {readline(prompt = "Press <Enter> to continue...")}}})

# Calculates eFa if provided with a contour list
setMethod(f="get.Nef", signature="Coo", definition=
  function(Coo,
		   nb.h=32,
		   smooth.it=0,
		   fromrt = FALSE) {
    col.n <- paste(rep(LETTERS[1:4], each = nb.h), rep(1:nb.h, times = 4), sep="")
    coeff <- matrix(NA, nc = 4 * nb.h, nr = length(Coo@coo),
                    dimnames = list(names(Coo@coo), col.n))
    for (i in seq(along = Coo@coo)) {
        nef <- eFa(Coo@coo[[i]], nb.h = nb.h, smooth.it = smooth.it, fromrt = fromrt)
        if (nef$A[1] < 0) {
            nef$A <- (-nef$A)
            nef$B <- (-nef$B)
            nef$C <- (-nef$C)
            nef$D <- (-nef$D)
            nef$lnef <- (-nef$lnef)}
        coeff[i, ] <- c(nef$A, nef$B, nef$C, nef$D)}
    return(Nef(coeff))})

# Plots qualitative deviations
setMethod(f="dev.qual", signature="Coo", definition=
  function(Coo,
           id        = 1:length(Coo@coo),
           nb.h      = 32,
           smooth.it = 0,
           range     = seq(1, nb.h, len = 4)){
  def.par <- par(no.readonly = TRUE)
  range   <- round(range)
  n       <- length(range)
  cols    <-colorRampPalette(c("dodgerblue4", "firebrick1"))(n)
  nr <- ceiling(sqrt(n))
  for (k in 1:length(id)) {
  par(oma=c(0,0,3,0), mar=rep(1,4))
  layout(matrix(1:(nr^2), nr, nr, byrow = TRUE))
  ef <- efourier(Coo@coo[[id[k]]], nb.h = nb.h, smooth.it = smooth.it)
  ef.max <- iefourier(ef$an, ef$bn, ef$cn, ef$dn, nb.h, 500, ef$ao, ef$co)
  for (i in seq(along = range)) {
      ie  <- iefourier(ef$an, ef$bn, ef$cn, ef$dn, range[i], 500, ef$ao, ef$co)
      ief <- closed.outline(ie)
      plot(ef.max, type = "l", asp = 1, axes = FALSE, ann = FALSE)
      title(paste(range[i], "harmonics"))
      polygon(ef.max, col = "gold")
      lines(ief, col = cols[i], lwd = 2)}
  par(def.par) # if the user escapes, par is safe
  mtext(names(Coo@coo)[id[k]], side=3, line=1, cex=1.3, outer=TRUE)
  if(length(id)>1){readline(prompt = "Press <Enter> to continue...")}}
  par(def.par)})

# Plots quantitative deviations
setMethod(f="dev.quant", signature="Coo", definition=
  function(Coo,
           id=1,
           nb.h = 32,
           smooth.it = 0,
           plot=TRUE){
  l2m <- function(l1) {return(matrix(c(l1$x, l1$y), nc=2))}
  e.dm <- function(m1, m2) {
    if (class(m1)=="list")  m1 <- l2m(m1)
    if (class(m2)=="list")  m2 <- l2m(m2)
    return(apply((m1-m2)^2,1,function(x) sqrt(sum(x))))}
  align2m  <- function(m1, m2) {
    if (class(m1)=="list")  m1 <- l2m(m1)
    if (class(m2)=="list")  m2 <- l2m(m2)
    n   <- nrow(m1)
    d.M <- matrix(c(0:(n-1),rep(NA,n)), n, 2) # matrix to save alig dist
    for (i in 0:(n-1)) { # compute the distance for each alignement combinations
      a    <- (1+i):(n+i)
      d.M[i+1, 2] <- sum(e.dm(m1, m2[c(a[a<=n], a[a>n]%%n),]))}
    x      <- d.M[which(d.M[,2]==min(d.M[,2])),1] # minimal alignment distance
    al     <- (1+x):(n+x)
    return(al=m2[c(al[al<=n], al[al>n]%%n),])}
  true.dist <- function(ref, nb.h=32) {
    calliper.length <- max(dist(l2m(ref)))
    ref <- cont.sample(coo=ref, n=nb.h*2)
    ef  <- efourier(ref, nb.h = nb.h, smooth.it = 0)
    rec <- iefourier(ef$an, ef$bn, ef$cn, ef$dn, nb.h, nb.h*2, ef$ao, ef$co)
    rec.aligned <- align2m(ref, rec)
    dev   <- e.dm(ref, rec.aligned)/calliper.length
    return(dev)}
  def.par <- par(no.readonly = TRUE)
  if (missing(id)) {id <- (1:length(Coo@coo))}
  range   <- 4:nb.h
  res <- matrix(NA, nr=length(range), nc=6,
                dimnames=list(range, c("5%", "25%", "50%","75%", "95%", "100%")))
  for (j in seq(along=range)) {
    cat("Computing deviations for", range[j], "over", nb.h, "harmonics\n")
    dev <- numeric()
    for (i in 1:length(id)) {
      dev <- c(dev, true.dist(Coo@coo[id[i]][[1]], nb.h=range[j]))}
      res[j,] <- as.numeric(quantile(dev, probs=c(0.05,0.25,0.5,0.75,0.95,1)))}
  if (plot){
    int.x   <- c(range,rev(range))
    mid.int <- c(res[,2],rev(res[,4]))
    ext.int <- c(res[,1],rev(res[,5]))
    plot(NA, xlim = range(range), ylim = c(0, max(res)), las=1, yaxs="i", 
         xlab = "Number of harmonics included", ylab = "Deviation from original outline",
         main="Deviations observed", sub=paste("(",length(id), "outlines )"))
    polygon(int.x, ext.int, col="grey80", border=NA)
    polygon(int.x, mid.int, col="grey60", border=NA)
    lines(range,  res[,3],  lty=1) ; points(range, res[,3], pch=1)
    lines(range, res[,6], lty=2) ; points(range, res[,6], pch=4) 
    w <- par("usr")
    legend(w[2] * 0.618, w[4]*0.9, lty = c(1,2,-1,-1),
           fill=c(NA,NA,"grey60","grey80"), border=rep(NA,4), 
          legend = c("Average deviation", "Maximal deviation",
                     "0.25-0.75 quantiles", "0.05-0.95 quantiles"), bty = "n")}
  par(def.par)
  return(res)})

# Plots Fourier power spectrum
setMethod(f="harm.pow", signature="Coo", definition=
  function(Coo,
          id,
          nb.h = 32,
          smooth.it = 0,
          plot = TRUE,
          max.h = nb.h,
          first = TRUE) {
    if(missing(id)) {range <- (1:length(Coo@coo))} else {range <- id}
    idx <- ifelse(first, 1, 2)
    res <- matrix(NA, length(range), nb.h - idx + 1,
                  dimnames = list((names(Coo@coo)[range]), idx:nb.h))
    for (i in seq(along=range)) {
        ef <- efourier(Coo@coo[[range[i]]], nb.h = nb.h, smooth.it = smooth.it)
        pow <- ((ef$an^2 + ef$bn^2 + ef$cn^2 + ef$dn^2)/2)[idx:nb.h]
        res[i, ] <-  (cumsum(pow)/sum(pow))}
    res <- apply(res, 2,
                 function(x) quantile(x,probs=c(0.05,0.25,0.5,0.75,0.95,1), na.rm=FALSE))
  if (plot){
    int.x   <- c((idx : nb.h), (nb.h : idx))
    mid.int <- c(res[2,],rev(res[4,]))
    ext.int <- c(res[1,],rev(res[5,]))
    plot(NA, xlim = c(0, nb.h), ylim = c(min(res), 1), las=1, yaxs="i", xaxs="i", 
         xlab = "Number of harmonics included", ylab = "Cumulative harmonic power",
         main="Fourier power spectrum", sub=paste("(",length(range), "outlines )"))
    polygon(int.x, ext.int, col="grey80", border=NA)
    polygon(int.x, mid.int, col="grey60", border=NA)
    lines((idx:nb.h), res[3,], lty=1) ; points((idx:nb.h), res[3,], pch=1)
    lines((idx:nb.h), res[6,], lty=2) ; points((idx:nb.h), res[6,], pch=4) 
    w <- par("usr")
    legend(w[2] * 0.618, w[4]*0.9, lty = c(1,2,-1,-1),
           fill=c(NA,NA,"grey60","grey80"), border=rep(NA,4), 
          legend = c("Average deviation", "Maximal deviation",
                     "0.25-0.75 quantiles", "0.05-0.95 quantiles"), bty = "n")}
    return(res)})       