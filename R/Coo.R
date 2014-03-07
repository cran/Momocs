#--------------------------------------------------------------------#
# Coo methods                                                        #
#--------------------------------------------------------------------#

# Coo declaration ####################################################
setClass(
  Class = "Coo",
  representation(coo        = "list",
                 names      = "character",
                 fac        = "data.frame",
                 ldk        = "list",
                 coo.nb     = "numeric",
                 coo.len    = "numeric",
                 coo.closed = "logical",
                 details    = "list"),
  validity=function(object){
    if(any(lapply(object@coo, class) != "matrix")) cat(" * A list of coordinates as 2-col matrices must be provided.\n")
    if(any(lapply(object@coo, ncol)  != 2))        cat(" * A list of coordinates as 2-col matrices must be provided.\n")
  }
)

# Coo builder ########################################################
# Builder
Coo <- function(coo, ...){
 if (class(coo)=="Coo") {
   return(new(Class="Coo",
       coo=coo@coo,
       names=coo@names,
       coo.nb=length(coo@coo),
       coo.len=as.numeric(lapply(coo@coo, nrow)),
       coo.closed=unlist(lapply(coo@coo, is.closed)),        
       fac=as.data.frame(apply(coo@fac, 2, as.factor)),
       ldk=coo@ldk, ...))
   } else {
  if (is.matrix(coo)) coo <- list(coo) #if only one shape is provided
  if(is.array(coo))   coo <- a2l(coo)
	if(is.null(names(coo))) {
    names(coo) <- paste("shp", 1:length(coo), sep="") }  
	new(Class="Coo",
		coo=coo,
        names=names(coo),
        coo.nb=length(coo),
        coo.len=as.numeric(lapply(coo, nrow)),
        coo.closed=!unlist(lapply(coo, is.closed)), 
      ...)}}

# Coo getters ########################################################
setMethod(f = "[", signature = "Coo", definition=function(x, i, j, drop){
  if (is.integer(i)) {
    return(x@coo[i])}
  if (is.numeric(i)) {
    return(x@coo[[i]])}
  if (missing(i)) {
    if (ncol(x@fac)!=0) {
      fac <- as.data.frame(factor(x@fac[j,]))
      names(fac) <- names(x@fac)
    } else {fac <- data.frame()}
    return(Coo(x@coo[j], fac=fac, ldk=x@ldk))}
  if (i=="coo")     return(x@coo)
  if (i=="names")   return(x@names[j])
  if (i=="fac")     return(x@fac[j,])
  if (i=="ldk")     return(x@ldk[j])
  if (i=="coo.nb")  return(x@coo.nb)
  if (i=="coo.len") return(x@coo.len[j])
  if (i=="details") return(x@details)
  else stop("Wrong slot name or subset specification.")})

# Coo setters ########################################################
setReplaceMethod(f = "[", signature = "Coo", definition =
  function(x, i, j, value){
    switch(EXPR=i,
           "coo"    = {x@coo[[j]] <- value},
           "names"  = {x@names[j] <- value},
           "fac"    = {x@fac      <- as.data.frame(value)},
           "ldk"    = {x@ldk   <- value },
           "details"= {x@details <- as.character(value)},
           stop(" * This attribute does not exist or cannot be modified."))
    validObject(x)
    return(x)})

# show(Coo) ##########################################################
setMethod(f="show", signature="Coo", definition=function(object){
  cat("A Coo object (see ?Coo) \n")
  cat(rep("*", 30),"\n", sep="")
  
  cat("\nGeneral\n")
  cat(rep("-", 20),"\n", sep="")
  cat(" -", object@coo.nb, ifelse(object@coo.nb<2, "outline\n", "outlines\n"))
  cat(" -", round(mean(object@coo.len)), "+/-", round(sd(object@coo.len)), "coordinates per outline\n")
  if (length(object@ldk)!=0) cat(" -", length(object@ldk[[1]]), "landmark(s) defined\n") else cat(" - No landmark defined\n")
  
  if (all(object@coo.closed)) {
    cat(" - All outlines are closed\n")
    } else {
      if (any(!object@coo.closed)) {
        cat(" - All outlines are unclosed\n")
      } else cat(" -", sum(object@coo.closed), "outlines are closed\n")}
  
  if (ncol(object@fac)!=0) cat(" -", ncol(object@fac), "grouping factor(s) defined\n") else cat(" - No groups defined\n")
  
  eg <- sample(object@coo.nb, 1)
  cat("\nCoordinates: @coo (e.g. '", object@names[eg], "')\n", sep="")
  cat(rep("-", 20),"\n", sep="")
  k <- ifelse (nrow(object@coo[[eg]]) < 4, nrow(object@coo[[eg]]), 4)
  print(object@coo[[eg]][1:k,])
  
  cat("\nNames: @names\n")
  cat(rep("-", 20),"\n", sep="")
  if(object@coo.nb>20) {
    cat(object@names[1:20], fill=120)
    cat("\n...") 
  } else {
    cat(object@names, fill=120)}
  cat("\n")
  
  if (ncol(object@fac)){
    cat("\nFactors: @fac\n")
    cat(rep("-", 20),"\n", sep="")
    for (i in 1:ncol(object@fac)){
      cat(names(object@fac[i]), ":", levels(object@fac[, i]), "\n")
    }}})

# plot(Coo) and other plotting methods ###############################
setMethod(f="plot", signature="Coo", definition=
  function(x, y, cols, borders, ldk=TRUE, ldk.pch=3, ldk.col="red", ldk.cex=1, ...){
  if (missing(cols)) {
    cols     <- rep(NA,    x@coo.nb)
  } else {
    if (length(cols)!=x@coo.nb){
      cols <- rep(cols[1], x@coo.nb)}}
  
  if (missing(borders)) {
    borders  <- rep("#708090", x@coo.nb)
  } else {
    if (length(borders)!=x@coo.nb){
      borders <- rep(borders[1], x@coo.nb)}}
  wdw <- apply(l2a(lapply(x@coo, function(x) apply(x, 2, range))), 2, range)
  coo.plot(xlim=wdw[, 1], ylim=wdw[, 2])
  for (i in 1:x@coo.nb) {
    coo.draw(x@coo[[i]], col=cols[i], border=borders[i], points=FALSE, ...)
    if (ldk & length(x@ldk)!=0) points(x@coo[[i]][x@ldk[[i]],], pch=ldk.pch, col=ldk.col, cex=ldk.cex)}})


setGeneric(name= "stack",   def=function(x, cols, borders,
                                         ldk=TRUE, ldk.pch=3, ldk.col="red", ldk.cex=1, ...){standardGeneric("stack")})

setMethod(f="stack", signature="Coo", definition=
  function(x, cols, borders,
           ldk=TRUE, ldk.pch=3, ldk.col="red", ldk.cex=1, ...){
  if (missing(cols)) {
    cols     <- rep(NA,    x@coo.nb)
  } else {
    if (length(cols)!=x@coo.nb){
      cols <- rep(cols[1], x@coo.nb)}}
  
  if (missing(borders)) {
    borders  <- rep("#708090", x@coo.nb)
  } else {
    if (length(borders)!=x@coo.nb){
      borders <- rep(borders[1], x@coo.nb)}}
  
  wdw <- apply(l2a(lapply(x@coo, function(x) apply(x, 2, range))), 2, range)
  coo.plot(xlim=wdw[, 1], ylim=wdw[, 2])
  for (i in 1:x@coo.nb) {
    coo.draw(x@coo[[i]], col=cols[i], border=borders[i], points=FALSE, ...)
    if (ldk & length(x@ldk)!=0) points(x@coo[[i]][x@ldk[[i]],], pch=ldk.pch, col=ldk.col, cex=ldk.cex)}})

setGeneric(name= "diapo",   def=function(Coo, id=1, col=NA, border="#708090",
                                         ldk=TRUE, ldk.pch=3, ldk.col="red", ldk.cex=1, ...){standardGeneric("diapo")})
setMethod(f="diapo", signature="Coo", definition= function(Coo, id=1, col=NA, border="#708090",
                                                           ldk=TRUE, ldk.pch=3, ldk.col="red", ldk.cex=1, ...) {
  for (i in id:Coo@coo.nb) {
    coo.plot(Coo@coo[[i]], col=col, border=border, ...)
    title(main=Coo@names[i], outer=TRUE, line=-2)
    if (ldk & length(Coo@ldk)!=0) points(Coo@coo[[i]][Coo@ldk[[i]],], pch=ldk.pch, col=ldk.col, cex=ldk.cex)
    readline(prompt = "Press <Enter> to continue...")}})

setGeneric(name= "panel",   def=function(Coo, cols, borders, names=NULL, ...){standardGeneric("panel")})
setMethod(f="panel", signature="Coo", definition= function(Coo, cols, borders, names=NULL, cex.names=0.6, ...){
  if (missing(cols)) {
    cols     <- rep(NA,    Coo@coo.nb)
  } else {
    if (length(cols)!=Coo@coo.nb){
      cols <- rep(cols[1], Coo@coo.nb)}}
  
  if (missing(borders)) {
    borders  <- rep("#708090", Coo@coo.nb)
  } else {
    if (length(borders)!=Coo@coo.nb){
      borders <- rep(borders[1], Coo@coo.nb)}}
  
  pos <- coo.list.panel(Coo@coo, cols=cols, borders=borders, ...)
  if (!is.null(names)){
    if (is.logical(names)) {
      text(pos[,1], pos[,2], labels=Coo@names, cex=cex.names)
    } else {    
      if (length(names)!=Coo@coo.nb) stop("The length of names provided and the number of outlines should be the same")
      text(pos[,1], pos[,2], labels=names, cex=cex.names)
    }}
})

# Methods on coordinates #############################################
setGeneric(name= "Coo.align",   def=function(Coo, ...){standardGeneric("Coo.align")})
setMethod(f="Coo.align", signature="Coo", definition= function(Coo){
  Coo2 <- Coo
  Coo2@coo <- lapply(Coo@coo, coo.align)
  return(Coo2)})

setGeneric(name= "Coo.center",   def=function(Coo, ...){standardGeneric("Coo.center")})
setMethod(f="Coo.center", signature="Coo", definition=function(Coo){
  Coo2 <- Coo
  Coo2@coo <- lapply(Coo@coo, coo.center)
  return(Coo2)})

setGeneric(name= "Coo.smooth",   def=function(Coo, ...){standardGeneric("Coo.smooth")})
setMethod(f="Coo.smooth", signature="Coo", definition=function(Coo, n){
  Coo2 <- Coo
  Coo2@coo <- lapply(Coo@coo, coo.smooth, n)
  return(Coo2)})

setGeneric(name= "Coo.sample",   def=function(Coo, ...){standardGeneric("Coo.sample")})
setMethod(f="Coo.sample", signature="Coo", definition=function(Coo, nb.pts=100){
  Coo2 <- Coo
  Coo2@coo <- lapply(Coo@coo, coo.sample, nb.pts)
  return(Coo2)})

setGeneric(name= "Coo.template",   def=function(Coo, size=1){standardGeneric("Coo.template")})
setMethod(f="Coo.template", signature="Coo", definition=function(Coo, size=1){
  Coo2 <- Coo
  Coo2@coo <- lapply(Coo@coo, coo.template, size)
  return(Coo2)})

setGeneric(name= "Coo.slide",   def=function(Coo, ...){standardGeneric("Coo.slide")})
setMethod(f="Coo.slide", signature="Coo", definition=function(Coo, ldk.id=1){
  if (length(Coo@ldk)==0) stop(" * First define landmarks on the Coo object. See ?defLandmarks.")
  Coo2 <- Coo
  for (i in 1:Coo@coo.nb) {
    Coo2@coo[[i]] <- coo.slide(Coo@coo[[i]], Coo@ldk[[i]][ldk.id])
    Coo2@ldk[[i]] <- (Coo@ldk[[i]] - (Coo@ldk[[i]][ldk.id] -1)) %% Coo@coo.len[i]
  }
  return(Coo2)})

setGeneric(name= "Coo.close",   def=function(Coo, ...){standardGeneric("Coo.close")})
setMethod(f="Coo.close", signature="Coo", definition= function(Coo){
  if (all(Coo@coo.closed)) return(cat(" * All outlines were already closed.\n"))
  coo2close <- !Coo@coo.closed
  lapply(Coo@coo[coo2close], coo.close)
  Coo@coo.closed[coo2close] <- TRUE
  cat(sum(coo2close),"outlines closed\n")
  return(Coo)})

setGeneric(name= "Coo.unclose",   def=function(Coo, ...){standardGeneric("Coo.unclose")})
setMethod(f="Coo.unclose", signature="Coo", definition= function(Coo){
  if (!all(Coo@coo.closed)) return(cat(" * All outlines were already unclosed.\n"))
  coo2unclose <- Coo@coo.closed
  lapply(Coo@coo[coo2unclose], coo.unclose)
  Coo@coo.closed[coo2unclose] <- FALSE
  cat(sum(coo2unclose),"outlines unclosed\n")
  return(Coo)})  

# Landmarks and cie
setGeneric(name= "defLandmarks",   def=function(Coo, nb.ldk){standardGeneric("defLandmarks")})
setMethod(f="defLandmarks", signature="Coo", definition=function(Coo, nb.ldk){
  if (missing(nb.ldk)) stop(" * 'nb.ldk' must be specified.")
  ldk <- list()
  for (i in seq(along=Coo@coo)){
    Coo@ldk[[i]] <- coo.ldk(Coo@coo[[i]], nb.ldk=nb.ldk)
  }
  return(Coo)})

# get array of coordinates from landmarks identified on Coo object
setGeneric(name= "cooLandmarks",   def=function(Coo){standardGeneric("cooLandmarks")})
setMethod(f="cooLandmarks", signature="Coo", definition=function(Coo){
  if (length(Coo@ldk)==0) stop(" * No landmarks defined. See ?defLandmarks.")
  nb.ldk <- length(Coo@ldk[[1]])
  cooA <- array(NA, dim=c(nb.ldk, 2, Coo@coo.nb),
                dimnames=list(paste("ldk", 1:nb.ldk, sep=""),
                              c("x", "y"),
                              Coo@names))
  for (i in seq(along=Coo@coo)){
    cooA[,,i] <- Coo@coo[[i]][Coo@ldk[[i]],]
  }
  return(cooA)})

setGeneric(name= "baseline",   def=
  function(Coo, ldk1=1, ldk2=2, t1=c(-0.5, 0), t2=c(0.5, 0)){standardGeneric("baseline")})
setMethod(f="baseline", signature="Coo", definition=
  function(Coo, ldk1=1, ldk2=2, t1=c(-0.5, 0), t2=c(0.5, 0)){
    for (i in 1:Coo@coo.nb){
      Coo@coo[[i]] <- coo.baseline(Coo@coo[[i]],
                                   ldk1=Coo@ldk[[i]][ldk1], ldk2=Coo@ldk[[i]][ldk2],
                                   t1=t1, t2=t2)}
    return(Coo)})

setGeneric(name= "procGPAlign",   def=
  function(Coo, tol=1e-30){standardGeneric("procGPAlign")})
setMethod(f="procGPAlign", signature="Coo", definition=
  function(Coo, tol=1e-30){
  if (length(Coo@ldk)==0) stop(" * No landmarks defined. See ?defLandmanarks.")
  Coo2 <- Coo
  Coo2@coo <- lapply(Coo@coo, function(x) coo.scale(coo.center(x)))
  ref  <- cooLandmarks(Coo2)
  tar <- procGPA(ref, tol1=tol, proc.output=TRUE)$rotated
  # would benefit to be handled by coo.baseline
  for (i in 1:Coo2@coo.nb) {
    tari <- tar[, , i]
    refi <- ref[, , i]
    t1x <- tari[1, 1]
    t1y <- tari[1, 2]
    t2x <- tari[2, 1]
    t2y <- tari[2, 2]
    r1x <- refi[1, 1]
    r1y <- refi[1, 2]
    r2x <- refi[2, 1]
    r2y <- refi[2, 2]
    # translation
    t <- tari[1, ] - refi[1, ]
    refi <- coo.trans(refi, t[1], t[2])
    # rotation  
    tx <- t2x - t1x
    ty <- t2y - t1y
    rx <- r2x - r1x
    ry <- r2y - r1y
    vi <- vecs.param(rx, ry, tx, ty)
    coo.i <- Coo2@coo[[i]]
    coo.i <- coo.trans(coo.i, t[1]-t1x, t[2]-t1y)
    coo.i <- coo.i / vi$r.norms
    coo.i <- coo.rotate(coo.i, -vi$d.angle)
    coo.i <- coo.trans(coo.i, t1x, t1y)
    Coo2@coo[[i]] <- coo.i
  }
  return(Coo2)})

# Fourier computation ################################################
# eFourier computation 
setGeneric(name= "eFourier",
           def=function(Coo,
           nb.h      = 32,
           smooth.it = 0,
           norm      = TRUE,
           start     = FALSE){standardGeneric("eFourier")})

setMethod(f="eFourier", signature="Coo", definition=
  function(Coo,
           nb.h      = 32,
           smooth.it = 0,
           norm       = TRUE,
           start     = FALSE) {
    q <- floor(min(Coo@coo.len)/2) - 1
    if (missing(nb.h))  {
      nb.h <- if (q >= 32) { 32 } else { q }
      cat(paste(" * 'nb.h' not provided and set to", nb.h))}
    if(nb.h  > (q+1)*2) {
      nb.h <- q # should not be 1
      warning(" * At least one outline has ", (q+1)*2, " coordinates. The number of harmonics has been set to: ", q)}
    coo <- Coo@coo
    col.n <- paste0(rep(LETTERS[1:4], each = nb.h), rep(1:nb.h, times = 4))
    coe <- matrix(ncol = 4 * nb.h, nrow = length(coo), dimnames = list(names(coo), col.n))
    for (i in seq(along = coo)) {
      ef <- efourier(coo[[i]], nb.h = nb.h, smooth.it = smooth.it, silent = FALSE)
      if (norm) {
        ef <- efourier.norm(ef, start=start)
        if (ef$A[1] < 0) {
          ef$A <- (-ef$A)
          ef$B <- (-ef$B)
          ef$C <- (-ef$C)
          ef$D <- (-ef$D)
          ef$lnef <- (-ef$lnef)}
        coe[i, ] <- c(ef$A, ef$B, ef$C, ef$D)
      } else {
        coe[i, ] <- c(ef$an, ef$bn, ef$cn, ef$dn)}}
    return(Coe(coe, fac=Coo@fac, method="eFourier"))})

# rFourier computation
setGeneric(name= "rFourier",   def=function(Coo,
                                            nb.h      = 40,
                                            smooth.it = 0,
                                            norm = TRUE){standardGeneric("rFourier")})

setMethod(f="rFourier", signature="Coo", definition=
  function(Coo,
           nb.h      = 40,
           smooth.it = 0,
           norm = TRUE) {
    q <- floor(min(Coo@coo.len)/2) - 1
    if (missing(nb.h))  {
      nb.h <- if (q >= 32) { 32 } else { q }
      cat(paste("  * nb.h not provided and set to", nb.h))}
    if(nb.h  > (q+1)*2) {
      nb.h <- q # should not be 1
      warning("At least one outline has ", (q+1)*2, " coordinates. The number of harmonics has been set to: ", q)}
    coo <- Coo@coo
    col.n <- paste0(rep(LETTERS[1:2], each = nb.h), rep(1:nb.h, times = 2))
    coe <- matrix(ncol = 2 * nb.h, nrow = length(coo), dimnames = list(names(coo), col.n))
    for (i in seq(along = coo)) {
      rf <- rfourier(coo[[i]], nb.h = nb.h, smooth.it = smooth.it, norm=norm)
      coe[i, ] <- c(rf$an, rf$bn)}
    return(Coe(coe, fac=Coo@fac, method="rFourier"))})


# tFourier computation
setGeneric(name= "tFourier",   def=function(Coo,
                                            nb.h      = 40,
                                            smooth.it = 0,
                                            norm=TRUE){standardGeneric("tFourier")})

setMethod(f="tFourier", signature="Coo", definition=
  function(Coo,
           nb.h      = 40,
           smooth.it = 0,
           norm=TRUE) {
    q <- floor(min(Coo@coo.len)/2) - 1
    if (missing(nb.h))  {
      nb.h <- if (q >= 32) { 32 } else { q }
      cat(paste("  * nb.h not provided and set to", nb.h))}
    if(nb.h  > (q+1)*2) {
      nb.h <- q # should not be 1
      warning("At least one outline has ", (q+1)*2, " coordinates. The number of harmonics has been set to: ", q)}
    coo <- Coo@coo
    col.n <- paste0(rep(LETTERS[1:2], each = nb.h), rep(1:nb.h, times = 2))
    coe <- matrix(ncol = 2 * nb.h, nrow = length(coo), dimnames = list(names(coo), col.n))
    for (i in seq(along = coo)) {
      tf <- tfourier(coo[[i]], nb.h = nb.h, smooth.it = smooth.it, norm=norm)
      coe[i, ] <- c(tf$an, tf$bn)}
    return(Coe(coe, fac=Coo@fac, method="tFourier"))})

# hqual ##########################################################
setGeneric(name= "hqual",     def=function(
   Coo,
   method = c("efourier", "rfourier", "tfourier"), 
   id        = sample(Coo@coo.nb, 1),
   smooth.it = 0,
   harm.range= c(1, 2, 4, 8, 16, 32),
   scale = TRUE,
   center = TRUE,
   align = TRUE,
   plot.method = c("stack", "panel")[1],
   legend = TRUE,
   legend.title = "# harmonics",
   palette = col.summer,
   shp.col="#70809033",
   shp.border="#708090EE"){standardGeneric("hqual")})

setMethod(f="hqual", signature="Coo", definition=
  function(Coo,
           method = c("efourier", "rfourier", "tfourier"),
           id        = sample(Coo@coo.nb, 1),
           smooth.it = 0,
           harm.range= c(1, 2, 4, 8, 16, 32),
           scale = TRUE,
           center = TRUE,
           align = TRUE,
           plot.method = c("stack", "panel")[1],
           legend = TRUE,
           legend.title = "# harmonics",
           palette = col.summer,
           shp.col="#70809033",
           shp.border="#708090EE"){
    # for one single outline
    if (missing(method)) {
      cat(" * Method not provided. efourier is used.\n")
      method   <- efourier
      method.i <- efourier.i 
      } else {
      p <- pmatch(tolower(method), c("efourier", "rfourier", "tfourier"))
      if (is.na(p)) { warning(" * Unvalid method. efourier is used.\n")
        } else {
        method   <- switch(p, efourier,   rfourier,   tfourier)
        method.i <- switch(p, efourier.i, rfourier.i, tfourier.i)}}
    
    # check for too ambitious harm.range
    if (max(harm.range) > (min(Coo@coo.len)/2 + 1)) {
      harm.range <- floor(seq(1, q/2 - 1, length=6))
      cat(" * harm.range was too high and set to: ", harm.range, ".\n")}
    coo <- Coo@coo[[id]]
    if (scale)  coo <- coo.scale(coo)
    if (center) coo <- coo.center(coo)
    if (align)  coo <- coo.align(coo)
    res <- list()
    for (i in seq(along=harm.range)) {
      res[[i]] <- method.i(method(coo, nb.h=max(harm.range), smooth.it=smooth.it), nb.h=harm.range[i])}
    # plotting
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    cols <- paste(palette(length(harm.range)), "88", sep="")
    if (plot.method=="stack") {
      coo.plot(coo.smooth(coo, smooth.it), main=Coo@names[id], col=shp.col, border=shp.border)
      for (i in seq(along=harm.range)) {lines(res[[i]], col=cols[i], lwd=1)}
      if (legend) {
        legend("topright", legend = as.character(harm.range), bty="o",
               col = cols, lty = 1, lwd=1, bg="#FFFFFFCC", cex=0.7,
               title = legend.title)}
      if (center) {points(0, 0, pch=3, col=shp.col)}
    } else {
      if (plot.method=="panel") {
        #par(oma=c(1, 1, 3, 0))
        pos <- coo.list.panel(res, cols=cols)
        if (legend) {text(x=pos[, 1], y=pos[, 2], as.character(harm.range))}
        mtext(Coo@names[id], side=3, line=0, cex=1.3, outer=TRUE)
      }}})

# hquant #########################################################
setGeneric(name= "hquant", def=function(
  Coo,
  method = c("efourier", "rfourier", "tfourier"),
  id        = 1,
  smooth.it = 0,
  harm.range = seq(4, 20, 4),
  norm.centsize = TRUE,
  dist.method = edm.nearest,
  dist.nbpts = 120,
  plot = TRUE,
  dev.plot=TRUE,
  title = "Deviations along the outline",
  legend = TRUE,
  legend.title = "# harmonics",
  palette = col.summer,
  lineat.y=c(0.5, 0.1, 0.01)){standardGeneric("hquant")})

setMethod(f="hquant", signature="Coo", definition=
  function(Coo,
           method = c("efourier", "rfourier", "tfourier"),
           id        = 1,
           smooth.it = 0,
           harm.range= seq(4, 20, 4),
           norm.centsize = TRUE,
           dist.method = edm.nearest,
           dist.nbpts = 120,
           plot = TRUE,
           dev.plot=TRUE,
           title = "Deviations along the outline",
           legend = TRUE,
           legend.title = "# harmonics",
           palette = col.summer,
           lineat.y=c(0.5, 0.1, 0.01)){
    if (missing(method)) {
      cat("  * Method not provided. efourier is used.\n")
      method   <- efourier
      method.i <- efourier.i 
    } else {
      p <- pmatch(tolower(method), c("efourier", "rfourier", "tfourier"))
      if (is.na(p)) { warning("Unvalid method. efourier is used.")
      } else {
        method   <- switch(p, efourier,   rfourier,   tfourier)
        method.i <- switch(p, efourier.i, rfourier.i, tfourier.i)}}
    # We define the highest possible nb.h along Coo@coo[id]
    min.nb.pts <- min(unlist(lapply(Coo@coo[id], nrow)))
    nb.h.best  <- floor(min.nb.pts/2)-1
    # we handle too ambitious harm.range
    if (max(harm.range) > nb.h.best) {
        harm.range <- floor(seq(4, nb.h.best, length=6))
      cat("  * 'harm.range' was too high and set to: ", harm.range, ".\n")}
    
    # we prepare the results array
    nb.pts <- ifelse(dist.nbpts == "max", 2*nb.h.best, dist.nbpts)
    nr <- length(harm.range)
    nc <- nb.pts
    nk <- length(id)
    res <- array(NA, dim=c(nr, nc, nk),
                 dimnames=list(paste0("h", harm.range),
                               paste("pt", 1:nb.pts),
                               Coo@names[id]))
    # progressbar
    if (nk > 5) {
      pb <- txtProgressBar(1, nk)
      t <- TRUE } else {t <- FALSE}
    # the core loops that will calculate deviations
    for (ind in seq(along=id)) {
      coo <- Coo@coo[[id[ind]]]
      # below, the best possible fit
      coo.best <- l2m(method.i(method(coo, nb.h=nb.h.best, smooth.it=smooth.it), nb.pts=nb.pts))
      for (i in seq(along=harm.range)) {
        # for each number of harmonics we calculate deviation with the FUN=method
        coo.i <- l2m(method.i(method(coo, nb.h=harm.range[i], smooth.it=smooth.it, silent=TRUE), nb.pts=nb.pts))
        res[i, , ind] <- dist.method(coo.best, coo.i)
      }
      # we normalize by the centroid size
      if (norm.centsize) {res[,,ind] <- res[,,ind]/coo.centsize(coo)}
      if (t) setTxtProgressBar(pb, ind)}
    # below we manage for single/several individuals
    if (nk > 1) { # if more than 1, we calculate median and sd
      m <- apply(res, 1:2, median)
      d <- apply(res, 1:2, sd)
    } else {
      m <- res[,,1]
      d <- NULL}
    # plotting stuff
    if (plot) {
      cols <- palette(nr)
      if (nk > 1) {ylim <- c(0, max(m+d, na.rm=TRUE))} else {ylim <- range(m)}
      if (norm.centsize) {
        ylab = "Deviation (in % of the centroid size)"
      } else {
        ylab = "Deviation (in original units)"}
      plot(NA, xlim=c(1, nc), ylim=ylim,
           xlab="Points sampled along the outline",
           ylab=ylab, main=title,
           xaxs="i", yaxs="i", axes=FALSE)
      axis(1, at=seq(0, dist.nbpts, length=5))
      axis(2)
      abline(h=lineat.y, lty=2, col="grey90")
      # if you want deviations, here they are
      if (dev.plot) {
        if (nk > 1) {dev.plot(m, d, cols=cols) } else {
          for (i in 1:nr) {
            lines(1:ncol(m), m[i, ], col=cols[i])}}}
      # same for legend
      if (legend) {
        legend("topright", legend = as.character(harm.range), bty="o",
               col = cols, lty = 1, lwd=1, bg="#FFFFFFCC", inset=0.005, cex=0.7,
               title = legend.title)}
      box() }
    return(list(res=res, m=m, d=d))})

# harm.pow ###########################################################
setGeneric(name= "hpow", def=function(    
  Coo,
  method = c("efourier", "rfourier", "tfourier"),
  id=1:Coo@coo.nb,
  probs=c(0, 0.5, 1),
  nb.h = 24,
  drop   = 1,
  smooth.it = 0,
  plot = TRUE,
  legend = FALSE,
  title="Fourier power spectrum",
  lineat.y=c(0.9, 0.95, 0.99, 0.999),
  bw=0.1){standardGeneric("hpow")})
  
setMethod(f="hpow", signature="Coo", definition=
  function(Coo,
           method = c("efourier", "rfourier", "tfourier"),
           id=1:Coo@coo.nb,
           probs=c(0, 0.5, 1),
           nb.h = 24,
           drop   = 1,
           smooth.it = 0,
           plot = TRUE,
           legend = FALSE,
           title="Fourier power spectrum",
           lineat.y=c(0.9, 0.95, 0.99, 0.999),
           bw=0.1) {
    # for one signle outline
    if (missing(method)) {
      cat("  * Method not provided. efourier is used.\n")
      method   <- efourier
    } else {
      p <- pmatch(tolower(method), c("efourier", "rfourier", "tfourier"))
      if (is.na(p)) { warning("Unvalid method. efourier is used.")
      } else {
        method   <- switch(p, efourier,   rfourier,   tfourier)}}
    res <- matrix(nrow=length(id), ncol=(nb.h-drop))
    x <- (drop+1) : nb.h
    for (i in seq(along=id)) {
      xf  <- method(Coo@coo[[id[i]]], nb.h = nb.h, smooth.it = smooth.it)
      pow <- harm.pow(xf)[x]
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
      abline(h=lineat.y, lty=2, col="grey90")
      segments(x,    res[1, ], x,    res[3, ], lwd=0.5)
      segments(x-bw, res[1, ], x+bw, res[1, ], lwd=0.5)
      segments(x-bw, res[3, ], x+bw, res[3, ], lwd=0.5)
      lines(x, res[2, ], type="o", pch=20, cex=0.6) 
      if (legend) {
        legend("topright", legend = as.character(probs), bty="o", lwd=1,
               bg="#FFFFFFCC", cex=0.7, inset = 0.005,
               title = "Quantiles")}
      box()
    }
    return(res)})

# Ptolemy ############################################################
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
               col=cols[j], angle=10, length=0.05, lwd=1.2)
      }
    }
    points(0, 0, pch=20, col=cols[1])
    if (legend) {
      legend("topright", legend = as.character(1:nb.h), bty="o",
             col = cols, lty = 1, lwd=1, bg="#FFFFFFCC", cex=0.7,
             title = "Number of harmonics")}
  })


# smooth.qual ########################################################
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
           col = cols, lty = 1, bty = "o", bg="#FFFFFFCC", cex=0.7, title = "Smooth iterations")})
