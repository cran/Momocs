# fac -------
.print_fac <- function(x, n=6){
  # # remove dim
  # pre <- format(x, n=n) %>% `[`(-1)
  # # remove 'more rows'
  # pre <- pre[-length(pre)]
  # # add dollars in front of column names
  # #pre[1] <- gsub("( )([[:alnum:]]{1})", " $\\2", pre[1])
  # # trim col numbers and add 2 spaces
  # c(paste0("- ", ncol(x), " classifiers (in $fac): "), pre) %>%
  #   gsub("(^\033.*\033\\[39m )", "  ", .) %>%
  #   paste0("  ", .) %>%
  #   # nicely cat
  #   cat(sep="\n")
  print(x, n=n)
}

# # Used in Coo/Coe printers
# .print.fac <- function(fac){
#   nf <- ncol(fac)
#   # here we print the number of classifiers
#   if (nf == 0) {
#     cat(" - $fac: No classifier defined in $fac\n")
#   } else {
#     if (nf<2) {
#       cat(" - $fac:", nf, "classifier:\n")
#     } else {
#       cat(" - $fac:", nf, "classifiers:\n")}
#     # here we print every classifier
#     for (i in 1:nf) {
#       if (is.numeric(fac[, i])){
#         xi <- fac[, i]
#         nas <- sum(is.na(xi))
#         xi  <- xi[!is.na(xi)]
#         xi.sum <- list(min=min(xi), med=median(xi), max=max(xi), mean=mean(xi), sd=sd(xi))
#         xi.sum <- lapply(xi.sum, signif, 3)
#         xi.sum$nas <- nas
#         cat("     '", colnames(fac)[i], "' (numeric): ",
#             #"min:", xi.sum$min,
#             #", med:", xi.sum$med,
#             #", max: ", xi.sum$max,
#             #"mean:", xi.sum$mean,
#             #", sd: ", xi.sum$sd,
#             "mean: ", xi.sum$mean, ", sd: ", xi.sum$sd,
#             ifelse(xi.sum$nas==0, ".\n", paste0(" (", xi.sum$nas, " NA).\n")), sep="")
#       } else {
#         # case where the column is a factor
#         lev.i <- levels(fac[, i])
#         # cosmectics below
#         if (sum(nchar(lev.i))>60){
#           maxprint <- which(cumsum(nchar(lev.i))>30)[1]
#           cat("     '", colnames(fac)[i], "' (factor ", nlevels(fac[, i]), "): ", paste(lev.i[1:maxprint], collapse=", "),
#               " ... + ", length(lev.i) - maxprint, " more.\n", sep="")
#         } else {
#           cat("     '", colnames(fac)[i], "' (factor ", nlevels(fac[, i]), "): ", paste(lev.i, collapse=", "), ".\n", sep="")
#         }
#       }
#     }
#   }
# }

.other_components <- function(x, not=c("coo", "fac")){
  ot <- ls(x)[!(ls(x) %in% not)]
  ot <- paste0("$", ot)
  if (length(ot)>1)
    ot <- paste(ot, sep=", ")
  paste0("  - also: ", ot) %>% cat
}

##### Various utilities

# check and errors -----
# check a condition. if not satisfied, stops with a message
.check <- function(cond_to_pass, msg_if_not) {
  if (!cond_to_pass)
    stop(msg_if_not, call. = FALSE)
}


# can be used is combination with try()
# TRUE/FALSE
.is.error <- function(x) {
  if (is.list(x))
    return(any(sapply(x, .is.error)))
  any(class(x)=="try-error")
}

# path ---------

.lf.auto <- function() {
  p <- file.choose()
  # damn ugly
  p <- strsplit(p, split = "/")
  p <- p[[1]][-length(p[[1]])]
  p <- paste0(p, collapse = "/")
  lf <- list.files(p, full.names = TRUE)
  return(lf)
}

.trim.ext <- function(lf) {
  gsub("\\.[[:alnum:]]+$", "", lf)
}

.trim.path <- function(lf) {
  lf %>%
    strsplit("/") %>%
    sapply(function(x) x[length(x)])
}

.trim.both <- function(lf) {
  lf %>% .trim.path() %>% .trim.ext()
}

# numbers misc ----
.normalize <- function(x, min.x, max.x) {
  # damn long but default arguments are not accepted
  if (missing(min.x))
    min.x <- min(x)
  x <- x - min(x)
  if (missing(max.x))
    max.x <- max(x)
  x <- x/max.x
  return(x)
}


.mprod <- function(m, s) {
  res <- m
  for (i in 1:ncol(m)) {
    res[, i] <- m[, i] * s[i]
  }
  return(res)
}

.which.out <- function(x, conf=1e-4){
  out <- which(dnorm(x, mean(x), sd(x))< conf)
  if(length(out)==0) {
    return(NA)
  } else {
    return(out)}}

# Some vector utilities.
#
# Returns ratio of norms and signed angle between two vectors provided as four
# numeric.
#
# @param r1 the 'real' part of the first vector, i.e. difference in
# x-coordinates.
# @param i1 the 'imaginary' part of the first vector, i.e. difference in
# y-coordinates.
# @param r2 the 'real' part of the second vector, i.e. difference in
# x-coordinates.
# @param i2 the 'imaginary' part of the second vector, i.e. difference in
# y-coordinates.
# @return A list with two components: \code{r.norms} the ratio of (norm of
# vector 1)/(norm of vector 2) and \code{d.angle} the signed angle 'from' the
# first 'to' the second vector.
# @examples
# vecs_param(1, 0, 0, 2)
#
# @export
.vecs_param <- function(r1, i1, r2, i2) {
  x <- c(r1, i1, r2, i2)
  if (!is.numeric(x)) {
    stop("4 numeric must be passed")
  }
  if (length(x) != 4) {
    stop("4 numeric must be passed")
  }
  r.norms <- sqrt((r2^2 + i2^2))/sqrt((r1^2 + i1^2))
  d1 <- sqrt(sum(r1^2 + i1^2))
  d2 <- sqrt(sum(r2^2 + i2^2))
  return(list(r.norms = d1/d2, d.angle = atan2(i2, r2) - atan2(i1,
                                                               r1)))
}


.normalize <- function(x, min.x, max.x) {
  # damn long but default arguments are not accepted
  if (missing(min.x))
    min.x <- min(x)
  x <- x - min(x)
  if (missing(max.x))
    max.x <- max(x)
  x <- x/max.x
  return(x)
}

# rcolors --------
.rcolors2hex <- function(x){
  x <- x %>% grDevices::col2rgb()
  if ("alpha" %in% rownames(x)){
    apply(x, 2, function(.) grDevices::rgb(.[1], .[2], .[3], alpha=.[4], maxColorValue = 255))
  } else {
    apply(x, 2, function(.) grDevices::rgb(.[1], .[2], .[3], maxColorValue = 255))
  }
}

# options like ------

.is_verbose <- function()
  ifelse(options("Momocs_verbose")[[1]], TRUE, FALSE)

# options("verbose"=F)
# .is_verbose()
# options("verbose"=T)
