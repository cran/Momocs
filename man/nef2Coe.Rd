\name{nef2Coe}
\alias{nef2Coe}

\title{
Imports \code{.nef} files and creates \code{Coe} objects.
}
\description{
\code{nef2Coe} converts \code{.nef} files to \code{Coe} objects. This function is intended to ease data exchange between Momocs and SHAPE suite (see Iwata in the \bold{References} below).
}

\usage{
nef2Coe(nef.path)
}

\arguments{
  \item{nef.path}{A \code{character} that indicates the path for the \code{.nef} file to convert.}
}

\value{
\code{nef2Coe} returns a \code{Coe} object.
}

\seealso{
\link{Coe2nef}, for the reverse operation.
}

\references{
Iwata H, Ukai Y. 2002. SHAPE: a computer program package for quantitative evaluation of biological shapes based on elliptic Fourier descriptors. The Journal of Heredity 93: 384-385.

You can also have a look to the SHAPE's manual distributed with the program suite, that gives a description of the \code{.nef} format.
}

\examples{
\dontrun{
# I (VB) just finished my post-doc in India and I wont have time for Momocs
# until mid march or maybe a bit after that.
# Ryan Felice found a bug in my own version. I provide hereafter the function
# he kindly sent to me and that worked fine for him.
# Hope this help. Thanks again to him.
NEF2COE<-function (nef.path){
  nef <- readLines(nef.path)
  HARMO.l <- grep(pattern = "HARMO", nef)
  nb.h <- as.numeric(substring(nef[HARMO.l], 8))
  nef <- nef[-(1:HARMO.l)]
  nb.coo <- length(nef)/(nb.h + 1)
  coo.i <- 1:nb.coo
  coo.beg <- (coo.i - 1) * (nb.h + 1) + 1
  coo.end <- coo.beg + nb.h
  res <- matrix(NA, nrow = nb.coo, ncol = nb.h * 4,
      dimnames = list(nef[coo.beg], paste(rep(LETTERS[1:4], each = nb.h),
      1:nb.h, sep = "")))
  for (i in seq(along = coo.i)) {
    nef.i <- nef[(coo.beg[i]+1):coo.end[i]]
    x <- as.numeric(unlist(strsplit(nef.i, " ")))
    x1<-x[!is.na(x)]
    a.i<-x1[seq(1,length(x1),4)]
    b.i<-x1[seq(2,length(x1),4)]
    c.i<-x1[seq(3,length(x1),4)]
    d.i<-x1[seq(4,length(x1),4)]
    res[i, ]<-c(a.i,b.i,c.i,d.i)
  }
  return(Coe(res,method="eFourier"))}
}
}

\keyword{Import}
