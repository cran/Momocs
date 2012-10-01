\name{meanShapes}
\alias{meanShapes}
\alias{meanShapes-methods}
\alias{meanShapes,Coe-method}

\title{Calculate groups mean shapes from \code{Coe}}
\description{
\code{meanShapes} calculates mean shapes of groups defined through the \code{@fac} slot and using the appropriate inverse Fourier method. 
}

\usage{
meanShapes(Coe,
           fac,
           nb.pts=300)
}

\arguments{
  \item{Coe}{A \code{Coe} object on which to define landmarks.}
  \item{fac}{One of the \code{factor} names of \code{Coe@fac}. See examples below.}
  \item{nb.pts}{The number of points to use for calculation shapes.}
}

\examples{
data(bot)
botF <- eFourier(bot)

meanShapes(botF) # average shape
shp.type <- meanShapes(botF, "type")
panel(Coo(shp.type), borders=col.gallus(2)) # use Coo plotting facilities !

# below we combine another Fourier approach, one more factor, and one panel plot with names.
bot@fac <- data.frame(type=bot@fac[, 1], plop=rep(letters[1:4], each=10))
botR <- rFourier(bot)
shp.plop <- meanShapes(botR, "plop")
panel(Coo(shp.plop), cols=col.summer(4), names=TRUE)
}
\keyword{Outline Analysis}

