\name{Utilities: dev.plot}
\alias{dev.plot}

\title{
Calculates and plots series with associated error bars.
}
\description{This function is used internally by methods based on deviations for one one many outlines. Yet, it provides a quick way to create plots of series, possibly with deviations, from scratch.
}
\usage{
dev.plot(mat, matv, cols, x = 1:ncol(mat), bw = 0.1)
}

\arguments{
  \item{mat}{A \code{matrix} containing one or many lines (as individuals) with the corresponding y values (as cols).}
  \item{matv}{A \code{matrix} of the same dimension as mat but containing the deviation from the \code{mat} matrix.}
  \item{cols}{A \code{vector} of \code{ncol(mat)} colors.}
  \item{x}{\code{integer}. The range on which to plot.}
  \item{bw}{\code{numeric}. The width of the errors bars to draw.}
}

\examples{
foo.m <- matrix(1:10, nr=10, nc=10, byrow=TRUE) + rnorm(100, sd=0.2)
plot(NA, xlim=range(foo.m), ylim=range(foo.m))
dev.plot(foo.m, foo.m*0.1, col.sari(10))
}
\keyword{Utilities}