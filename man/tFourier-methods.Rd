\name{tFourier}
\alias{tFourier}
\alias{tFourier-methods}
\alias{tFourier,Coo-method}

\title{
Calculates tangent angle analysis on \code{Coo} objects.
}

\description{
\code{tFourier} performs a tangent angle Fourier analysis on a \code{Coo}-class object, with the specified parameters (number of harmonics and number of smoothing iterations) and returns a \code{Coe}-class object containing harmonic coefficients (usually normalized). It accepts the same arguments as the function \link{efourier}.
}

\usage{
tFourier(Coo, nb.h= 40, smooth.it = 0, norm=TRUE)
}

\arguments{
	\item{Coo}{The \code{Coo} object.}
	\item{nb.h}{\code{integer}. The number of harmonics to calculate/use.}
	\item{smooth.it}{\code{integer}. The number of smoothing iterations to perform.}
	\item{norm}{\code{logical}. Whether to scale the outlines so that the mean length of the radii used equals 1.}
}

\seealso{
See \link{tfourier} for the mathematical background.
}

\examples{
data(bot)
tFourier(bot)
}

\keyword{tangent angle Fourier analysis}
