\name{rFourier}
\alias{rFourier}
\alias{rFourier-methods}
\alias{rFourier,Coo-method}

\title{
Calculates radii variation analysis on \code{Coo} objects.
}

\description{
This method performs a radii variation Fourier analysis on a \code{Coo}-class object, with the specified parameters (number of harmonics and number of smoothing iterations) and returns a \code{Coe}-class object containing harmonic coefficients (usually normalized). It accepts the same arguments as the function \link{rfourier}.
}

\usage{
rFourier(Coo, nb.h= 40, smooth.it = 0, norm=TRUE)
}

\arguments{
	\item{Coo}{The \code{Coo} object.}
	\item{nb.h}{\code{integer}. The number of harmonics to calculate/use.}
	\item{smooth.it}{\code{integer}. The number of smoothing iterations to perform.}
	\item{norm}{\code{logical}. Whether to scale the outlines so that the mean length of the radii used equals 1.}
}

\seealso{
See \link{rfourier} for the mathematical background.
}

\examples{
data(bot)
rFourier(bot)
}

\keyword{radii variation Fourier analysis}


