\name{Outline: tfourier}
\alias{tfourier}
\alias{tfourier.i}
\alias{tfourier.shape}

\title{Fourier Tangent angle core functions.}
\description{
Provide facilities for tangent angle Fourier analysis. \code{tfourier} computes Fourier coefficients from a set of coordinates; \code{tfourier.i} calculates a shape given a set of Fourier coefficients computed with \code{tfourier}; \code{tfourier.shape} calculates a "Fourier tangent angle shape" given Fourier coefficients (see \code{Details}).
}

\usage{
tfourier(coo, nb.h=32, smooth.it=0)
tfourier.i(tf, nb.pts, nb.h=length(tf$an), thetao=0)
tfourier.shape(an, bn, ao = 0, nb.h, nb.pts = 80, alpha = 2, plot = TRUE)
}

\arguments{
	\item{coo}{A list or matrix of coordinates}
	\item{nb.h}{\code{integer}. The number of harmonics to calculate/use}
	\item{smooth.it}{\code{integer}. The number of smoothing iterations to perform}
	\item{tf}{a list with ao, an and bn components, typically as returned by tfourier}
	\item{an}{\code{numeric}. A \code{vector} of harmonic coefficient(s).}
	\item{bn}{\code{numeric}. A \code{vector} of harmonic coefficient(s).}
	\item{ao}{\code{numeric}. radius angle to the reference (in radians)}
	\item{nb.pts}{\code{integer}. The number of points to calculate}
	\item{thetao}{\code{numeric}. radius angle to the reference (in radians)}
	\item{alpha}{\code{numeric}. The power coefficient associated with the (usually decreasing) amplitude of the Fourier coefficients (see \link{efourier} - Details)}
	\item{plot}{\code{logical}. Whether to plot the shape or only return the list of coordinates.}
}

\details{
Given a closed outline which the outline has been scaled to \eqn{2\pi}, \eqn{\phi(t)} can be expressed as follows:
\deqn{
\phi(t) = \theta(t) - \theta(0) - t
}
where \eqn{t} is the distance along the outline, \eqn{\theta(t)} the angle of the tangent vector at \eqn{t} and \eqn{\theta(0)} the angle of the tangent vector taken for the first point. It can be removed for normalizing the coefficients obtained. Two coefficients per harmonics can be estimated as follow:

\deqn{
a_n = \frac{2}{p}\sum\limits_{n=1}^{p}\phi(t)\cos n \theta_i
}
\deqn{
b_n = \frac{2}{p}\sum\limits_{n=1}^{p}\phi(t)\sin n \theta_i
}
with
\deqn{
a_0 = \sqrt{\frac{2}{p}}\sum\limits_{n=1}^{p}\phi(t)
}
}

\value{
Returns a list of coordinates.
}
\references{
Zahn CT, Roskies RZ. 1972. Fourier Descriptors for Plane Closed Curves. IEEE Transactions on Computers C-21: 269-281.

Claude J. 2008. Morphometrics with R. Springer Verlag.
}

\seealso{
\link{efourier}, \link{rfourier}.
}

\examples{
data(bot)
tfourier(bot@coo[[1]])
tfourier.shape()
}

\keyword{Outline Analysis}

