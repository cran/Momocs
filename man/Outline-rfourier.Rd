\name{Outline: rfourier}
\alias{rfourier}
\alias{rfourier.i}
\alias{rfourier.shape}

\title{
Fourier Radius variation core functions.
}
\description{
Provide facilities for radius angle Fourier analysis. \code{rfourier} computes Fourier coefficients from a set of coordinates; \code{rfourier.i} calculates a shape given a set of Fourier coefficients computed with \code{rfourier}; \code{rfourier.shape} calculates a "Fourier radius variation shape" given Fourier coefficients (see \code{Details}).
}

\usage{
rfourier(coo, nb.h = 32, smooth.it = 0)
rfourier.i(rf, nb.pts=300, nb.h=length(rf$an))
rfourier.shape(an, bn, nb.h, nb.pts=80, alpha=2, plot=TRUE)
}

\arguments{
	\item{coo}{A \code{list} or \code{matrix} of coordinates.}
	\item{nb.h}{\code{integer}. The number of harmonics to calculate/use.}
	\item{smooth.it}{\code{integer}. The number of smoothing iterations to perform.}
	\item{rf}{A \code{list} with \code{ao}, \code{an} and \code{bn} components, typically as returned by \code{rfourier}.}
	\item{an}{\code{numeric}. A \code{vector} of harmonic coefficient(s).}
	\item{bn}{\code{numeric}. A \code{vector} of harmonic coefficient(s).}
	\item{nb.pts}{\code{integer}. The number of points to calculate.}
	\item{alpha}{\code{numeric}. The power coefficient associated with the (usually decreasing) amplitude of the Fourier coefficients (see the \bold{Details} section in \link{efourier}).}
	\item{plot}{\code{logical}. Whether to plot the shape or only return the list of coordinates.}
}

\details{
Given a closed outline, the radius \eqn{r}, taken as the distance from the outline barycentre and a given point of the outline, can be expressed as a periodic function of the angle \eqn{\theta}. Harmonics from \eqn{0} to \eqn{k} approximate the function \eqn{r(\theta)}:

\deqn{r(\theta)= \frac{1}{2}a_0 +  \sum\limits_{n=1}^{k}a_k\cos(w_k\theta + b_k\sin(w_k\theta)}
with:
\deqn{
a_n = \frac{2}{p}\sum\limits_{n=1}^{p}r_i\cos(n\theta_i)
}
\deqn{
b_n = \frac{2}{p}\sum\limits_{n=1}^{p}r_i\sin(n\theta_i)
}
with
\deqn{
a_0 = \sqrt{\frac{2}{p}}\sum\limits_{n=1}^{p}r_i
}

The \eqn{a_n} and \eqn{b_n} harmonic coefficients, extracted for every individual shape, are then used for multivariate analyses.

}
\value{
Returns a list of coordinates.
}
\references{
Claude, J. (2008) \emph{Morphometrics with R}, Use R! series, Springer 316 pp.
}
\seealso{
\link{efourier}, \link{tfourier}.
}

\examples{
rfourier.shape()
}

\keyword{Outline Analysis}
