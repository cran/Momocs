\name{Outline: efourier}
\alias{efourier}
\alias{efourier.i}
\alias{efourier.shape}
\alias{efourier.norm}

\title{
Elliptical Fourier Analysis core functions.
}

\description{
These functions provide facilities for elliptical Fourier analysis. \code{efourier} computes Fourier coefficients from a set of coordinates; \code{efourier.i} calculates a shape given a set of Fourier coefficients computed with \code{efourier}; \code{efourier.shape} calculates a "Fourier elliptical shape" given Fourier coefficients (see \code{Details}) ; \code{efourier.norm} normalizes Fourier coefficients for rotation, tranlation, size and orientation of the first ellipse.
}

\usage{
efourier(coo, nb.h = 32, smooth.it = 0)
efourier.i(ef, nb.h, nb.pts = 300)
efourier.shape(an, bn, cn, dn, nb.h, nb.pts=80, alpha=2, plot=TRUE)
efourier.norm(ef, start = FALSE)
}

\arguments{
	\item{coo }{A \code{list} or a \code{matrix} of coordinates.}
	\item{nb.h }{\code{integer}. The number of harmonics to use}
	\item{smooth.it }{\code{integer}. The number of smoothing iterations to perform.}
	\item{nb.pts }{\code{integer}. The number of points to calculate.}
	\item{ef }{\code{list}. A list containing \eqn{a_n}, \eqn{b_n}, \eqn{c_n} and \eqn{d_n} Fourier coefficients, such as returned by \code{efourier}.}
	\item{an }{\code{numeric}. The \eqn{a_n} Fourier coefficients on which to calculate a shape.}
	\item{bn }{\code{numeric}. The \eqn{b_n} Fourier coefficients on which to calculate a shape.}
	\item{cn }{\code{numeric}. The \eqn{c_n} Fourier coefficients on which to calculate a shape.}
	\item{dn }{\code{numeric}. The \eqn{d_n} Fourier coefficients on which to calculate a shape.}
	\item{alpha }{\code{numeric}. The power coefficient associated with the (usually decreasing) amplitude of the Fourier coefficients (see \bold{Details}).}
	\item{plot }{\code{logical}. Whether to plot or not the shape.}
	\item{start }{\code{logical}. Whether to conserve the position of the first point of the outline.}
}
\details{
These functions and their mathematical background detailed below are here detailed to ease their use in new methods but are used internally by methods on \code{Coo}-objects.

Elliptic Fourier analysis and normalization are calculated as follows. Let \eqn{T} be the perimeter of a given closed outline, here considered as the period of the signal. One sets \eqn{\omega = 2\pi/T} to be the pulse. Then, the curvilinear abscissa, \eqn{t} varies from \eqn{0} to \eqn{T}. One can express \eqn{x(t)} and \eqn{y(t)} as:
\deqn{
x(t) = \frac{a_0}{2}+\sum\limits_{n=1}^{+\infty} a_n\cos n\omega t + b_n\sin n\omega t
}
with
\deqn{
a_n = \frac{2}{T}+ \int\limits_{0}^{T} x(t)\cos (n\omega t) \mathrm{d} t
}
\deqn{
b_n = \frac{2}{T}+ \int\limits_{0}^{T} x(t)\sin (n\omega t) \mathrm{d} t }

similarly,
\deqn{
y(t) = \frac{c_0}{2}+\sum\limits_{n=1}^{+\infty} c_n\cos n\omega t + d_n\sin n\omega t
}
with
\deqn{
c_n = \frac{2}{T}+ \int\limits_{0}^{T} y(t)\cos (n\omega t) \mathrm{d} t
}
\deqn{d_n = \frac{2}{T}+ \int\limits_{0}^{T} y(t)\sin (n\omega t) \mathrm{d} t
}

Since the outline contains a \eqn{k} finite number of points, one can therefore calculate discrete estimators for every harmonic coefficient of the \eqn{n^{th}} harmonics:
\deqn{
a_n=\frac{T}{2\pi^2n^2}\sum\limits_{p=1}^k \frac{\Delta x_p}{\Delta t_p}(\cos\frac{2\pi nt_p}{T}-\cos\frac{2\pi nt_{p-1}}{T}) }
\deqn{
b_n=\frac{T}{2\pi^2n^2}\sum\limits_{p=1}^k \frac{\Delta x_p}{\Delta t_p}(\sin\frac{2\pi nt_p}{T}-\sin\frac{2\pi nt_{p-1}}{T})
}

\eqn{\Delta x_1=x_1-x_k} and \eqn{c_n} and \eqn{d_n} are calculated similarly. \eqn{a_0} and \eqn{c_0} correspond to the estimate of the coordinates of the centroid of original outline and are estimated by:
\deqn{
a_0=\frac{2}{T}\sum\limits_{i=1}^p x_i
}
and
\deqn{
c_0=\frac{2}{T}\sum\limits_{i=1}^p y_i
}

Intuitively, for all positive integers \eqn{n}, the sum of a cosine curve and a sine curve represent the \eqn{n^{th}} harmonic content of the \eqn{x} and \eqn{y} projections of the \eqn{k}-edged polygon, and for any \eqn{n}, these two curves define an ellipse in the plane. Ferson and colleagues noticed that in the "time" it takes the \eqn{n^{th}} harmonic to traverse its ellipse \eqn{n} times, the \eqn{(n+1)^{th}} harmonic has traversed its own ellipse \eqn{n+1} times. The reconstruction of the original polygon is done by vector adding these ellipses for all harmonics, which echoes astronomical Ptolemy's epicycles (see \link{Ptolemy}), and the reconstruction obtained from \eqn{N} harmonics is the best possible fit in a least-squares sense.

\code{efourier.shape} can be used to draw "elliptical Fourier shapes" by specifying \code{nb.h} and \code{alpha}. The coefficients are then sampled in an uniform distribution \eqn{(-\pi ; \pi)} and this amplitude is then divided by \eqn{harmonic_rank^alpha}. If \code{alpha} is lower than 1, consecutive coefficients will thus increase.

}

\value{
	\item{an }{\code{vector} of \eqn{A_{1->n}} harmonic coefficients.}
	\item{bn }{\code{vector} of B_{1->n} harmonic coefficients.}
	\item{cn }{\code{vector} of C_{1->n} harmonic coefficients.}
	\item{dn }{\code{vector} of D_{1->n} harmonic coefficients.}
	\item{ao }{\code{ao} Harmonic coefficient.}
	\item{co }{\code{co} Harmonic coefficient.}
	\item{x }{\code{vector} of \code{x}-coordinates.}
	\item{y }{\code{vector} of \code{y}-coordinates.}
	\item{size }{Magnitude of the semi-major axis of the first fitting ellipse.}
	\item{theta }{Angle, in radians, between the starting point and the semi-major axis of the first fitting ellipse.}
	\item{psi }{Orientation of the first fitting ellipse.}
	\item{A }{\code{vector} of \code{numeric} A_{1->n} \emph{normalized} harmonic coefficients.}
	\item{B }{\code{vector} of \code{numeric} B_{1->n} \emph{normalized} harmonic coefficients.}
	\item{C }{\code{vector} of \code{numeric} C_{1->n} \emph{normalized} harmonic coefficients.}
	\item{D }{\code{vector} of \code{numeric} D_{1->n} \emph{normalized} harmonic coefficients.}
	\item{lnef }{A \code{list} with A, B, C and D concatenated in a \code{vector} that may be convenient for some uses.}
}

\references{
Claude, J. (2008) \emph{Morphometrics with R}, Use R! series, Springer 316 pp.

Ferson S, Rohlf FJ, Koehn RK. 1985. Measuring shape variation of two-dimensional outlines. \emph{Systematic Biology} \bold{34}: 59-68.
}
\seealso{
\link{Ptolemy} for an implementation of the Ptolemaic ellipses.
}
\examples{
data(bot)
q <- (efourier(bot@coo[[1]]))
efourier.i(q) # equivalent to efourier.shape(q$an, q$bn, q$cn, q$dn)
efourier.norm(q)
efourier.shape(nb.h=5, alpha=1.2)
efourier.shape(nb.h=12, alpha=0.9)
}

\keyword{Outline Analysis}