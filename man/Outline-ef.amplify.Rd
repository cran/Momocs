\name{Outline: ef.amplify}
\alias{coo.amplify}
\alias{ef.amplify}

\title{
Dilation of shape based on elliptical Fourier decomposition.
}

\description{
These two functions calculate dilated and eroded shapes based on elliptical Fourier  decomposition \emph{i.e.} taking into account the shape as a whole. \code{ef} objects can be passed to \code{ef.amplify} or a list or matrix of coordinates to \code{coo.amplify}.
}

\usage{
ef.amplify(ef, amp=rep(0.5, 4))
coo.amplify(coo, amp=rep(0.5, 4), nb.h, draw=FALSE, ...)
}

\arguments{
	\item{ef }{\code{list}. A list containing \eqn{a_n}, \eqn{b_n}, \eqn{c_n} and \eqn{d_n} Fourier coefficients, such as returned by \code{efourier}.}
	\item{coo }{A \code{list} or a \code{matrix} of coordinates.}
	\item{amp }{A vector of \code{numeric}. If \code{amp} is of length 4, the value specify the multiplication factor for \eqn{a_1}, \eqn{b_1}, \eqn{c_1} and \eqn{d_1} ; if only one value is provided, then the multiplication factor will be the same for the four coefficients \eqn{abcd_1}.}
	\item{draw }{\code{logical}. Whether to draw (see \link{coo.draw}) the amplified \code{coo}.}
	\item{nb.h }{\code{integer}. The number of harmonics to calculate.}
	\item{...}{Additional parameters to be passed to \code{coo.draw}.}
}

\value{
\code{ef.amplify} returns the \code{ef} provided but with "amplified" coefficients for the first harmonics. \code{coo.amplify} returns the \code{coo} provided but "amplified" based on elliptical Fourier.
}

\seealso{
\link{efourier} for a description of the elliptical Fourier analysis and \link{Ptolemy} for an illustration of the first ellipse/harmonic defining the shape "amplitude".
}
\examples{
# same amp factor for every coeff.
data(bot)
bot1 <- coo.sample(bot@coo[[1]], 50)
coo.plot(bot1, col=NA)
amp <- seq(0.9, 0.3, -0.1)
for (i in seq(along=amp)) {
  coo.amplify(bot1, amp[i], draw=TRUE, first=FALSE,
              col=NA, border=col.summer(length(amp))[i])}

# random shape, separate amplification
poly <- efourier.shape(nb.h=5, alpha=3, plot=FALSE)
poly2 <- coo.sample(poly, 20)
layout(matrix(1:4, nc=2, byrow=TRUE))
par(oma=rep(0.5, 4), mar=c(2, 2, 3, 2))
coo.plot(poly2, col=NA, main="an")
for (i in seq(-1, 1, 0.5)) {
  coo.amplify(poly2, c(i, 1, 1, 1), draw=TRUE)}

coo.plot(poly2, col=NA, main="bn")
for (i in seq(-1, 1, 0.5)) {
  coo.amplify(poly2, c(1, i, 1, 1), draw=TRUE)}

coo.plot(poly2, col=NA, main="cn")
for (i in seq(-1, 1, 0.5)) {
  coo.amplify(poly2, c(1, 1, i, 1), draw=TRUE)}

coo.plot(poly2, col=NA, main="dn")
for (i in seq(-1, 1, 0.5)) {
  coo.amplify(poly2, c(1, 1, 1, i), draw=TRUE)}

}

\keyword{Outline Analysis}