\name{efourier}
\alias{efourier}
\title{Computes the Fourier coefficients on a list of coordinates.}
\description{
Computes the Fourier coefficients \eqn{a_o}, \eqn{a_n}, \eqn{b_n}, \eqn{c_o}, \eqn{c_n}, \eqn{d_n} from a list of \eqn{(x; y)} coordinates of the sampled
points.
}
\usage{
efourier(coo, nb.h = 32, smooth.it = NULL)
}
\arguments{
  \item{coo}{the \code{Coo}-object}
  \item{nb.h}{code{integer} indicating how many harmonics to calculate}
  \item{smooth.it}{code{integer} indicating how many smoothing iterations to perform}
}

\value{
  \item{ao}{numeric: the \eqn{a_o} harmonic coefficient}
  \item{co}{numeric: the \eqn{c_o} harmonic coefficient}
  \item{an}{a vector of numeric indicating the \eqn{a_(1 \to n)} harmonic coefficients}
  \item{bn}{a vector of numeric indicating the \eqn{b_(1 \to n)} harmonic coefficients}
  \item{cn}{a vector of numeric indicating the \eqn{c_(1 \to n)} harmonic coefficients}
  \item{dn}{a vector of numeric indicating the \eqn{d_(1 \to n)} harmonic coefficients}
}
\author{
Originally written by Julien Claude.
Claude, J. (2008) \emph{Morphometrics Using R}, Use R! series, Springer 330 pp.
}
\examples{
data(bottles.cont)
efourier(bottles.cont@coo[[1]])
}

\keyword{Misc functions}
