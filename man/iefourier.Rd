\name{iefourier}
\alias{iefourier}
\title{
Calculates inverse Fourier Elliptical.
}
\description{
Calculates inverse Fourier Elliptical if passed with harmonic coefficients.
}
\usage{
iefourier(an, bn, cn, dn, k, n, ao = 0, co = 0)
}
\arguments{
  \item{an}{a vector of \code{numeric} indicating the \eqn{a_{1 \to n}} harmonic coefficient}
  \item{bn}{a vector of \code{numeric} indicating the \eqn{b_{1 \to n}} harmonic coefficient}
  \item{cn}{a vector of \code{numeric} indicating the \eqn{c_{1 \to n}} harmonic coefficient}
  \item{dn}{a vector of \code{numeric} indicating the \eqn{d_{1 \to n}}harmonic coefficient}
  \item{k}{\code{integer} indicating the number of harmonics to claculate}
  \item{n}{\code{integer} indicating the number of points to sample on the calculated outline}
  \item{ao}{\code{numeric}: the \eqn{a_0} harmonic coefficient}
  \item{co}{\code{numeric}: \eqn{c_0} harmonic coefficient}
}
\value{
a \code{list} of \code{(x; y)} coordinates.
}
\author{
Entirely written by Julien Claude.
Claude, J. (2008) \emph{Morphometrics Using R}, Use R! series, Springer 330 pp.
}
\keyword{Misc functions}
