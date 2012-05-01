\name{draw.Fell}
\alias{draw.Fell}
\title{Draws "Fourier Ellipses".}
\description{Calculates and draws a "Fourier Ellipse" corresponding to the harmonic coefficients provided}
\usage{
draw.Fell(an = pi, bn = -pi, cn = pi, dn = pi,
			n = 200, cols = topo.colors, title = FALSE)
}
\arguments{
  \item{an}{a \code{numeric} corresponding to the \eqn{a_n} harmonic coefficient}
  \item{bn}{a \code{numeric} corresponding to the \eqn{b_n} harmonic coefficient}
  \item{cn}{a \code{numeric} corresponding to the \eqn{c_n} harmonic coefficient}
  \item{dn}{a \code{numeric} corresponding to the \eqn{d_n} harmonic coefficient}
  \item{n}{\code{integer} indicating how many points to retrieve from outline reconstruction}
  \item{cols}{a color palette such as \code{\link{topo.colors}}, or those produced by \code{colorRampPalette}}
  \item{title}{\code{integer} indicating whether to add a title to the plot}
}
\examples{
draw.Fell()
draw.Fell(2*pi, -pi, pi, 3*pi, title=TRUE)
}
\keyword{Misc functions}
