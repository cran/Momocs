\name{panel.lm}
\alias{panel.lm}
\title{Calculates and plots confidence ellipses.}
\description{
Given a set of \eqn{(x; y)} coordinates, calculates confidence ellipses and plots them.
}
\usage{
panel.lm(x, y, r = 1, col = "black", lwd = 1, lty = 1)
}
\arguments{
  \item{x}{a vector of numeric \eqn{x} coordinates}
  \item{y}{a vector of numeric \eqn{y} coordinates}
  \item{r}{a numeric indicating the number of standard deviations to calculates confidence ellipses}
  \item{col}{code{integer} or a character indicating the color of the ellipses to draw}
  \item{lwd}{a numeric indicating the lwd to use when drawing ellipses}
  \item{lty}{code{integer} indicating the lty to use when drawing ellipses}
}
\examples{
plot(x <- rnorm(50), y <- rnorm(50))
panel.lm(x,y)
}
\keyword{Misc functions}
