\name{tps2d}
\alias{tps2d}
\title{Returns the position of interpolated coordinates.}
\description{Returns the position of interpolated coordinates using Thin Plate Splines.}
\usage{
tps2d(M, matr, matt)
}
\arguments{
  \item{M}{original coordinates to be mapped by TPS}
  \item{matr}{Reference configuration matrix}
  \item{matt}{Target configuration matrix}
}

\value{
Interpolated coordinates arranged in a matrix object.
}
\author{
Entirely written by Julien Claude.
Claude, J. (2008) \emph{Morphometrics Using R}, Use R! series, Springer 330 pp.
}

\keyword{Misc functions}
