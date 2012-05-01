\name{tps}
\alias{tps}
\title{Produces deformation grids.}
\description{
Produces deformation grids using Thin Plate Splines
}
\usage{
tps(matr, matt, n, plot = TRUE, col = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matr}{the reference configuration \code{matrix}}
  \item{matt}{the target configuration \code{matrix}}
  \item{n}{\code{integer} indicating the number of displayed column cells}
  \item{plot}{\code{logical} indicating whether to plot the grid}
  \item{col}{\code{integer} or a character indicating the grid color}
}
\author{
Entirely written by Julien Claude.
Claude, J. (2008) \emph{Morphometrics Using R}, Use R! series, Springer 330 pp.
}
\value{
Deformation grid obtained by Thin Plate Splines interpolation.
}

\keyword{Misc functions}
