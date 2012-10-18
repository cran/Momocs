\name{coo.slide}
\alias{coo.slide}

\title{
"Slides" a list or a matrix of coordinates.
}

\description{
\code{coo.slide} "slides" \code{coo}, \emph{i.e.} makes the \code{id1} become the first and change the others accordingly.
}

\usage{
coo.slide(coo, id1)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{id1}{\code{integer}. Specifies the index of the coordinates that has to be defined as the first coordinate.}
 }

\examples{
m <- matrix(1:10, 5, 2)
m
coo.slide(m, 3)
}
\keyword{Utilities}