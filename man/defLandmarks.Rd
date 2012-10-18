\name{defLandmarks}
\alias{defLandmarks}
\alias{defLandmarks-methods}
\alias{defLandmarks,Coo-method}

\title{Define landmarks on a \code{Coo} object}
\description{
\code{defLandmarks} helps to define landmarks on a \code{Coo} object. The number of landmarks must be specified and rows indices that correspnd to the nearest points clicked on every outlines are stored in the \code{@ldk} slot of the \code{Coo} object.
}

\usage{
defLandmarks(Coo, nb.ldk)
}

\arguments{
	\item{Coo}{A \code{Coo} object on which to define landmarks.}
	\item{nb.ldk}{An \code{integer} taht indicates the number of landmarks to define along the \code{Coo} object.}
}

\seealso{
\link{coo.ldk}.
}

\keyword{Outline Analysis}

