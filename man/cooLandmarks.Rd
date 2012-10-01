\name{cooLandmarks(Coo)}
\alias{cooLandmarks,Coo-method}
\alias{cooLandmarks}
\alias{cooLandmarks-methods}

\title{
Retrieves landmarks coordinates from a \code{Coo} object.
}

\description{
Retrieves landmarks (if any) coordinates from a \code{Coo} object.
}
\usage{
cooLandmarks(Coo)
}
\arguments{
  \item{Coo}{The \code{Coo} object.}}

\value{
Returns an array of landmark coordiantes such as described in \link{A.plot}.
}


\seealso{
\link{procGPAlign}, \link{defLandmarks}, \link{A.plot}.
}


\examples{
data(hearts)
ldks <- cooLandmarks(hearts)
A.plot(ldks)

}
\keyword{Utilities}