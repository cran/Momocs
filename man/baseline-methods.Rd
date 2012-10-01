\name{baseline(Coo)}
\alias{baseline}
\alias{baseline-methods}
\alias{baseline,Coo-method}

\title{
Re-register a \code{Coo} on a new baseline.
}

\description{
\code{edm} returns the euclidean distances between points \deqn{1 -> n} of two 2-col matrices of the same dimension. This function is used internally but may be of interest for other analyses.
}
\usage{
baseline(Coo, ldk1=1, ldk2=2, t1=c(-0.5, 0), t2=c(0.5, 0))
}

\arguments{
  \item{Coo}{The \code{Coo} object from which to retrieve coordinates.}
  \item{ldk1}{index of the first reference landmark.}
  \item{ldk2}{index of the second reference landmark. See Examples.}
  \item{t1}{\eqn{(x; y)} coordinates of the first target point.}
  \item{t2}{\eqn{(x; y)} coordinates of the second target point.}
  }

\value{
Returns a \code{Coo} registered in the new baseline.
}

\seealso{
\link{ed}, \link{edm.nearest}, \link{dist}.
}


\examples{
def.par <- par(no.readonly = TRUE)
layout(matrix(1:2, 1, 2))
data(hearts)
stack(hearts, border="#1A1A1A22")
hearts2 <- baseline(hearts, 2, 4)
stack(hearts2, border="#1A1A1A22")
par(def.par)

}
\keyword{Utilities}