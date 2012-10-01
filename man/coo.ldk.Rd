\name{coo.ldk}
\alias{coo.ldk}

\title{
Utility to identify landmarks on a matrix or list of coordinates.
}

\description{
\code{coo.ldk} allows user to click to identify the closest points that belong to the outline. Typically used to define landmarks.
}

\usage{
coo.ldk(coo, nb.ldk)
}

\arguments{
  \item{coo}{A \code{list} or a \code{matrix} of coordinates.}
  \item{nb.ldk}{\code{integer}. Indicates the number of landmarks to define.}
  }

\value{
A vector of the \code{coo} rows indices that coorespond to the landmarks identified.
}

\seealso{
\link{defLandmarks}.
}

\examples{
\dontrun{
data(bot)
coo <- bot@coo[[1]]
coo.ldk(coo, 2) # you have to click now !
}
}

\keyword{coo Utilities}