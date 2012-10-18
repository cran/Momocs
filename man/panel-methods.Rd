\name{panel(Coo)}
\alias{panel,Coo-method}
\alias{panel}
\alias{panel-methods}

\title{
Displays a one-page plot of a \code{Coo}-object.
}

\description{
\code{panel(Coo)} displays a one-page plot of a \code{Coo}-object.
}

\usage{
panel(Coo, cols, borders, names=NULL, ...)
}

\arguments{
  \item{Coo}{The \code{Coo} object to plot.}
  \item{cols}{A \code{vector} of colors for drawing the outlines. Either a single value or of length exactly equals to the number of coordinates.}
  \item{borders}{A \code{vector} of colors for drawing the borders. Either a single value or of length exactly equals to the number of coordinates.}
  \item{names}{\code{logical}. Whether to display names.}
  \item{...}{Additional arguments to be passed to \link{coo.list.panel}.}
}

\value{
No returned value.
}

\seealso{
 \link{coo.list.panel} is behind the scene. See also \link{stack}, \link{diapo}.
}


\examples{
data(mosquito)
panel(mosquito, borders="grey20", col="grey90")

data(hearts)
panel(hearts, cols=rep(col.summer(8), each=30), names=TRUE)
}
\keyword{Utilities}