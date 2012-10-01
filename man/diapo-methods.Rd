\name{diapo(Coo)}
\alias{diapo}
\alias{diapo-methods}
\alias{diapo,Coo-method}

\title{
Plots a slideshow of each outline in a \code{Coo} object.
}

\description{
\code{diapo(Coo)} plots a slideshow of each outline in a \code{Coo} object with graphical options.
}

\usage{
diapo(Coo, id=1, col=NA, border="#708090",
      ldk=TRUE, ldk.pch=3, ldk.col="red", ldk.cex=1, ...)
}

\arguments{
  \item{Coo}{The \code{Coo} object to plot.}
  \item{id}{The id of the Coo to start with.}
  \item{col}{A color for drawing the outline.}
  \item{border}{A color for drawing the border.}
  \item{ldk}{\code{logical}. Whether to display landmarks (if any).}
  \item{ldk.pch}{A \code{pch} for these landmarks.}
  \item{ldk.col}{A color for these landmarks.}
  \item{ldk.cex}{A \code{cex} for these landmarks.}
  \item{...}{Additional arguments to be passed to \link{coo.draw}.}
}

\value{
No returned value.
}

\seealso{
\link{stack}, \link{diapo}.
}

\examples{
\dontrun{
data(mosquito)
diapo(hearts, 230, col="pink", border="black")
}
}
\keyword{Utilities}