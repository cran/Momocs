\name{get.cont}
\alias{get.cont}
\alias{Conte}
\title{Extract \eqn{(x; y)} coordinates and create a \code{Coo}-object.}
\description{
Extracts from a set of black and white images or a list of coordinates written in a set of 2-columns ("x" and "y") \code{.txt} files.
}
\usage{
get.cont(path)
}
\arguments{
  \item{path}{a path to indicate where are the images or the \code{.txt} files to use}
}
\details{
If no \code{path} is provided, the user is interactively asked to choose a folder.

\code{get.cont} uses the \code{Conte} algorithm that starts on the center of every outline in the \code{Coo}-object provided (or the \code{imagematrix} provided if \code{Conte()} is directly used). If this point does not correspond to a black pixel, \emph{i.e.} not contained within the shape, the user is interactively asked to select interactively a point within the shape.
}
\value{
a \code{Coo}-object is returned.
}
\author{
\code{Conte} was originally written by Julien Claude.
Claude, J. (2008) \emph{Morphometrics Using R}, Use R! series, Springer 330 pp.
}
\examples{
\dontrun{
data(bottles.cont)
get.cont()
}
}

\keyword{Misc functions}
