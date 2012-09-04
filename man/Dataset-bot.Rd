\name{Dataset: bot}
\alias{bot}
\docType{data}
\title{Bottles dataset}
\description{
The \code{bot} dataset is a \code{Coo}-class object that contains outline coordinates of 20 beer and 20 whisky bottles.
}
\usage{
data(bot)
}
\format{
A \code{Coo}-object that contains in the slot \code{@coo}, the lists of \eqn{(x; y)} coordinates and in the slot \code{@fac} the corresponding \code{factor}s.
}
\source{
Images have been grabbed on the internet and prepared by the package's authors. No particular choice has been made on the dimension of the original images or the brands cited here.
}

\examples{
data(bot)
bot
}
\keyword{Datasets}
