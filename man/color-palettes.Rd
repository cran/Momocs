\name{Color palettes}
\alias{col.summer}
\alias{col.sari}
\alias{col.gallus}
\alias{col.blackgallus}
\alias{col.bw}
\alias{col.bcol}
\alias{col.wcol}
\alias{col.india}

\title{
Some color palettes.
}
\description{
\code{col.summer}, \code{col.india}, \code{col.sari} and \code{col.gallus} are standard color palettes, obtained with \code{colorRampPalette}. \code{col.bw} generates a palette from black to white. \code{col.bcol} and \code{col.wcol} are used to generate a palette from one color to black or white, respectively.
}
\usage{
col.summer(n)
col.sari(n)
col.gallus(n)
col.blackgallus(n)
col.india(n)
col.bw(n)
col.bcol(col.hex)
col.wcol(col.hex)
}
\arguments{
  \item{n}{\code{integer}. The number of colors to create.}
  \item{col.hex}{An hexadecimal \code{character} string coding for a color.}
}
\value{
A \code{vector} of hexadecimal color that can be passed to many methods and functions of \code{Momocs} and also used as any palette created with \link{colorRampPalette}.
}

\examples{
pie(rep(1, 12), col=col.summer(12), main="col.summer")
pie(rep(1, 12), col=col.india(12),  main="col.india")
pie(rep(1, 12), col=col.sari(12),   main="col.sari")
pie(rep(1, 12), col=col.gallus(12), main="col.gallus")
pie(rep(1, 12), col=col.blackgallus(12), main="col.blackgallus")
pie(rep(1, 12), col=col.bw(12), main="col.bw")
pie(rep(1, 12), col=col.bcol(col.hex="#FF6600")(12), main="col.bcol and some orange")
pie(rep(1, 12), col=col.wcol(col.hex="#000080")(12), main="col.bcol and some blue")
barplot(1:33, col=col.summer(33))
}
\keyword{Utilities}