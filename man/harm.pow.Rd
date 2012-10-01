\name{harm.pow}
\alias{harm.pow}
\title{Calculates harmonic power given a list from e/t/rfourier
}
\description{
Given a list with \code{an, bn (and eventually cn and dn)}, returns the harmonic power.
}

\usage{
harm.pow(xf)
}
\arguments{
  \item{xf}{A list with an, bn (and cn, dn) components, typically from a e/r/tfourier passed on coo.}
}
\value{
Returns a \code{vector} of harmonic power
}

\examples{
data(bot)
ef <- efourier(bot@coo[[1]], 24)
rf <- efourier(bot@coo[[1]], 24)
harm.pow(ef)
harm.pow(rf)

plot(cumsum(harm.pow(ef)[-1]), type="o",
  main="Cumulated harmonic power without the first harmonic",
  ylab="Cumulated harmonic power", xlab="Harmonic rank")

}

\keyword{Utilities}