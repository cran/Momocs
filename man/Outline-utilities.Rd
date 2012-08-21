\name{Utilities: coeff.sel}
\alias{coeff.sel}
\alias{coeff.split}

\title{
Selection of harmonic coefficients.
}

\description{
Harmonic coefficients are arranged in the \code{@coeff} slot of \code{Nef}-objects in that way:
\deqn{A_1, \dots, A_n, B_1, \dots, B_n, C_1, \dots, C_n, D_1, \dots, D-n}. \code{coeff.sel} helps to select a given number of harmonics and \code{coeff.split} returns a named list of coordinates.
}
\usage{
coeff.sel(retain = 8, drop = 0, nb.h = 32, cph = 4)
coeff.split(cs, nb.h = 8, cph = 4)
}

\arguments{
  \item{cs}{A \code{vector} of harmonic coefficients.}
  \item{retain}{\code{numeric}. The number of harmonics to retain.}
  \item{drop}{\code{numeric}. The number of harmonics to drop}
  \item{nb.h}{\code{numeric}. The maximum harmonic rank.}
  \item{cph}{\code{numeric}. The number of coefficients per harmonic.}
}

\value{
\code{coeff.sel} returns indices that can be used to select columns from an harmonic coefficient matrix.
\code{coeff.split} returns a named list of coordinates.
}
\examples{
coeff.sel(retain=8, drop=0, nb.h=24)
coeff.sel(retain=8, drop=0, nb.h=24, cph=2) # for tangent angle / radius variation matrices.
coeff.split(1:128, nb.h=32, cph=4)
}
\keyword{Utilities}

