\name{col.sel}
\alias{col.sel}
\title{
Helps to select the columns indices of an harmonic coefficient matrix.
}
\description{
Returns the columns' indices of a matrix which colums are in the format:
\deqn{(A_1, ..., A_n, B_1, ..., B_n, C_1, ..., C_n, D_1, ..., D_n)}
}
\usage{
col.sel(h.fr = 1, h.to = 8, h.max = 32, drop = FALSE)
}
\arguments{
  \item{h.fr}{\code{integer} indicating the first harmonic to retain}
  \item{h.to}{\code{integer} indicating the last harmonic to retain}
  \item{h.max}{\code{integer} indicating the total number of harmonics (usually number of columns/4)}
  \item{drop}{\code{logical} indicating whether to drop or not the first \eqn{A_1} harmonic}
}
\value{
Returns a vector of \code{integer} indicating the selected columns indices.
}
\examples{
col.sel(1, 8, 32)
col.sel(1, 4, 8, TRUE)
}

\keyword{Misc functions}
