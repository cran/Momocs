% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cl-handling.R
\name{rename}
\alias{rename}
\title{Renames (ala dplyr) on Momocs objects}
\usage{
rename(.data, ...)
}
\arguments{
\item{.data}{a \code{Coo}, \code{Coe}, \code{PCA} object}

\item{...}{comma separated list of unquoted expressions}
}
\value{
a Momocs object of the same class.
}
\description{
Rename variables by name, from the \code{$fac}. See examples and \code{?dplyr::rename}.
}
\details{
dplyr verbs are maintained.
}
\examples{
olea
rename(olea, Ind=ind, View=view)
}
\seealso{
Other handling functions: \code{\link{arrange}},
  \code{\link{at_least}}, \code{\link{chop}},
  \code{\link{combine}}, \code{\link{dissolve}},
  \code{\link{filter}}, \code{\link{mutate}},
  \code{\link{rm_uncomplete}}, \code{\link{rw_rule}},
  \code{\link{sample_frac}}, \code{\link{sample_n}},
  \code{\link{select}}, \code{\link{slice}},
  \code{\link{subset.Coo}}, \code{\link{transmute}}
}
