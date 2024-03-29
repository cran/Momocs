% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cl-def-Out.R
\name{OutCoe}
\alias{OutCoe}
\title{Builds an OutCoe object}
\usage{
OutCoe(coe = matrix(), fac = dplyr::tibble(), method, norm)
}
\arguments{
\item{coe}{\code{matrix} of harmonic coefficients}

\item{fac}{(optional) a \code{data.frame} of factors,
specifying the grouping structure}

\item{method}{used to obtain these coefficients}

\item{norm}{the normalisation used to obtain these coefficients}
}
\value{
an \code{OutCoe} object
}
\description{
In Momocs, \code{OutCoe} classes objects are wrapping around
lists of morphometric coefficients, along with other informations,
on which generic methods such as plotting methods (e.g. \link{boxplot})
and specific methods can be applied.
\code{OutCoe} objects are primarily \code{\link{Coe}} objects.
}
\details{
These methods can be applied on \code{Out} objects:
}
\examples{
# all OutCoe methods
methods(class='OutCoe')
}
\seealso{
Other classes: 
\code{\link{Coe}()},
\code{\link{Coo}()},
\code{\link{Ldk}()},
\code{\link{OpnCoe}()},
\code{\link{Opn}()},
\code{\link{Out}()},
\code{\link{TraCoe}()}
}
\concept{classes}
