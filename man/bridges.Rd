% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/babel-bridges.R
\name{bridges}
\alias{bridges}
\alias{l2m}
\alias{m2l}
\alias{d2m}
\alias{m2d}
\alias{l2a}
\alias{a2l}
\alias{a2m}
\alias{m2a}
\alias{m2ll}
\title{Convert between different classes}
\usage{
l2m(l)

m2l(m)

d2m(d)

m2d(m)

l2a(l)

a2l(a)

a2m(a)

m2a(m)

m2ll(m, index = NULL)
}
\arguments{
\item{l}{\code{list} with x and y coordinates as components}

\item{m}{\code{matrix} of (x; y) coordinates}

\item{d}{\code{data.frame} with two columns}

\item{a}{\code{array} of (x; y) coordinates}

\item{index}{\code{numeric}, the number of coordinates for every slice}
}
\value{
the data in the required class
}
\description{
Convert between different classes
}
\note{
\code{a2m}/\code{m2a} change, by essence, the dimension of the data.
\code{m2ll} is used internally to hanle coo and cur in \code{Ldk} objects but may be
useful elsewhere
}
\examples{
# matrix/list
wings[1] \%>\% coo_sample(4) \%>\%
   m2l() \%T>\% print \%>\%        # matrix to list
   l2m()                       # and back

# data.frame/matrix
wings[1] \%>\% coo_sample(4) \%>\%
   m2d() \%T>\% print \%>\%        # matrix to data.frame
   d2m                         # and back

 # list/array
 wings \%>\% slice(1:2) \%$\%
 coo \%>\% l2a \%T>\% print \%>\%    # list to array
 a2l                           # and back

 # array/matrix
 wings \%>\% slice(1:2) \%$\%
 l2a(coo) \%>\%                  # and array (from a list)
 a2m \%T>\% print \%>\%            # to matrix
 m2a                           # and back

 # m2ll
m2ll(wings[1], c(6, 4, 3, 5)) # grab slices and coordinates
}
\seealso{
Other bridges functions: 
\code{\link{as_df}()},
\code{\link{complex}},
\code{\link{export}()}
}
\concept{bridges functions}
