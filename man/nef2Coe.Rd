\name{nef2Coe}
\alias{nef2Coe}

\title{
Imports \code{.nef} files and creates \code{Coe} objects.
}
\description{
\code{nef2Coe} converts \code{.nef} files to \code{Coe} objects. This function is intended to ease data exchange between Momocs and SHAPE suite (see Iwata in the \bold{References} below).
}

\usage{
nef2Coe(nef.path)
}

\arguments{
  \item{nef.path}{A \code{character} that indicates the path for the \code{.nef} file to convert.}
}

\value{
\code{nef2Coe} returns a \code{Coe} object.
}

\seealso{
\link{Coe2nef}, for the reverse operation.
}

\references{
Iwata H, Ukai Y. 2002. SHAPE: a computer program package for quantitative evaluation of biological shapes based on elliptic Fourier descriptors. The Journal of Heredity 93: 384-385.

You can also have a look to the SHAPE's manual distributed with the program suite, that gives a description of the \code{.nef} format.
}

\keyword{Import}
