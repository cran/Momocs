\name{Coe2nef}
\alias{Coe2nef}

\title{
Exports \code{Coe} objects to \code{.nef} files.
}
\description{
\code{Coe2nef} converts \code{Coe} objects to \code{.nef} files. This function is intended to ease data exchange between Momocs and SHAPE suite (see Iwata in the \bold{References} below).
}

\usage{
Coe2nef(Coe, file="nef.nef")
}

\arguments{
  \item{Coe}{A \code{Coe} object.}
  \item{file}{A \code{character} to specify where to write the \code{.nef} file.}
}

\value{
Writes a \code{.nef} file.
}

\seealso{
\link{nef2Coe}, for the reverse operation.
}

\references{
Iwata H, Ukai Y. 2002. SHAPE: a computer program package for quantitative evaluation of biological shapes based on elliptic Fourier descriptors. The Journal of Heredity 93: 384-385.

You can also have a look to the SHAPE's manual distributed with the program suite, that gives a description of the \code{.nef} format.
}

\keyword{Import}
