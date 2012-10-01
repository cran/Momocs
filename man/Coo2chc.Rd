\name{Coo2chc}
\alias{Coo2chc}

\title{
Exports \code{Coo} objects to \code{.chc} files that can be used in SHAPE.
}

\description{
\code{Coo2chc} converts \code{Coo} objects into \code{.chc} files that can then be used in the SHAPE suite (see Iwata & Ukai in the \bold{References} below).
}

\usage{
Coo2chc(Coo, file="chc.chc")
}

\arguments{
  \item{Coo}{A \code{Coo} object.}
  \item{file}{A \code{character} to specify where to write the \code{.chc} file.}
  }

\details{
See \link{pix2chc} for an illustration of the functions behind and chain-coding of outlines.}

\value{
Writes a \code{.chc} file.
}

\seealso{
\link{Coe2nef} and \link{nef2Coe}.
}

\references{
Freeman H. 1974. Computer processing of line-drawing images. ACM Computing Surveys (CSUR) 6: 57-97.

Iwata H, Ukai Y. 2002. SHAPE: a computer program package for quantitative evaluation of biological shapes based on elliptic Fourier descriptors. The Journal of Heredity 93: 384-385.

Kuhl FP, Giardina CR. 1982. Elliptic Fourier features of a closed contour. Computer Graphics and Image Processing 18: 236-258.

You can also have a look to the SHAPE's manual distributed with the program suite, that gives a description of the \code{.chc} format.
}


\keyword{Import}
