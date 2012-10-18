\name{chc2Coo}
\alias{chc2Coo}

\title{
Imports \code{.chc} files from SHAPE to \code{Coo} objects that can be used in Momocs.
}

\description{
\code{chc2Coo} converts \code{.chc} files, \emph{e.g.} from the SHAPE suite (see Iwata & Ukai in the \bold{References} below) that can then be used in Momocs.
}

\usage{
chc2Coo(chc.path)
}

\arguments{
  \item{chc.path}{A \code{character} that indicates the path for the \code{.chc} file to convert.}
  }

\details{
See \link{pix2chc} for an illustration of the functions behind and chain-coding of outlines.}

\value{
Returns a \code{Coo} object.
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
