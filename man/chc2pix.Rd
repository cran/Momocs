\name{chc2pix}
\alias{chc2pix}

\title{
Converts chain-coded coordinates to a matrix of coordinates.
}
\description{
\code{chc2pix} converts chain-coded coordinates such as those used in the SHAPE suite (see Iwata in the \bold{References} below) to a matrix of coordinates.
}

\usage{
chc2pix(chc)
}
\arguments{
  \item{chc}{A \code{numeric} that corresponds to a chain-coded outline.}
  }
\details{
Chain-code is a coding system for describing outlines in numbers from 0 to 7. Given the \eqn{i^{th}} pixel taken on an outline and considering its eight neighbors, the chain-code equivalent for describing the relative position of the \eqn{(i+1)^{th}} pixel is a number from \eqn{0} (the next cell is eastward) to, counting counter clockwise, \eqn{7} (the next cell is south-eastward). See \bold{References} and \bold{Examples} below. 
}

\value{
Returns a 2-col \code{matrix} that corresponds to the \eqn{(x; y)} coordinates encoded in the chain provided.
}

\seealso{
\link{pix2chc}, for the reverse operation.
}

\references{
Freeman H. 1974. Computer processing of line-drawing images. ACM Computing Surveys (CSUR) 6: 57-97.

Iwata H, Ukai Y. 2002. SHAPE: a computer program package for quantitative evaluation of biological shapes based on elliptic Fourier descriptors. The Journal of Heredity 93: 384-385.

Kuhl FP, Giardina CR. 1982. Elliptic Fourier features of a closed contour. Computer Graphics and Image Processing 18: 236-258.

You can also have a look to the SHAPE's manual distributed with the program suite, that gives a description of the \code{.chc} format.
}
\examples{
data(bot)
coo <- bot@coo[[1]]
chc <- pix2chc(coo)
coo.plot(chc2pix(chc))

# Illustration of chain coding
plot(NA, xlim=c(0, 3), ylim=c(0, 3), axes=FALSE, ann=FALSE, xaxs="i", yaxs="i")
title("Position of the next pixel and corresponding chain-code")
abline(h=0:3, v=0:3)
rect(1, 1, 2, 2, col="grey80")
text(1.5, 1.5, "Starting\npixel", cex=2)
text(x=c(2.5, 2.5, 1.5, rep(0.5, 3), 1.5, 2.5),
	y=c(1.5, rep(2.5, 3), 1.5, rep(0.5, 3)), labels=0:7, cex=2)
}

\keyword{Import}
