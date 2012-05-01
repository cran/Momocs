\name{Momocs-package}
\alias{Momocs-package}
\alias{Momocs}
\docType{package}
\title{
Outline Analysis using Elliptical Fourier Analysis.
}
\description{
Momocs is intended to ease and popularize shape analysis of outlines (especially using elliptical Fourier analysis). It mostly hinges on the functions developped in  Morphometrics with R (Claude, 2008). From outline extraction of images and elliptical Fourier calculation to multivariate analysis and the visualization of transformations within the morphological space, Momocs provides a complete and convenient toolkit to specialists within every field that are, or may be, interested in morphological comparisons of outlines.

It comes with its vignette that details step by step how to perform Elliptical Fourier Analysis on a set of shapes, whether starting from images or coordinates in \code{.txt} files.
}

\author{
\enumerate{
\item Vincent Bonhomme, French Institute of Pondicherry, India.
\item Sandrine Picq, UMR CBAE, Montpellier, France.
\item Julien Claude UMR, ISEM, Universite de Montpellier II, France.
}}

\references{
Claude, J. (2008) \emph{Morphometrics with R}, Use R! series, Springer 316 pp.
}

\examples{
\dontrun{
data(bottles.cont)
plot(bottles.cont)
dev.qual(bottles.cont)
dev.quant(bottles.cont)
harm.pow(bottles.cont)
nef <- get.Nef(bottles.cont)
fac <- factor(rep(c("beer","whisky"), each=20))
pca(nef, fac)
pca3(nef, fac)
pca(tps(nef, fac)
from <- c(-0.5, 0.25)
to <- -from
tps.iso(nef, fr=from, to=to)
tps.vf(nef, fr=from, to=to)
tps.grid(nef, fr=from, to=to)
}
}
\keyword{Package}