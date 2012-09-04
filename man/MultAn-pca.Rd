\name{Multivariate Analysis: pca}
\alias{pca}
\alias{pca-methods}
\alias{pca,Nef-method}

\title{
A wrapper for \code{dudi.pca} on \code{Nef}-objects.
}
\description{A wrapper for \code{dudi.pca} in \code{ade4} on \code{Nef}-objects.
}
\usage{
pca(Nef, subset=NULL,
		 row.w,
		 col.w,
		 center=TRUE,
		 scale=FALSE,
		 scannf=FALSE,
		 nf=3)
}
\arguments{
	\item{Nef}{The \code{Nef} object}
	\item{subset}{If non \code{NULL}, then the subset on which to perform the PCA. Must be specified as one of the levels of the \code{factor} slot (\code{Nef@fac}.}
	\item{row.w}{A \code{vector} giving the rows weights (see \bold{Details}).}
	\item{col.w}{A \code{vector} giving the cols weights (see \bold{Details}).}
	\item{center}{\code{logical}. Whether to center or not the harmonic coefficients.}
	\item{scale}{\code{logical}. Whether to scale or not the harmonic coefficients.}
	\item{scannf}{\code{logical}. Whether to ask or not the number of PC to retain.}
	\item{nf}{\code{integer}. If \code{scannf} is FALSE, then the number of PC to retain.}
  }

\details{The default parameters should be satisfying for most of the cases. For instance, amplitude of coefficients is not rescaled since first harmonics capture most of the information. See the reference below for further technical discussion on PCA on harmonic coefficients.
}
\value{
Return a \code{dudi.pca} object
}
\references{
See the papers below that introduce \code{ade4} and also the package's homepage: \url{http://pbil.univ-lyon1.fr/ADE-4/}

Dray, S. and Dufour, A.B. (2007): The ade4 package: implementing the duality diagram for ecologists. \emph{Journal of Statistical Software}. \bold{22}(4): 1-20.

Chessel, D. and Dufour, A.B. and Thioulouse, J. (2004): The ade4 package-I- One-table methods. \emph{R News}. \bold{4}: 5-10.

Dray, S. and Dufour, A.B. and Chessel, D. (2007): The ade4 package-II: Two-table and K-table methods. \emph{R News}. \bold{7}(2): 47-52.

}
\seealso{
\link{dudi.pca}.
}
\examples{
data(bot)
botD <- pca(eFourier(bot))
botD # note the new components
}
\keyword{Multivariate Analysis}