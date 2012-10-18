\name{pca}
\alias{pca}
\alias{pca-methods}
\alias{pca,Coe-method}

\title{
A wrapper for \code{dudi.pca} on \code{Coe}-objects.
}
\description{
\code{pca} is a wrapper for \code{dudi.pca} in \code{ade4} on \code{Coe}-objects.
}

\usage{
pca(Coe, subset.fac=NULL,
         subset.lev=NULL,
		 row.w,
		 col.w,
		 center=TRUE,
		 scale=FALSE,
		 scannf=FALSE,
		 nf=3)
}
\arguments{
	\item{Coe}{The \code{Coe} object}
	\item{subset.fac}{If non \code{NULL}, a \code{character} to specify the grouping factor on which to perform the PCA.}
	\item{subset.lev}{If non \code{NULL} for \code{subset.fac}, then a \code{character} to specify the level of \code{subset.fac} on which to perform the PCA.}
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
Returns a \code{dudi.pca} object
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
botE <- eFourier(bot)
botP <- pca(botE)
dudi.plot(botP, "type")
beer <- pca(botE, "type", "beer")
dudi.plot(beer)
whisky <- pca(botE, "type", "whisky")
dudi.plot(whisky)

}
\keyword{Multivariate Analysis}