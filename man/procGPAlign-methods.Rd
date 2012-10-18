\name{procGPAlign(Coo)}
\alias{procGPAlign}
\alias{procGPAlign-methods}
\alias{procGPAlign,Coo-method}

\title{
Performs a Generalized Procrustes Alignement on a \code{Coo} object.
}
\description{
Performs a Generalized Procrustes Alignement on a \code{Coo} object, providing that landmarks have previously been defined.}

\usage{
procGPAlign(Coo, tol=1e-30)
}

\arguments{
  \item{Coo}{The \code{Coo} object.}
  \item{tol}{thresold for GPA. See \link{procGPA}.}
  }

\details{
If landmarks are defined, their coordinates are used for the GPA alignment, with the \link{procGPA} in \code{shapes}. Then the original outlines are rotated and scaled accordingly.
}

\value{
Returns a \code{Coo} with centered, rotated and scaled outlines.
}

\seealso{
\link{defLandmarks}, \link{cooLandmarks}, \link{procGPA}.
}


\examples{
\dontrun{
data(hearts)
# we add a lot of perturbations below
for (i in 1:hearts@coo.nb){
  coo.i <- hearts@coo[[i]] * runif(1, 0, 20)
  coo.i <- coo.rotate(coo.i, runif(1, 0, 2*pi))
  coo.i <- coo.trans(coo.i, runif(1, 0, 200), runif(1, 0, 200))
  hearts@coo[[i]] <- coo.i
}

def.par <- par(no.readonly = TRUE)
layout(matrix(1:3, 1, 3))
stack(hearts, borders="black", ldk.pch=20, ldk.cex=0.5)
title("Troubled hearts")
# we perform GPA
hearts2 <- procGPAlign(hearts)
stack(hearts2, borders="#1A1A1A22", ldk.pch=20, ldk.cex=0.5)
title("After GP alignement")
# we register a new baseline
hearts3 <- baseline(hearts2, 2, 4)
stack(hearts3, borders="#1A1A1A22", ldk.pch=20, ldk.cex=0.5)
title("And with a new baseline")
par(def.par)
}
}
\keyword{Utilities}