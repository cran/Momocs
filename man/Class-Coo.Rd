\name{Class: Coo}
\Rdversion{1.1}
\docType{class}
\alias{Coo}
\alias{Coo-class}
\alias{show,Coo-method}
\alias{sample}
\alias{sample-methods}
\alias{Coo.template}
\alias{Coo.template-methods}
\alias{Coo.template,Coo-method}
\alias{Coo.center}
\alias{Coo.center-methods}
\alias{Coo.center,Coo-method}
\alias{Coo.sample}
\alias{Coo.sample-methods}
\alias{Coo.sample,Coo-method}
\alias{Coo.align}
\alias{Coo.align-methods}
\alias{Coo.align,Coo-method}
\alias{Coo.slide}
\alias{Coo.slide-methods}
\alias{Coo.slide,Coo-method}
\alias{Coo.smooth}
\alias{Coo.smooth-methods}
\alias{Coo.smooth,Coo-method}
\alias{Coo.close}
\alias{Coo.close-methods}
\alias{Coo.close,Coo-method}
\alias{Coo.unclose}
\alias{Coo.unclose-methods}
\alias{Coo.unclose,Coo-method}

\alias{[<-,Coo,ANY,ANY,ANY-method}
\alias{[,Coo,ANY,ANY,ANY-method}


\title{Class \code{"Coo"}}
\description{A class for handling morphometric datasets, \emph{i.e.} lists of \eqn{(x; y)} outline coordinates with or without a grouping structure.
}
\section{Objects from the Class}{
New \code{Coo}-objects can be created with \code{new("Coo", ...)} or more directly using the \code{Coo} builder.
}
\section{Slots}{
  \describe{
    \item{\code{coo}:}{\code{list} a list of shapes arranged in two columns coordinate matrices and that can have different lengths.}
    \item{\code{names}:}{\code{character} \code{vector} to name every shape.}
    \item{\code{fac}:}{\code{data.frame} defining the grouping structure, if any. Columns must be \code{factor}s.}
    \item{\code{ldk}:}{\code{list} defining the coordinates indices that will be considered as homologous points (not yet used in current version).}
    \item{\code{coo.nb}:}{\code{numeric}. The number of shapes.}
    \item{\code{coo.len}:}{a \code{vector} of \code{numeric} that contains the number of coordinates per shape.}
    \item{\code{coo.closed}:}{a \code{vector} of \code{logical} that indicates if outlines are closed.}
    \item{\code{details}:}{any additional details, such as licence and authors can be added in this \code{list}.}
  }
}
\section{Methods}{

See \code{showMethods(class="Coo")}.

}
\examples{
## Working zone...
showClass("Coo")
# an example of the Coo builder
data(bot)
list_of_coordinates <- bot@coo
Coo(list_of_coordinates)
}
\keyword{Classes}