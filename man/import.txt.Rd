\name{import.txt}
\alias{import.txt}

\title{
Everything to convert \code{.txt} files to a list of coordinates.

}
\description{
Takes a path list to \code{.txt} files and returns a list of coordinates arranged as matrices that can be passed to the \code{Coo} builder.
}
\usage{
import.txt(txt.list, ...)
}
\arguments{
  \item{txt.list}{a list of full/relative path to \code{.txt} files. (see \bold{Details}).}
  \item{...}{Complementary arguments to be passed to \code{read.table} within the function.}
}
\details{
\code{.txt} files must have coordinates arranged in two columns (\emph{i.e.} x and y) with no header. The list of paths pointing to \code{.txt} files must be \emph{full} if your working directory is different from your folder containing images, see \bold{Examples} below.
}
\value{A \code{list} of matrices of coordinates.
}
\seealso{
\link{import.jpg}.
}
\examples{
\dontrun{
txt.list <- list.files(path_to_your_folder_containing_.txt_files, full=TRUE)
# altenartively :
# setwd(path_to_your_folder_containing_.txt_files)
# txt.list <- list.files()
I <- import.txt(txt.list)
Coo(I)
I # that should be a Coo object
}
}
\keyword{Import}