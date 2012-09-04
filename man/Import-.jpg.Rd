\name{Import: .jpg}
\alias{import.jpg}
\alias{import.multi1.jpg}
\alias{import.img.prepare}
\alias{import.img.Conte}

\title{
Everything to convert images to a list of coordinates.
}
\description{Functions to import \code{.jpg} images and convert them to list of coordinates.
}
\usage{
import.jpg(jpg.list)
import.multi1.jpg(path)
import.img.prepare(path)
import.img.Conte(img, x, auto=TRUE, plot=TRUE)
}
\arguments{
  \item{jpg.list}{A \code{vector} of \code{character} containing the path to your \code{.jpg} images.}
  \item{path}{\code{character}. A single path to a \code{.jpg} image.}
  \item{img}{an \code{imagematrix} object.}
  \item{x}{A \code{vector} of \eqn{(x; y)} coordinates from where to start Conte algorithm.}
  \item{auto}{\code{logical}. Whether to try or not to start at the center of the image(s) before asking the user to click within the shape.}
  \item{plot}{\code{logical}. Whether to plot or not the image. Used internally by \code{import.multi1.jpg} to not reload the same image.}
}
\details{
Typically, an object returned by \link{list.files} on a folder containing your images is passed to \code{import.jpg}. \code{import.img.prepare} and \code{import.img.Conte} are typically not used by front-user but internally by \code{import.jpg}. They clean, binarize, threshold, etc. raw .jpg images and extract a list of coordinates from a black and white \code{imagematrix}, respectively. The best option is to work with black and white \code{.jpg} image with the black mask of the outline overlapping the center of the image.
\code{import.multi1.jpg} helps to extract several outlines from the same \code{.jpg} image.
}
\value{
\code{import.jpg} returns a list of (from 1 to thousands) \eqn{(x; y)} coordinates arranged as matrices and that can be then converted to a \code{Coo}-object. \code{import.multi1.jpg} returns a \code{list} of \eqn{(x; y)} coordinates. \code{import.img.prepare} returns an \code{imagematrix} object, \code{import.img.Conte} returns a matrix of \eqn{(x; y)} coordinates.
}
\seealso{
link{import.txt}.
}
\examples{
\dontrun{
jpg.list <- list.files(path_to_your_folder_containing_.txt_files)
I <- import.jpg(your path)
I # that should be a Coo object
}
}
\keyword{Import}