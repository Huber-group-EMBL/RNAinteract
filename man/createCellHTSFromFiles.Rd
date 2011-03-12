\name{createCellHTSFromFiles}
\alias{createCellHTSFromFiles}
\title{
create cellHTS2 object from text files
}
\description{
A cellHTS2 object is created from a set of text files.
}
\usage{
createCellHTSFromFiles(filePlatelist = "Platelist.txt",
                       name = "anonymous",
		       path = ".", pdim = NULL)
}
\arguments{
  \item{filePlatelist}{The platelist (See vignette("RNAinteract")) }
  \item{name}{A (arbtrary) string providing the name for the screen}
  \item{path}{The path were the data files are located}
  \item{pdim}{Giving the plate dimensions, e.g. pdim=c(8,12) or pdim=c(16,24).}
}
\details{
See vignette("RNAinteract") for an example how to create an RNAinteract object.
}
\value{
An object of class \code{\link{cellHTS2}}.
}
\references{
 ~put references to the literature/web site here ~
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}, \code{\link{createRNAinteractFromFiles}}, \code{\link{createRNAinteract}}
}
\keyword{ manip }

