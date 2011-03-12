\name{createRNAinteract}
\alias{createRNAinteract}
\title{
create a RNAinteract object
}
\description{
Creates a RNAinteract object given data matrices, annotation, query and template design.
}
\usage{
createRNAinteract(data, well, plate, pdim,
          	  Reagents, Targets,
                  TemplateDesign, QueryDesign,
                  Transformation = NULL)
}
\arguments{
  \item{data}{An array with dimensions features x screens x channels.}
  \item{well}{A vector of length #features with well names.}
  \item{plate}{A vector of length #features with plate numbers.}
  \item{pdim}{A vector of length 2 with plate dimensions (e.g. pdim=c(12,8)).}
  \item{Reagents}{A data.frame describing the reagents.}
  \item{Targets}{A data.frame describing the targets.}
  \item{TemplateDesign}{A data.frame with the layout of the template plates.}
  \item{QueryDesign}{A data.frame with the layout of the query plates.}
  \item{Transformation}{A Transformation that is applied to the data. If NULL the data is log2 transformed.}
}
\details{
See vignette("RNAinteract") for an example how to create an RNAinteract object.
}
\value{
An object of class \code{\link{RNAinteract}}.
}
\references{
 ~put references to the literature/web site here ~
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}, \code{\link{createRNAinteractFromFiles}}, \code{\link{createCellHTSFromFiles}}
}
\keyword{ manip }

