\name{createRNAinteractFromFiles}
\alias{createRNAinteractFromFiles}
\title{
create an RNAinteract object from text files
}
\description{
Reads text files with annotation, query and template design, and data. Creates a RNAinteract object.
}
\usage{
createRNAinteractFromFiles(name = "anonymous",
			  filePlatelist = "Platelist.txt",
		   	  fileReagents = "Reagents.txt",
		   	  fileTargets = "Targets.txt",
		   	  fileTemplateDesign = "TemplateDesign.txt",
		   	  fileQueryDesign = "QueryDesign.txt",
		   	  path = ".",
		   	  pdim = NULL,
		   	  Transformation = "log2")
}
\arguments{
  \item{name}{A name for the screen.}
  \item{filePlatelist}{The filename of the text file containing the plate list.}
  \item{fileReagents}{The filename of the text file containing the reagent annotation.}
  \item{fileTargets}{The filename of the text file containing the target annotation.}
  \item{fileTemplateDesign}{The filename of the text file containing the template design.}
  \item{fileQueryDesign}{The filename of the text file containing the query design.}
  \item{path}{The system directory were the textfiles \code{filePlatelist}, \code{fileReagents}, \code{fileTargets}, \code{fileTemplateDesign}, and \code{fileQueryDesign} are located.}
  \item{pdim}{The dimensions of the multi-well plates (e.g. \code{pdim = c(nrow=24,ncol=16)}). if \code{pdim =NULL} (default), the plate dimension will be estimated from the input data.}
  \item{Transformation}{The transformation that is applied to the data. All calculations are done on additive scale.}
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
\code{\link{RNAinteract-package}}, \code{\link{createRNAinteract}}, \code{\link{createCellHTSFromFiles}}
}
\keyword{ manip }

