\name{getMain}
\alias{getMain}
\alias{getMainNeg}
\title{
get main effects
}
\description{
Returns the main effects.
}
\usage{
getMain(sgi, type = "main", design = "template", summary = "none",
        QueryNr = NULL, TemplatePlate = NULL, 
	do.inv.trafo = FALSE, format = "plain", withoutgroups = c(), 
	screen = NULL, channel = NULL, normalized = TRUE, drop = TRUE)
getMainNeg(sgi, type = "all", do.inv.trafo = FALSE, 
	screen = NULL, channel = NULL, drop = TRUE)
}
\arguments{
  \item{sgi}{An object of class \code{\link{RNAinteract}}}
  \item{type}{always "main"}
  \item{design}{Either "template" or "query" defining if template or query main effects are returned.}
  \item{summary}{If summary is "targets" the main effects are summarized per target gene.}
  \item{QueryNr,TemplatePlate}{Onle main effects of one query nr or one template plate are returned. }
  \item{format}{targetmatrix}
  \item{withoutgroups}{The genes within this group are not shown in the heatmap. It is convenient to hide screen controls.}
  \item{do.inv.trafo}{If TRUE, the data will be back-transformed for original scale of data. In the case of log-transformed data, the main effects are returned as factors, otherwise the main effects are returned as log values.}
  \item{screen}{The screen from which the main effects should be returned.}
  \item{channel}{The channel from which the main effects should be returned.}
  \item{drop}{Does return a drop array dimensions, even if only one screen or one channel is selected.}
  \item{normalized}{If true the normalized main effects are returned.}
}
\value{
An array containing the main effects.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
getMain(sgi)
getMainNeg(sgi)
}
\keyword{ attribute }

