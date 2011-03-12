\name{sgisubsetQueryDesign}
\alias{sgisubsetQueryDesign}
\title{
Subsetting query genes in a RNAinteract object.
}
\description{
A \code{\link{RNAinteract}} object with a subset of query genes is returned.
}
\usage{
sgisubsetQueryDesign(sgi, query.targets = NULL, query.reagents = NULL)
}
\arguments{
  \item{sgi}{An \code{\link{RNAinteract}} object.}
  \item{query.targets}{The query target ids to be selected.}
  \item{query.reagents}{The query reagent ids to be selected.}
}
\value{
Returns an oject of class \code{\link{RNAinteract}}.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
% \examples{
% data("sgi")
% sgi
% sginew <- sgisubsetQueryDesign(sgi, query.targets = c("Jra", "Sos", "Ras85D", "Ctrl_Fluc"))
% sginew
% }
\keyword{ manip }

