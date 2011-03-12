\name{getReplicateData}
\alias{getReplicateData}
\alias{getIndDesignData}
\title{
Extract replicates measurements from the screen.
}
\description{
A genetic interaction screen can contain within screen replicates, if some reagent pairs are measured at least twice. Usually this appears when measuring reagent pairs once as template-query and once as query-template. \code{getReplicateData} returns a list of these technical replicates.

If multiple reagents are used to target the same gene, different reagent pairs that target the same gene pair are extracted from the screen. These pairs are returned by \code{getIndDesignData}.
}
\usage{
getReplicateData(sgi, screen, channel,
                 type = "data", design = "template",
                 do.trafo = TRUE, do.inv.trafo = FALSE,
                 normalized = FALSE)
getIndDesignData(sgi, screen, channel,
                 type = "data", design = "template",
		 do.trafo = TRUE, do.inv.trafo = FALSE,
		 normalized = FALSE)
}
\arguments{
  \item{sgi}{An object of class \code{\link{RNAinteract}}.}
  \item{screen}{The screen name from which the replicates will be extracted.}
  \item{channel}{The channel name from which the replicates will be extracted.}
  \item{type}{The type of data that is extracted. It is the type argument of the \code{\link{getData}} function.}
  \item{design, do.trafo, do.inv.trafo, normalized}{See the \code{\link{getData}} documentation for details.}
}
\value{Returns a data.frame with columns x and y. }
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
res <- getIndDesignData(sgi, screen="1", channel="nrCells", type = "data")
plot(res$x, res$y)
}
\keyword{ attribute }

