\name{plotHeatmap}
\alias{plotHeatmap}
\title{
plots a heatmap for an interaction screen.
}
\description{
A heatmap of an interaction screen is plotted.
}
\usage{
plotHeatmap(sgi, screen, channel, pi.max = NULL,
            main = expression(paste(pi, "-score")),
	    hc.row = NULL, hc.col = NULL, 
	    withoutgroups = c("neg", "pos"))
}
\arguments{
  \item{sgi}{An object of class \code{\link{RNAinteract}}}
  \item{screen}{The screen name of which the interaction matrix is plotted.}
  \item{channel}{The channel name of which the interaction matrix is plotted.}
  \item{pi.max}{The pairwise interaction score that is represented at the top of the color scale. All interaction scores above this value can not be distinguished any more.}
  \item{main}{The title of the plot.}
  \item{hc.row}{A hierarchical clustering (hclust) for the rows.}
  \item{hc.col}{A hierarchical clustering (hclust) for the columns.}
  \item{withoutgroups}{The genes within this group are not shown in the heatmap. It is convinient to hide screen controls in the heatmap.}
}
\details{
A heatmap for one screen and one channel is plotted. Positive interactions are marked blue, negative ones are marked yellow. A colorbar is shown on the left hand side.
}
\value{Returns a grob.}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
plotHeatmap(sgi, screen="1", channel="nrCells")
}
\keyword{ hplot }

