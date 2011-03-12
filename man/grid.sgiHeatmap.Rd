\name{grid.sgiHeatmap}
\alias{grid.sgiHeatmap}
\title{
A heatmap grob
}
\description{
A grob is created and printed for a matrix PI which is intended to represent pairwise interaction scores. 
}
\usage{
grid.sgiHeatmap(PI, pi.max = NULL, main = expression(paste(pi, "-score")), 
		hc.row = NULL, hc.col = NULL)
}
\arguments{
  \item{PI}{A matrix of pairwise interactions.}
  \item{pi.max}{The interaction score at the top end of the colorbar. pairwise interaction score larger than this value can not be distinguished anymore.}
  \item{main}{A title for the plot.}
  \item{hc.row}{An hierarchical clustering as produced by hclust of the rows.}
  \item{hc.col}{Clustering of the columns.}
}
\details{
A heatmap is plotted with positive interaction represented in yellow and negative interactions represented in blue. A colorbar is plotted on the left and dendrograms are added. This function can be used to integrate the plot in other grid objects. It is recommended to use the function \code{\link{plotHeatmap}} to plot heatmaps of an \code{\link{RNAinteract}} object.
}
\value{A grob is returned.}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
PI = getData(sgi, type="pi", format="targetMatrix", screen="1", channel="nrCells")
grid.sgiHeatmap(PI)
}
\keyword{ hplot }

