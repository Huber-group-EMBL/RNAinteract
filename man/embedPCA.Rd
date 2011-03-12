\name{embedPCA}
\alias{embedPCA}
\title{
Computes a PCA for a pairwise interaction matrix.
}
\description{
A principal component analysis is performed for a pairwise interaction matrix. The low-dimensional embedding is returned.
}
\usage{
embedPCA(sgi, screen, channel,
         dim = 4, embed = "template",
         withoutgroups = c())
}
\arguments{
  \item{sgi}{An object of class \code{\link{RNAinteract}}}
  \item{screen}{The screen name whose interaction matrix will be embedded.}
  \item{channel}{The channel name whose interaction matrix will be embedded.}
  \item{dim}{The embedding dimension.}
  \item{embed}{Either "template" (default) or "query" denotes if the embedding is done for rows or columns. }
  \item{withoutgroups}{Genes annotated with these groupnames are not considered for embedding.}
}
\value{
Returns a matrix with dimensions genes x dim.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
X <- embedPCA(sgi, screen="1", channel="nrCells", dim=2)
plot(X[,1], X[,2], pch=20, cex=0.01)
text(X[,1], X[,2], row.names(X))
}
\keyword{ manip }

