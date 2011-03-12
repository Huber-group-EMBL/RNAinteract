\name{computePI}
\alias{computePI}
\title{
compute pairwise interaction
}
\description{
Compute the pairwise interaction term for each single experiments.
}
\usage{
computePI(sgi)
}
\arguments{
  \item{sgi}{
  An object of class \code{\link{RNAinteract}}.
}
}
\details{
Computes the pairwise interaction term for each single experiment. Multiple values for each gene pair are not yet summarized.
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
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
sgi <- computePI(sgi)
PI <- getData(sgi, type="pi", format="targetMatrix")
}
\keyword{ manip }

