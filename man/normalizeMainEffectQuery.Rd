\name{normalizeMainEffectQuery}
\alias{normalizeMainEffectQuery}
\title{
normalize query main effect
}
\description{
Normalize for a time effect of the query genes.
}
\usage{
normalizeMainEffectQuery(sgi, batch = NULL, time = NULL)
}
\arguments{
  \item{sgi}{
  An object of class \code{\link{RNAinteract}}.
}
  \item{batch}{
  batch is a vector if integers with length equal to the number of queries. It assigns each query to a batch. Within each batch a linear regression is estimated assuming a linear effect between the order of queries and the main effects.
}
  \item{time}{
  batch is a vector of numbers. A linear regression is estimated fitting the main effect as a function of the time.
}
}
\details{
Normalizing the query main effect does not influence the estimation of the pairwise interaction term.
}
\value{
  An object of class \code{\link{RNAinteract}}.
}

\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
sgi <- normalizeMainEffectQuery(sgi)
}
\keyword{ manip }

