\name{estimateMainEffect}
\alias{estimateMainEffect}
\title{
estimate main effect
}
\description{
estimates the main effects in an additive model. 
}
\usage{
estimateMainEffect(sgi, use.query = NULL)
}
\arguments{
  \item{sgi}{
  An object of class \code{\link{RNAinteract}}.
}
  \item{use.query}{
A list of reagent identifiers as annotated in the RID field of the reagent list. For the estimation of the template main effects only these queries are used.
}
}
\details{
The main effect is the single RNAi knockdown effect. When use.query is not specified, the main effect is estimated by minimizing the L1 distance from the non-interacting model to the double RNAi measurements. The implemented non-interacting model is the additive model (sum of single main effects). If the screen does not contain a lot of query genes with no or ery small main effect, it is recommended to estimate the template main effects only by using selected query genes. This can be obtained by specifying use.query. To estimate main effects in a multiplicative model define \code{Transformation="log2"} when creating the RNAinteract object (See \code{\link{createRNAinteractFromFiles}}), which is already the default.
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
sgi <- estimateMainEffect(sgi)
getMain(sgi)
}
\keyword{ manip }

