\name{computePValues}
\alias{computePValues}
\title{
compute p-values
}
\description{
Compute p-values for genetic interactions terms. Assess if genetic interaction term is different from zero.
}
\usage{
computePValues(sgi,
               method = "pooled.ttest",
               mixTemplateQuery = TRUE,
               p.adjust.function = function(x) { p.adjust(x, method = "BH")},
               verbose = 0)
}
\arguments{
  \item{sgi}{
  An object of class \code{\link{RNAinteract}}.
}
  \item{method}{
  The method used to compute p-values. One of "pooled.ttest","ttest", "limma", "HotellingT2".

  For "ttest" a Student t-test is applied for each gene pair. The variance is estimated locally for each gene pair. 
  For "pooled.ttest", a pooled variance is estimated from all gene pairs. The variance applied for each gene pair is the maximum of the pooled and the local variance estimate. This method obtains conservative p-values.
  For "limma" mediates between the local and the global variance estimation in a Bayesian framework. The \link{limma-package} is applied to compute the p-values.
  For "HotellingT2" Hotelling-T^2 statistics is computed jointly for all dimensions. It results in a single p-value summarizing all channels. For simplification the p-values are stored in a matrix of dimension genes x genes x screens x channels and the p-values are repeated for each channel. The same holds for q-values.
}
  \item{mixTemplateQuery}{
  If a gene-pair is measured twice as template-query and as query-template, a single p-value is computed by combining all measurements, if \code{mixTemplateQuery =  TRUE}.
  Else a p-value is computed independently for both cases.
}
  \item{p.adjust.function}{
  A function that corrects the p-values for multiple testing. Default method is the Benjamini-Hochberg method.
}
  \item{verbose}{Either 0 (default, no output), 1 (minimum output), or 2 (outout)}
}
\details{
 Computes p-values from a t-test, using the bioconductor package limma, or with a multidimensional Hotelling T^2 test.  
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
sgi <- computePValues(sgi, method = "HotellingT2")
# Hotelling T^2 test will provide one p-value for all channels, PV will be the same
# for all channels in this case
PV <- getData(sgi, type="p.value", format="targetMatrix", channel="nrCells")
}
\keyword{ manip }
\keyword{ htest }
\keyword{ multivariate }

