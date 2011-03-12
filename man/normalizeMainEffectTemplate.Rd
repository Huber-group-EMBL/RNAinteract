\name{normalizeMainEffectTemplate}
\alias{normalizeMainEffectTemplate}
\title{
normalize template main effect
}
\description{
Normalize for a spatial main effect of the template genes.
}
\usage{
normalizeMainEffectTemplate(sgi, screen = NULL, channel = NULL)
}
\arguments{
  \item{sgi}{
  An object of class \code{\link{RNAinteract}}.
}
  \item{screen}{
The name of the screen in which the normalization should be applied. If \code{screen = NULL}, the normalization is applied on all screens.
}
  \item{channel}{
The name of the channel in which the normalization should be applied. If \code{channel = NULL}, the normalization is applied on all channels.
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
sgi <- normalizeMainEffectTemplate(sgi)
}
\keyword{ manip }

