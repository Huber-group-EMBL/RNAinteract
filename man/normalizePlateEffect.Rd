\name{normalizePlateEffect}
\alias{normalizePlateEffect}
\title{
Normalization of plate effects
}
\description{
Normalization of plate effects in the screen.
}
\usage{
normalizePlateEffect(sgi, type = "Bscore", maxit = 20, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sgi}{An object of class \code{\link{RNAinteract}}}
  \item{type}{If type is "Bscore" (default) a Bscore-normalization is performed. If type is "spatial", a locfit regression is estimated that accounts for spatial effects.}
  \item{maxit}{Maximum number of iterations for locfit.}
  \item{verbose}{Either 0 (default, no output), 1 (minimum output), or 2 (outout).}
}
\details{
The Bscore normalization estimates row and column effects for each plate. It returns the residuals to the sum of row and column effects. The spatial normalization estimates a non-linear 2D regression for each plate and returns the residuals.
}
\value{
An object of class \code{\link{RNAinteract}}. The returned object contains the normalization information.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
normalizePlateEffect(sgi)
}
\keyword{ manip }

