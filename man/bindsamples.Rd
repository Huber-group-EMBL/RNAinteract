\name{bindscreens}
\alias{bindscreens}
\title{
bind RNAinteract objects along screens
}
\description{
Bind two RNAinteract objects along screens.
}
\usage{
bindscreens(sgi1, sgi2)
}
\arguments{
  \item{sgi1}{
  An object of class \code{\link{RNAinteract}}.
}
  \item{sgi2}{
  An object of class \code{\link{RNAinteract}}.
}
}
\details{
This function binds two double interaction screens along screens.
}
\value{
  An object of class \code{\link{RNAinteract}} with all screens in sgi1 and sgi2.
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
sgi
sginew <- summarizeScreens(sgi, screens=c("1","2"), newscreenname = "m")
sginew
sgibind <- bindscreens(sgi, sginew)
sgibind
}
\keyword{ manip }

