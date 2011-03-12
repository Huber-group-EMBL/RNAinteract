\name{summarizeScreens}
\alias{summarizeScreens}
\title{
summarize screens
}
\description{
Creates a new object of class \code{\link{RNAinteract}} with one screen. The new screen is the mean of all screens in the input object.
}
\usage{
summarizeScreens(sgi, screens, newscreenname = "mean")
}
\arguments{
  \item{sgi}{
  An object of class \code{\link{RNAinteract}}.
}
  \item{screens}{
The screen names to be summarized.
}
  \item{newscreenname}{
The name of the new summary screen.
}
}
\details{
If multiple screens with the same layout are stored in the same \code{\link{RNAinteract}} object, then these screens are summarized by averaging to a new screen. The returned object contains one screen.
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
sgi
sginew <- summarizeScreens(sgi, screens=c("1","2"), newscreenname = "m")
sginew
}
\keyword{ manip }

