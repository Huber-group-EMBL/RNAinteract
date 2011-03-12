\name{sgisubset}
\alias{sgisubset}
\title{
subset of an RNAinteract object.
}
\description{
A new object of class \code{\link{RNAinteract}} is created that contains a subset of screens and channels.
}
\usage{
sgisubset(sgi, screen = getScreenNames(sgi), channel = getChannelNames(sgi))
}
\arguments{
  \item{sgi}{An object of class \code{\link{RNAinteract}}.}
  \item{screen}{
Names of the selected screens.
}
  \item{channel}{
Names of the selected channels.
}
}
\details{
This function returns a \code{\link{RNAinteract}} object that only contains the selected screens and channels.
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
sgi1 <- sgisubset(sgi, screen = "1")
sgi1
sgi2 <- sgisubset(sgi, channel = "nrCells")
sgi2
}
\keyword{ manip }

