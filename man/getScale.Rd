\name{getScale}
\alias{getScale}
\title{
get scale information for a channel.
}
\description{
Returns a character string with the scale of each channel.
}
\usage{
getScale(sgi, channel)
}
\arguments{
  \item{sgi}{A \code{\link{RNAinteract}}} object.
  \item{channel}{A channel name.}
}
\value{Returns a character string with scale information for each channel.}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
getScale(sgi, channel="nrCells")
}
\keyword{ attribute }

