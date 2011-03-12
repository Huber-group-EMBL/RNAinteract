\name{getScreenNames}
\alias{getScreenNames}
\alias{getChannelNames}
\title{get names of screens and channels}
\description{
Returns the names of all screens or all channels.
}
\usage{
getScreenNames(sgi)
getChannelNames(sgi)
}
\arguments{
  \item{sgi}{\code{\link{RNAinteract}}}
}
\value{
A vector of screen or channel names.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")
getScreenNames(sgi)
getChannelNames(sgi)
}
\keyword{ attribute }

