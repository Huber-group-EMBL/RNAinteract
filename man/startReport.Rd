\name{startReport}
\alias{startReport}
\alias{endReport}
\title{
start and end a RNAinteract report
}
\description{
\code{startReport} will open a html page and starts writing an html report for a RNAinteract screen. \code{endReport} finishes the report and closes the html-file.
}
\usage{
startReport(outputpath)
endReport(report)
}
\arguments{
  \item{outputpath}{the path to the output directory were the report is written to.}
  \item{report}{An report object as returned by startReport or any report... function.}
}
\details{
~~ details ~~
}
\value{
startReport returns an report object. It is handed over to each report-function.
}

\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}},
\code{\link{reportAnnotation}}, \code{\link{reportStatistics}}, \code{\link{reportGeneLists}}, \code{\link{reportNetworks}}, \code{\link{reportScreenData}}, \code{\link{reportDoublePerturbation}}, \code{\link{reportMainEffects}}
}
\examples{
data("sgi")
report <- startReport("report")
reportAnnotation(sgi, report = report)
endReport(report)
browseURL(file.path("report","index.html"))
}
\keyword{ print }

