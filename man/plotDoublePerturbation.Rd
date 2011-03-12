\name{plotDoublePerturbation}
\alias{plotDoublePerturbation}
\title{
Double Perturbation Plot
}
\description{
These function draws a double perturbation plot for interaction screens. All interactions of one gene are displayed in one panel. The double perturbation readout level is plotted against the single perturbation level.
}
\usage{
plotDoublePerturbation( sgi, screen, channel, target,
                       withoutgroups = c("neg", "pos"), design,
                       main, xlab, ylab, range,
                       show.labels = "none", label.par, label,
		       avoid.overlap, col, fill,
                       D , MT, MQ, PV, QV, PI, ...)
}
\arguments{
  \item{sgi}{
An object of class \code{RNAinteract}.
}
  \item{target}{
A character name of the target gene.
}
  \item{screen}{
The character name of the screen to display. If not specified, the first screen is used. Does not have to be specified, if \code{sgi} contains only one screen.
}
  \item{channel}{
The character name of the channel to display. If not specified, the first channel is used. Does not have to be specified, if \code{sgi} contains only one channel.
}
  \item{withoutgroups}{
Interactions to genes from these groups (as specified in the reagent or target annotation) are excluded from the plot, e.g. positive and negative controls.
}
  \item{design}{
  The Either \code{"template"} (default) or \code{"query"}. The single perturbation effects are either the template main effects or the query main effects.
}
  \item{main}{An overall title of the plot.}
  \item{xlab}{A title of the x-axis.}
  \item{ylab}{A title of the y-axis.}
  \item{range}{A numeric vector of length two. range equals the \code{xlim}, \code{ylim} argument in \link{plot}.}
  \item{show.labels}{Automatically select text labels for the points. \code{'all'} shows a text label for all genes, \code{"q.value"} and \code{"p.value"} show a text label for all genes with a q.value (p.value) larger than \code{label.par}, \code{"none"} does not show any text label. This argument has no effect, if \code{label} is specified.}
  \item{label.par}{
Cut-off value for q.value or p.value for displaying text labels (See \code{show.labels}).
}
  \item{label}{
Either a character vector with gene names, or a named vector of text labels. The names of the vector represent the gene names.
}
  \item{avoid.overlap}{If \code{TRUE} (default), text is moved such that text labels are not overlapping.}
  \item{col}{
A named vector with colors. The names of \code{col} define which points are colored (See also \code{fill}).
}
  \item{fill}{
A list up to four values.
\code{colors} defines a set of colors from which a colorramp is created. If \code{colramp} is specified, colors has no effect.
\code{colramp} directly specifies the colorramp.
\code{values} define the values that are color coded. If \code{values} is not specified, the pairwise interaction term is used instead.
\code{at} is a numeric vector defining the breakpoints along the values. If not specified, breakpoints are selected to range three times the standard deviation of the values around zero.
\code{fill} has no effect, if \code{col} is specified.
}
  \item{D,MT,MQ,PV,QV,PI}{Internal usage.}
  \item{\dots}{Further argument passed to \code{\link{grid.doublePerturbation}} or \code{\link{doublePerturbationGrob}}.}
}
\details{
Plots a double perturbation plot. It shows the interaction profile for one (query) gene.
}
\value{
  A grob is returned.
}

\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}, \code{\link{grid.doublePerturbation}}, \code{\link{reportDoublePerturbation}}
}
\examples{
data("sgi")
plotDoublePerturbation( sgi, screen="1", channel="nrCells", target="rl", show.labels="p.value")
}
\keyword{ hplot }

