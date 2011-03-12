\name{grid.doublePerturbation}
\alias{grid.doublePerturbation}
\alias{doublePerturbationGrob}
\alias{postDrawDetails.doublePerturbation}
\title{
Double Perturbation Plot Grob
}
\description{
These functions create a double perturbation grob for interaction screens. All interactions of one gene are displayed in one panel. The double perturbation readout level is plotted against the single perturbation level.
}
\usage{
doublePerturbationGrob( mainEffect, dpEffect, mainEffectTarget,
                       range=NULL, main=NULL, xlab=NULL, ylab=NULL,
                       text=NULL, avoid.overlap=TRUE,
                       axisOnOrigin = FALSE,
                       drawBox = TRUE,
		       pch = 21, size=unit(1, "char"), fill = NULL,
                       gpMain = gpar(lty="dashed", lwd=3, col="cyan"),
                       gpNI = gpar(lty="dashed", lwd=3, col="orange"),
                       gpPoints = gpar(pch=21),
                       gpText = NULL,
                       gpAxis = NULL,
                       gpWTLines=NULL,
                       name=NULL, gp=NULL, vp=NULL )
grid.doublePerturbation(..., draw = TRUE)

# a helper function for doublePerturbationGrob:
postDrawDetails.doublePerturbation(x, recording)
}
\arguments{
  \item{mainEffect}{
A numeric vector of main effects.
}
  \item{dpEffect}{
A numeric vector of double perturbation effects.
}
  \item{mainEffectTarget}{
The main effect of the target gene (A single numeric value).
}
  \item{range}{
The range of the plot. Equals the \code{xlim}, \code{ylim} arguments of \code{\link{plot}}.
}
  \item{main}{An overall title of the plot.}
  \item{xlab}{A title of the x-axis.}
  \item{ylab}{A title of the y-axis.}
  \item{text}{A character vector of text. Has to have the same length as \code{mainEffect}.}
  \item{avoid.overlap}{If \code{TRUE} (default) the text labels are moved such that the text is not overlapping.}
  \item{axisOnOrigin}{If \code{TRUE}, the x- and y-axis are draw on the origin of the data. If \code{FALSE} (default), the axis are drawn on the left and on the bottom.}
  \item{drawBox}{If \code{TRUE} (default), a box is drawn around the plot.}
  \item{pch}{Either an integer specifying a symbol or a single character to be used in plotting points. See \code{\link{points}} for possible values.}
  \item{size}{A unit object specifying the size of the plotting symbols.}
  \item{fill}{
A list containing (some of) the following elements:
\code{col} defines a fill color for the points. Either a single value or a vector of the same length as \code{mainEffect}. If \code{col} is defined, all other elements of fill have no effect.
\code{values} is a numeric vector of the same length as \code{mainEffcet} that contains values that are mapped to colors. 
\code{at} is a numeric vector indicating breakpoints along the values. If not specified will be equally spaced on the range of the values.
\code{colors} defines a set of colors from which a colorramp is created.
\code{colramp} defines a colorramp directly. \code{colramp} has no effect, if \code{colors} is defined.
}
  \item{gp}{
An object of class \code{gpar}, typically the output from a call to the function \code{\link{gpar}}. This is basically a list of graphical parameter settings. Overall settings for the plot are set in \code{gp}. 
}
  \item{gpMain,gpNI}{An object of class \code{gpar} (See \code{gp}). \code{gpMain} and \code{gpNI} indicate the graphics parameter for the main effect lines and the non-interacting line. }
  \item{gpPoints,gpText,gpAxis,gpWTLines}{An object of class \code{gpar} (See \code{gp}). These arguments define graphical parameters for single compartments of the plot. }
  \item{name}{A character identifier.}
  \item{vp}{A Grid viewport object (or NULL).}
  \item{draw}{If \code{TRUE} the grob is drawn on the current device.}
  \item{\dots}{Further arguments passed to \code{\link{doublePerturbationGrob}}.}
  \item{x, recording}{Internal usage only.}
}
\details{
This function creates a grob for a double perturbation plot. It is probably more convenient to use the function \code{\link{plotDoublePerturbation}}.
}
\value{
  A grob is returned.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}, \code{\link{plotDoublePerturbation}}, \code{\link{reportDoublePerturbation}}
}
\keyword{ hplot }

