\name{RNAinteract-class}
\Rdversion{1.1}
\docType{class}
\alias{RNAinteract-class}
\alias{show,RNAinteract-method}

\title{Class "RNAinteract"}
\description{
A class for double perturbation experiments (genetic interaction screens, drug-drug interaction screens). There are functions for creation, analysis, and display of interaction screens.
}
\section{Objects from the Class}{
Objects can be created by calls of \code{\link{createRNAinteractFromFiles}}. See vignette("RNAinteract") for an example of creating an \code{RNAinteract} object.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"array"} with dimension \code{sgi@F x sgi@S x sgi@C}. The raw data of the screen. }
    \item{\code{screenNames}:}{Object of class \code{"character"} with length \code{sgi@S}. }
    \item{\code{channelNames}:}{Object of class \code{"character"} with length \code{sgi@C}.}
    \item{\code{well}:}{Object of class \code{"character"} with length \code{sgi@F}. Well name (e.g. F04) for each measurement.}
    \item{\code{plate}:}{Object of class \code{"integer"} with length \code{sgi@F}. Number of the plate for each measurement}
    \item{\code{pdim}:}{Object of class \code{"integer"} of length 2. Plate dimensions.}
    \item{\code{NT}:}{Object of class \code{"integer"} of length 1. Number of template reagents.}
    \item{\code{NQ}:}{Object of class \code{"integer"} of length 1. Number of query reagents.}
    \item{\code{C}:}{Object of class \code{"integer"} of length 1. Number of readout channels.}
    \item{\code{S}:}{Object of class \code{"integer"} of length 1. Number of screens.}
    \item{\code{F}:}{Object of class \code{"integer"} of length 1. Number of measurements or single experiments per screen.}
    \item{\code{reagents}:}{Object of class \code{"data.frame"} describing each reagents. Obligatory columns: \code{RID} and \code{TID}.}
    \item{\code{targets}:}{Object of class \code{"data.frame"} describing each target gene. Obligatory columns: \code{TID}, \code{Symbol}, \code{group}, \code{GID}.}
    \item{\code{templateDesign}:}{Object of class \code{"data.frame"} with \code{sgi@NT} rows  describing the template design. Obligatory columns: \code{TemplatePlate}, \code{Well}, \code{RID}, \code{QueryNr}. }
    \item{\code{queryDesign}:}{Object of class \code{"data.frame"} with \code{sgi@NQ} rows describing the query design. Obligatory columns: \code{Plate}, \code{TemplatePlate}, \code{QueryNr}, \code{RID}. }
    \item{\code{transformation}:}{Object of class \code{"character"} of length \code{sgi@C}. The transformation applied to the input data. }
    \item{\code{mainTemplate}:}{Object of class \code{"array"} with dimension \code{sgi@NT x sgi@S x sgi@C}. The main effect of the template reagents.}
    \item{\code{mainQuery}:}{Object of class \code{"array"} with dimension \code{sgi@NQ x sgi@S x sgi@C}. The main effect of the query reagents.}
    \item{\code{mainSderrTemplate}:}{Object of class \code{"array"} with dimension \code{sgi@NT x sgi@S x sgi@C}. The standard error of the main effect of the template reagents.}
    \item{\code{mainSderrQuery}:}{Object of class \code{"array"} with dimension \code{sgi@NQ x sgi@S x sgi@C}. The standard error of the main effect of the query reagents.}
    \item{\code{mainSdTemplate}:}{Object of class \code{"array"} with dimension \code{sgi@NQ x sgi@S x sgi@C}. The standard deviation of the main effect of the query reagents.}
    \item{\code{mainSdQuery}:}{Object of class \code{"array"} with dimension \code{sgi@NQ x sgi@S x sgi@C}. The standard deviation of the main effect of the query reagents.}
    \item{\code{mainTimeEffect}:}{Object of class \code{"array"} with dimension \code{sgi@NQ x sgi@S x sgi@C}. The systematic changes of the query main effects, e.g. decreasing cell number over time.}
    \item{\code{mainSpatialEffect}:}{Object of class \code{"array"} with dimension \code{sgi@F x sgi@S x sgi@C}. The systematic spatial plate effects.}
    \item{\code{mainSpatialEffectRow}:}{Object of class \code{"array"}. Spatial effects per row (as computed by Bscore).}
    \item{\code{mainSpatialEffectCol}:}{Object of class \code{"array"}. Spatial effects per column (as computed by Bscore).}
    \item{\code{mainNeg}:}{Object of class \code{"array"} with dimension \code{sgi@S x sgi@C}. The main effect of the negative control.}
    \item{\code{mainNegTemplate}:}{Object of class \code{"array"} with dimension \code{sgi@S x sgi@C}. The template main effect of the negative control.}
    \item{\code{mainNegQuery}:}{Object of class \code{"array"} with dimension \code{sgi@S x sgi@C}. The query main effect of the negative control.}
    \item{\code{data2mainTemplate}:}{Object of class \code{"integer"} with dimension \code{sgi@F}. Mapping of single experiments to template reagents.}
    \item{\code{data2mainQuery}:}{Object of class \code{"integer"} with dimension \code{sgi@F}. Mapping of single experiments to query reagents.}
    \item{\code{ni.model}:}{Object of class \code{"array"} with dimension \code{sgi@F x sgi@S x sgi@C}. The expected values of the non-interacting model.}
    \item{\code{pi}:}{Object of class \code{"array"} with dimension \code{sgi@F x sgi@S x sgi@C}. The pairwise interaction score.}
    \item{\code{plateeffect}:}{Object of class \code{"array"}.}
    \item{\code{p.value}:}{Object of class \code{"array"} with dimension \code{sgi@NT x sgi@NQ x sgi@S x sgi@C} describing the p.values.}
    \item{\code{q.value}:}{Object of class \code{"array"} with dimension \code{sgi@NT x sgi@NQ x sgi@S x sgi@C} describing the q.values.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "RNAinteract")}: ... }
	 }
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
showClass("RNAinteract")
}
\keyword{classes}
