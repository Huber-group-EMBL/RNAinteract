\name{RNAinteract-package}

\alias{RNAinteract-package}
\alias{RNAinteract}

\docType{package}

\title{Analysis of Pairwise Interaction Screens.}

\description{
The package contains functions to organize the data from (single- and multi-parametric) genetic interaction screens. Methods to estimate main effects (single perturbation effects) and pairwise interactions. p-values are computed. Furthermore a comprehensive html-report is generated.
}

\section{Package content}{

  Class \code{RNAinteract} (Documentation: \code{\link{RNAinteract-class}})

  Data input and creating of an object of class \code{RNAinteract}.
  \itemize{
    \item \code{\link{createCellHTSFromFiles}}
    \item \code{\link{createRNAinteract}}, \code{\link{createRNAinteractFromFiles}}
  }
  
  Data access
  \itemize{
	\item \code{\link{getData}} Primary data access function for multiple types of screen data.
	\item \code{\link{getMain}}, \code{\link{getMainNeg}} access to main effects.
	\item \code{\link{getReplicateData}}, \code{\link{getIndDesignData}} Comparing replicate data.
	\item \code{\link{getChannelNames}}, \code{\link{getScreenNames}}, \code{\link{getScale}} 
  }

  Subsetting, summarizing, and binding screens
  \itemize{
	\item \code{\link{sgisubset}}, \code{\link{sgisubsetQueryDesign}}
	\item \code{\link{bindscreens}}
	\item \code{\link{summarizeScreens}}
  }

  Main effects and pairwise interactions
  \itemize{
	\item \code{\link{estimateMainEffect}} 
	\item \code{\link{normalizeMainEffectQuery}}, \code{\link{normalizeMainEffectTemplate}}, \code{\link{normalizePlateEffect}}
	\item \code{\link{computePI}}, \code{\link{computePValues}}
	\item \code{\link{embedPCA}}
  }

  Plotting
  \itemize{
	\item \code{\link{plotDoublePerturbation}}, \code{\link{plotHeatmap}} standard plot functions
	\item \code{\link{doublePerturbationGrob}}, \code{\link{grid.doublePerturbation}}, \code{\link{grid.sgiHeatmap}} specialized grid plotting functions for experts
  }

  HTML report
  \itemize{
	\item \code{\link{startReport}}, \code{\link{endReport}} starting and finalizing a report
	\item \code{\link{reportAnnotation}}, \code{\link{reportStatistics}} global reports
	\item \code{\link{reportDoublePerturbation}}, \code{\link{reportGeneLists}}, \code{\link{reportHeatmap}}, \code{\link{reportMainEffects}}, \code{\link{reportNetworks}}, \code{\link{reportScreenData}} reports specific for each screen and each channel
  }
}

\details{
See vignette("RNAinteract") for details.
}

\author{
Bernd Fischer

Maintainer: Bernd Fischer <bernd.fischer@embl.de>
}
\references{
T. Horn, T. Sandmann, B. Fischer, W. Huber, M. Boutros. Mapping of Signalling Networks through Synthetic 
Genetic Interaction Analysis by RNAi. Nature Methods, 2011.
}
\keyword{ package }


