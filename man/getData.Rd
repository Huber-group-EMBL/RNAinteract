\name{getData}
\alias{getData}
\title{
Primary access function for all screen data.
}
\description{
This function is the primary access function for a wide range of data from the screen. It does perform normalization, transformation, and reshaping if specified.
}
\usage{
getData(sgi, type = "data", format = "plain",
	design = "template", mixTemplateQuery = TRUE,
	screen = NULL, channel = NULL,
	do.trafo = TRUE, do.inv.trafo = FALSE,
	normalized = FALSE, withoutgroups = c(),
	drop = TRUE)
}
\arguments{
  \item{sgi}{An object of class \code{\link{RNAinteract}}}
  \item{type}{Specifies which data is returned. Possible values are:
\itemize{
\item "p.value", "q.value": returns the p-value or q-value as computed by \code{\link{computePValues}}.
\item "data": returns the input data.
\item "pi": returns the pairwise interaction score.		
\item "plateeffect": returns the plate effect estimated by \code{\link{normalizePlateEffect}}.
\item "ni.model" returns the non-interacting model as estimated by \code{\link{estimateMainEffect}}.
\item "main": returns the main effects.
\item "mainsderr": returns the std error of the main effects.
\item "mainsd": returns the std deviation of the main effects.
\item "maintime": returns the estimated time effect as estimated by \code{\link{normalizeMainEffectQuery}}
\item "mainspatial": returns the estimated spatial effect as estimated by \code{\link{normalizePlateEffect}}
}
}
  \item{format}{The output format. Possible values:
\itemize{
\item "plain": The data can be returned as a plain vector
\item "platelist": a list of plate matrices that can be passed to \code{\link{plotScreen}}
\item "reagentMatrix": All values for the same reagent pair are summarized in a matrix of dimension reagents x reagents
\item "targetMatrix": All values for the same gene pair are summarized in a matrix of dimension genes x genes
}}
  \item{design}{If type is one of the main effect types, the design can be specified to state if the "template" or "query" main effect is returned. }
  \item{mixTemplateQuery}{If TRUE, The template-query and query-template entries in the matrix are symmetrized.}
  \item{screen}{The screen names of which data should be returned.}
  \item{channel}{The channel names of which data should be returned.}
  \item{do.trafo}{Only effective, if type is "data". If TRUE, the data is transformed.}
  \item{do.inv.trafo}{Not effective if type is "data", "p.value", or "q.value". If TRUE, the values are back-transformed to the original scale.}
  \item{normalized}{If TRUE, the normalization data is returned.}
  \item{withoutgroups}{The genes from the specified groups are not returned in the data.}
  \item{drop}{If FALSE, the returned array is reduced in dimensions, whenever there is a dimension 1.}
}
\value{
An array containing the specified values is returned. In the case, the format is chosen to be "platelist", a list of matrices is returned.
}
\author{
Bernd Fischer
}
\seealso{
\code{\link{RNAinteract-package}}
}
\examples{
data("sgi")

# get the original data, as plain file, reshaped in plate layout,
# reshaped and summarized as target matrix
D <- getData(sgi, type="data", do.inv.trafo = TRUE)
Dplatelayout <- getData(sgi, type="data", 
	     format="platelist", do.inv.trafo = TRUE)
splots::plotScreen(Dplatelayout[["1"]][["nrCells"]],
	     nx=sgi@pdim[2], ny=sgi@pdim[1], ncol=3)
Dmatrix <- getData(sgi, type="data", 
	     format="targetMatrix", do.inv.trafo = TRUE)

# get main effects as plate layout with specified transformation
# (usually log-transformed)
Mplatelayout <- getData(sgi, type="main", design="template", 
	     screen="1", channel="nrCells", format="platelist")
splots::plotScreen(Mplatelayout, nx=sgi@pdim[2], ny=sgi@pdim[1],
	     ncol=3)

# get non-interacting model and pairwise interaction scores as matrix
NImatrix <- getData(sgi, type="ni.model", format="targetMatrix")
PImatrix <- getData(sgi, type="pi", format="targetMatrix")
PIplatelayout <- getData(sgi, type="main", design="query",
	     screen="1", channel="nrCells", format="platelist")
splots::plotScreen(PIplatelayout, nx=sgi@pdim[2], ny=sgi@pdim[1],
             ncol=3)

# get p-values and q-values
PVmatrix <- getData(sgi, type="p.value", format="targetMatrix")
QVmatrix <- getData(sgi, type="q.value", format="targetMatrix")
}
\keyword{ attribute }

