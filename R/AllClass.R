
############################################
### RNAinteract
############################################

## validityRNAinteract <- function(object) {
##   return(TRUE)
## }

setClass("RNAinteract",
   representation(data = "array",
                  screenNames = "character",
                  channelNames = "character",
                  well = "character",
                  plate = "integer",
                  pdim = "integer",
                  NT   = "integer",
                  NQ   = "integer",
                  C    = "integer",
                  S    = "integer",
                  F    = "integer",
                  reagents          = "data.frame",
                  targets           = "data.frame",
                  templateDesign    = "data.frame",
                  queryDesign       = "data.frame",
                  transformation    = "character",
                  mainTemplate      = "array",
                  mainQuery         = "array",
                  mainSderrTemplate = "array",
                  mainSderrQuery    = "array",
                  mainSdTemplate    = "array",
                  mainSdQuery       = "array",
                  mainTimeEffect    = "array",
                  mainSpatialEffect = "array",
                  mainSpatialEffectRow = "array",
                  mainSpatialEffectCol = "array",
                  mainNeg           = "array",
                  mainNegTemplate   = "array",
                  mainNegQuery      = "array",
                  data2mainTemplate = "integer",
                  data2mainQuery    = "integer",
                  ni.model          = "array",  # values of non-interacting model
                  pi                = "array",  # pairwiese interaction
                  plateeffect       = "array",
                  p.value           = "array",
                  q.value           = "array")
)

setMethod("show",signature="RNAinteract", function(object) {
  message("RNA interaction screen")
  message("Nr of template reagents:     ", object@NT)
  message("Nr of query reagents:        ", object@NQ)
  message("Nr of experiments in screen: ", object@F)
  message("Nr of channels:              ", object@C)
  message("Nr of screens:               ", object@S)
})

