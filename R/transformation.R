
transform <- function(sgi, x, channel) {
  fct = getTransformation(sgi, channel = channel)
  x = fct(x)
  return(x)
}

invtransform <- function(sgi, x, channel) {
  fct = getInvTransformation(sgi, channel = channel)
  x = fct(x)
  return(x)
}

getTransformation <- function(sgi, channel) {
  stopifnot( is( sgi, "RNAinteract" ) )

  fct = switch (sgi@transformation[channel],
    "log"   = function(x) { return(log(x)) },
    "log2"  = function(x) { return(log2(x)) },
    "log10" = function(x) { return(log10(x)) },
    "sqrt"  = function(x) { return(sqrt(x)) },
    "sqr"   = function(x) { return(x^2) },
    function(x) { return(x) } )

  return(fct)
}

getInvTransformation <- function(sgi, channel) {
  stopifnot( is( sgi, "RNAinteract" ) )

  fct = switch (sgi@transformation[channel],
    "log"   = function(x) { return(exp(x)) },
    "log2"  = function(x) { return(2^x) },
    "log10" = function(x) { return(10^x) },
    "sqrt"  = function(x) { return(x^2) },
    "sqr"   = function(x) { return(sqrt(x)) },
    function(x) { return(x) } )

  return(fct)
}

getScale <- function(sgi, channel) {
  stopifnot( is( sgi, "RNAinteract" ) )

  fct = switch (sgi@transformation[channel],
    "log"   = "log-scale",
    "log2"  = "log2-scale",
    "log10" = "log10-scale",
    "sqrt"  = "sqrt-scale",
    "sqr"   = "square-scale",
    "linear-scale" )

  return(fct)
}
