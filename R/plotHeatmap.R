

plotHeatmap <- function(sgi, screen, channel, pi.max=NULL, main=expression(paste(pi,"-score")), hc.row = NULL, hc.col = NULL, withoutgroups=c("neg", "pos")) {
  stopifnot( is( sgi, "RNAinteract" ) )

  PI <- getData(sgi, type="pi", format="targetMatrix",screen=screen,channel=channel, withoutgroups=withoutgroups)
  
  grid.sgiHeatmap(PI)
}

