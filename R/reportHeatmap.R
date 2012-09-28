
reportHeatmap <- function(sgi, verbose = 0, path = ".", dir = "heatmap", prefix = "heatmap", png.args = list(width = 1000, height = 1000), pdf.args = list(width = 15, height = 15), report = NULL, withoutgroups=c("neg", "pos")) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (verbose > 0) {
    time1 <- Sys.time()
    message("Write Heatmap. \r", appendLF = FALSE)
  }
  dir.create(sprintf("%s/%s", path, dir), showWarnings = FALSE)

  res <- list()
  if (!is.null(report)) {
    hwrite(sprintf("<b>%s</b>",prefix), br=TRUE, page=report)
    for (c in getChannelNames(sgi)) {
      hwrite(paste(c," ",collapse=""), br=FALSE, page=report)
      for (s in getScreenNames(sgi)) {
        hwrite(sprintf("[%s]",s), link = sprintf("%s/index-%s-%s-%s.html",dir, prefix, s, c), target = "showframe", br=FALSE, page=report)
      }
      hwrite("<br>", br=FALSE, page=report)
    }
    hwrite("<br>", br=FALSE, page=report)
  }

  zz <- 0
  ZZ <- sgi@S * sgi@C
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      if (verbose > 1) {
        zz <- zz + 1
        message(sprintf("Write Heatmap. Nr %d from %d.\r", zz, ZZ), appendLF = FALSE)
      }

      p = openPage(sprintf("%s/%s/index-%s-%s-%s.html",path, dir, prefix, s, c))
      
      PI <- getData(sgi, type="pi", format="targetMatrix",screen=s,channel=c, withoutgroups=withoutgroups)
      fn.png = sprintf("%s/%s/%s-%s-%s.png", path, dir, prefix, s, c)
      fn.pdf = sprintf("%s/%s/%s-%s-%s.pdf", path, dir, prefix, s, c)
      png(width=png.args$width,height=png.args$height, filename=fn.png)
      grid.sgiHeatmap(PI)
      dev.off()
      pdf(width=png.args$width,height=png.args$height, file=fn.pdf)
      grid.sgiHeatmap(PI)
      dev.off()
      fn.png = sprintf("%s-%s-%s.png", prefix, s, c)
      fn.pdf = sprintf("%s-%s-%s.pdf", prefix, s, c)
      hwrite("<h3>Heatmap</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
      Img <- hwriteImage(fn.png,table=FALSE)
      hwrite(Img,link=fn.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)
 
      closePage(p,splash=FALSE)
    }
  }

  if (verbose > 0) {
    time2 = Sys.time()
    message(sprintf("Write Heatmap. Done. Time elapsed: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }
  invisible(res)
}

