
plotScreenData <- function(sgi, screen, channel, type="data", design = "template", do.trafo = TRUE, do.inv.trafo = FALSE, normalized = FALSE, plotScreen.args=list(ncol=6L, do.legend=TRUE, fill = c("red","white","blue"))) {
  stopifnot( is( sgi, "RNAinteract" ))

  D <- getData(sgi, type = type, design = design, do.trafo = do.trafo, do.inv.trafo = do.inv.trafo, format = "platelist", screen = screen, channel = channel, normalized=normalized)
  legend.label <- attr(D, "axislabel")[channel]
  
  zrange <- range(unlist(D),finite=TRUE)
  if (type %in% c("pi", "maindiff")) {
    m <- 0
  } else {
    if (type == "data") {
      m <- sgi@mainNeg[screen,channel]
    } else {
      m <- median(unlist(D),na.rm=TRUE)
    }
  }
  r <- max(abs(zrange - m))
  zrange <- c(m-r,m+r)

  if (!plotScreen.args$do.legend) { legend.label <- FALSE }
  main = sprintf("screen %s (%s)", screen, channel)
  plotScreen(D, zrange = zrange, ncol = plotScreen.args$ncol, main = main, do.legend = TRUE, legend.label = legend.label, fill=plotScreen.args$fill, nx=sgi@pdim[2], ny = sgi@pdim[1])
}

plotReplicateScatter <- function(sgi, screen, channel, type="data", design = "template", do.trafo = TRUE, do.inv.trafo = FALSE, normalized = FALSE) {
  stopifnot( is( sgi, "RNAinteract" ) )

  res <- getReplicateData(sgi, screen, channel, type, design, do.trafo, do.inv.trafo, normalized = normalized)
  main = sprintf("within-screen-replicates screen %s (%s)", screen, channel)
  smoothScatter(res$x, res$y, pch=20, main=main, xlab=res$lab,ylab=res$lab)
}

plotIndDesignScatter <- function(sgi, screen, channel, type="data", design = "template", do.trafo = TRUE, do.inv.trafo = FALSE, normalized = FALSE) {
  stopifnot( is( sgi, "RNAinteract" ) )

  res <- getIndDesignData(sgi, screen, channel, type, design, do.trafo, do.inv.trafo, normalized = normalized)

  main = sprintf("ind. reagents, screen %s (%s)", screen, channel)
  smoothScatter(res$x,res$y,pch=20,main=main, xlab=res$lab,ylab=res$lab)
}

reportScreenData <- function(sgi, type="data", design = "template", do.trafo = TRUE, do.inv.trafo = FALSE, verbose = 0, path = ".", dir = "screenplots", prefix = "screenplot", png.args = list(width = 960, height=960), pdf.args = list(width=7, height=7), plotScreen.args=list(ncol=6L, do.legend=TRUE, fill = c("red","white","blue")), png.scatter.args = list(width=400, height=400), pdf.scatter.args = list(width=7, height=7), report=NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (verbose > 0) {
    time1 <- Sys.time()
    message("Plot screen. \r", appendLF = FALSE)
  }
  dir.create(sprintf("%s/%s", path, dir), showWarnings = FALSE)

  if (type == "data") {
    prefix2 = "data"
  } else {
    if (type %in% c("pi", "ni.model")) {
      prefix2 <- type
    } else {
      prefix2 <- paste(type,"-",design,collapse = "",sep="")
    }
  }

  if (!is.null(report)) {
    hwrite(sprintf("<b>%s-%s</b>",prefix,prefix2), br=TRUE, page=report)
    for (c in getChannelNames(sgi)) {
      hwrite(paste(c," ",collapse=""), br=FALSE, page=report)
      for (s in getScreenNames(sgi)) {
        hwrite(sprintf("[%s]",s), link = sprintf("%s/index-%s-%s-%s-%s.html", dir, prefix, prefix2, s, c), target = "showframe", br=FALSE, page=report)
      }
      hwrite("<br>", br=FALSE, page=report)
    }
    hwrite("<br>", br=FALSE, page=report)
  }

  zz <- 0
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      if (verbose > 1) {
        zz <- zz + 1
        ZZ <- (sgi@S*sgi@C)
        message(sprintf("Plot screen. Nr %d from %d.\r", zz, ZZ), appendLF = FALSE)
      }
      fn.png = sprintf("%s/%s/%s-%s-%s-%s.png", path, dir, prefix, prefix2, s, c)
      fn.pdf = sprintf("%s/%s/%s-%s-%s-%s.pdf", path, dir, prefix, prefix2, s, c)
      png(width=png.args$width, height=png.args$height, filename=fn.png)
      plotScreenData(sgi, type = type, design = design, do.trafo = do.trafo, do.inv.trafo = do.inv.trafo, screen = s, channel = c,normalized=TRUE, plotScreen.args = plotScreen.args)
      dev.off()
      pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
      plotScreenData(sgi, type = type, design = design, do.trafo = do.trafo, do.inv.trafo = do.inv.trafo, screen = s, channel = c,normalized=TRUE, plotScreen.args = plotScreen.args)
      dev.off()

      fn.scatter.rep.png = sprintf("%s/%s/%s-repscatter-%s-%s-%s.png", path, dir, prefix, prefix2, s, c)
      fn.scatter.rep.pdf = sprintf("%s/%s/%s-repscatter-%s-%s-%s.pdf", path, dir, prefix, prefix2, s, c)
      png(width=png.scatter.args$width, height=png.scatter.args$height, filename=fn.scatter.rep.png)
      plotReplicateScatter(sgi, type = type, design = design, do.trafo = do.trafo, do.inv.trafo = do.inv.trafo, screen = s, channel = c,normalized=TRUE)
      dev.off()
      pdf(width=pdf.scatter.args$width, height=pdf.scatter.args$height, file=fn.scatter.rep.pdf)
      plotReplicateScatter(sgi, type = type, design = design, do.trafo = do.trafo, do.inv.trafo = do.inv.trafo, screen = s, channel = c,normalized=TRUE)
      dev.off()

      fn.scatter.ind.png = sprintf("%s/%s/%s-indscatter-%s-%s-%s.png", path, dir, prefix, prefix2, s, c)
      fn.scatter.ind.pdf = sprintf("%s/%s/%s-indscatter-%s-%s-%s.pdf", path, dir, prefix, prefix2, s, c)
      png(width=png.scatter.args$width, height=png.scatter.args$height, filename=fn.scatter.ind.png)
      plotIndDesignScatter(sgi, type = type, design = design, do.trafo = do.trafo, do.inv.trafo = do.inv.trafo, screen = s, channel = c,normalized=TRUE)
      dev.off()
      pdf(width=pdf.scatter.args$width, height=pdf.scatter.args$height, file=fn.scatter.ind.pdf)
      plotIndDesignScatter(sgi, type = type, design = design, do.trafo = do.trafo, do.inv.trafo = do.inv.trafo, screen = s, channel = c,normalized=TRUE)
      dev.off()

      p = openPage(sprintf("%s/%s/index-%s-%s-%s-%s.html",path, dir, prefix, prefix2, s, c))

      fn.png = sprintf("%s-%s-%s-%s.png", prefix, prefix2, s, c)
      fn.pdf = sprintf("%s-%s-%s-%s.pdf", prefix, prefix2, s, c)
      fn.scatter.rep.png = sprintf("%s-repscatter-%s-%s-%s.png", prefix, prefix2, s, c)
      fn.scatter.rep.pdf = sprintf("%s-repscatter-%s-%s-%s.pdf", prefix, prefix2, s, c)
      fn.scatter.ind.png = sprintf("%s-indscatter-%s-%s-%s.png", prefix, prefix2, s, c)
      fn.scatter.ind.pdf = sprintf("%s-indscatter-%s-%s-%s.pdf", prefix, prefix2, s, c)
      hwrite("<h3>Screen Plot</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<h4>(%s; screen=%s; channel=%s)</h4>",prefix2, s, c),table=FALSE,br=FALSE,page=p)
      hwriteImage(fn.png,link=fn.pdf, cellpadding=3, cellspacing=0,border=1,table=FALSE,br=TRUE,page=p)
      hwrite("<h3>With-in-Screen Replicates</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<h4>(%s; screen=%s; channel=%s)</h4>",prefix2, s, c),table=FALSE,br=FALSE,page=p)
      hwrite("For each reagent pair there is one dot in the scatter plot. It compares two measurements of the same reagent pair within the same screen.<br>",table=FALSE,br=TRUE,page=p)
      hwriteImage(fn.scatter.rep.png,link=fn.scatter.rep.pdf, cellpadding=3, cellspacing=0,border=1,table=FALSE,br=TRUE,page=p)
      hwrite("<h3>Independent Reagents</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<h4>(%s; screen=%s; channel=%s)</h4>",prefix2, s, c),table=FALSE,br=FALSE,page=p)
      hwrite("For each target pair there is at most one dot in the scatter plot. It compares measurements of independent reagents targeting the same targets.<br>",table=FALSE,br=TRUE,page=p)
      hwriteImage(fn.scatter.ind.png,link=fn.scatter.ind.pdf, cellpadding=3, cellspacing=0,border=1,table=FALSE,br=TRUE,page=p)

      closePage(p,splash=FALSE)
    }
  }
  if (verbose > 0) {
    time2 = Sys.time()
    message(sprintf("Plot screen. Done. Time: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }
  invisible(NULL)
}




