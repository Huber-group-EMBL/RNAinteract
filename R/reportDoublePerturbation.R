
reportDoublePerturbation <- function(sgi, verbose = 0, path = ".", dir = "doublePerturbations",
                                    prefix = "doublePerturbationPlots", report = NULL,
                                    withoutgroups=c("neg", "pos"),
                                    png.args = list(width = 500, height=500),
                                    pdf.args = list(width=7, height=7), ...) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (verbose > 0) {
    time1 <- Sys.time()
    message("Plot double perturbations. \r", appendLF=FALSE)
  }
  dir.create(sprintf("%s/%s", path, dir), showWarnings = FALSE)

  if (!is.null(report)) {
    hwrite(sprintf("<b>%s</b>",prefix), br=TRUE, page=report)
    for (c in getChannelNames(sgi)) {
      hwrite(paste(c," ",collapse=""), br=FALSE, page=report)
      for (s in getScreenNames(sgi)) {
        hwrite(sprintf("[%s]",s), link = sprintf("%s/index-%s-%s-%s.html",dir, prefix, s, c), target="showframe", br=FALSE, page=report)
      }
      hwrite("<br>", br=FALSE, page=report)
    }
    hwrite("<br>", br=FALSE, page=report)
  }

  symbols <- colnames(getData(sgi, screen = getScreenNames(sgi)[1], channel = getChannelNames(sgi)[1], type="data", format="targetMatrix", mixTemplateQuery = FALSE, normalized = TRUE, withoutgroups=withoutgroups))
  zz <- 0
  ZZ <- sgi@S * sgi@C * length(symbols)
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      D <- getData(sgi, screen = s, channel = c, type="data", format="targetMatrix", mixTemplateQuery = FALSE, normalized = TRUE, withoutgroups=withoutgroups)
      MainNeg <- getMainNeg(sgi, screen = s, channel = c)
      D <- D - MainNeg
      MT <- getData(sgi, screen = s, channel = c, type="main", format="targetMatrix", mixTemplateQuery = FALSE, normalized = TRUE, design="template", withoutgroups=withoutgroups)
      MQ <- getData(sgi, screen = s, channel = c, type="main", format="targetMatrix", mixTemplateQuery = FALSE, normalized = TRUE, design="query", withoutgroups=withoutgroups)
      QV <- getData(sgi, screen = s, channel = c,type="q.value", format="targetMatrix", mixTemplateQuery = FALSE, withoutgroups=withoutgroups)
      PV <- getData(sgi, screen = s, channel = c,type="p.value", format="targetMatrix", mixTemplateQuery = FALSE, withoutgroups=withoutgroups)
      PI <- getData(sgi, screen = s, channel = c,type="pi", format="targetMatrix", mixTemplateQuery = FALSE, 
                        withoutgroups=withoutgroups)

      scale = getScale(sgi,channel=c)

      for (tid in symbols) {
        if (verbose > 1) {
          zz <- zz + 1
          message(sprintf("Plot double perturbation. Nr %d from %d.\r", zz, ZZ), appendLF = FALSE)
        }

        itg <- plotDoublePerturbation(sgi=sgi, screen=s, channel=c, target=tid,
                                     withoutgroups=withoutgroups,
                                     D=D, MT=MT, MQ=MQ, PV=PV, QV=QV, PI=PI, draw=FALSE, ...)

        fn.png = sprintf("%s/%s/%s-%s-%s-%s.png", path, dir, prefix, s, c, tid)
        fn.pdf = sprintf("%s/%s/%s-%s-%s-%s.pdf", path, dir, prefix, s, c, tid)
        png(width=png.args$width, height=png.args$height, file=fn.png)
        grid.draw(itg)
        dev.off()
        pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
        grid.draw(itg)
        dev.off()
      }

      I <- order(-MQ[1,])
      fn.png = sprintf("%s-%s-%s-%s.png", prefix, s, c, symbols[I])
      fn.pdf = sprintf("%s-%s-%s-%s.pdf", prefix, s, c, symbols[I])

      p = openPage(sprintf("%s/%s/index-%s-%s-%s.html",path, dir, prefix, s, c))

      hwrite("<h3>Double Perturbation Plots</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
      hwrite("Each panel shows all interactions of one target. On the x-axis the single mutant level of the second target plotted. The y-axis shows the double mutant level. The two solid lines show the single mutant level of the negative control (e.g. wild type). The dotted gray lines show the single mutant levels. The horizonatal gray lines is the single mutant level of the first target (The name is written in the title of the panel). The diagonal gray line is the single mutant level of the second target. The orange line is the sum of the two gray lines. This is the level of the non-interaction model. The target pairs on the diagonal gray line are the non-interacting target pairs. The target pairs with a q-value < 5% (Thus the interactions) are marked with a text label.",table=FALSE,br=TRUE,page=p)
      Img <- hwriteImage(matrix(fn.png,nr=length(fn.png),nc=1),table=FALSE)
      row.names(Img) <- symbols[I]
      hwrite(Img,link=fn.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

      closePage(p,splash=FALSE)
    }
  }
  if (verbose > 0) {
    time2 = Sys.time()
    message(sprintf("Plot double perturbations. Done. Time: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }
  invisible(NULL)
}





