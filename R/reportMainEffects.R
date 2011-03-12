
plotScatterErrorbars <- function(main1, main2, mad1, mad2, controls=NULL, controls.col=rainbow(length(controls)),main="", same.range = TRUE, ...) {
  N = length(main1)
  m1 = mad1
  m2 = mad2
  if (is.null(mad1)) { m1 = 0 }
  if (is.null(mad2)) { m2 = 0 }
  if (same.range) {
    range = c(0,0)
    range[1] = min(main1-m1,main2-m2,na.rm=TRUE)
    range[2] = max(main1+m1,main2+m2,na.rm=TRUE)
    xrange = yrange = range
  } else {
    xrange = c(0,0)
    xrange[1] = min(main1-m1,na.rm=TRUE)
    xrange[2] = max(main1+m1,na.rm=TRUE)
    yrange = c(0,0)
    yrange[1] = min(main2-m2,na.rm=TRUE)
    yrange[2] = max(main2+m2,na.rm=TRUE)
  }
  if (same.range) {
    plot(range,range,xlim=xrange,ylim=yrange,type="l",col="gray", main=main, ...)
    points(main1,main2,pch=20,xlim=range,ylim=range, main=main, ...)
  } else {
    plot(main1,main2,pch=20,xlim=xrange,ylim=yrange,type="p", main=main, ...)
  }
  for (i in 1:N) {
    if (!is.null(mad1)) { lines(c(main1[i]-mad1[i],main1[i]+mad1[i]),c(main2[i],main2[i]),type="l") }
    if (!is.null(mad2)) { lines(c(main1[i],main1[i]),c(main2[i]-mad2[i],main2[i]+mad2[i]),type="l") }
  }
  if (!is.null(controls)) {
    I <- match(names(main1),controls)
    J <- which(!is.na(I))
    I <- I[J]
    if (length(J) > 0) {
      points(main1[J],main2[J],pch=20,col=controls.col[I])
      for (i in 1:length(J)) {
        if (!is.null(mad1)) { lines(c(main1[J[i]]-mad1[J[i]],main1[J[i]]+mad1[J[i]]),c(main2[J[i]],main2[J[i]]),type="l",col=controls.col[I[i]]) }
        if (!is.null(mad2)) { lines(c(main1[J[i]],main1[J[i]]),c(main2[J[i]]-mad2[J[i]],main2[J[i]]+mad2[J[i]]),type="l",col=controls.col[I[i]]) }
      }
      I <- sort(unique(I))
      legend("topleft", controls[I], fill = controls.col[I],inset=0.02)
    }
  }
  invisible(NULL)
}

plotMainEffects <- function(sgi, screen, channel, design="template", TemplatePlate=1, QueryNr=c(1,2), plot.args=list(), do.inv.trafo=FALSE, type="main", errorbar = "mainsderr", compare.targets=FALSE,controls=NULL, controls.col=rainbow(length(controls)), same.range = TRUE) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (length(screen) == 1)        screen <- rep(screen,2)
  if (length(channel) == 1)       channel <- rep(channel,2)
  if (length(design) == 1)        design <- rep(design,2)
  if (length(TemplatePlate) == 1) TemplatePlate <- rep(TemplatePlate,2)
  if (length(QueryNr) == 1)       QueryNr <- rep(QueryNr,2)
  if (length(type) == 1)          type <- rep(type,2)
  if (length(errorbar) == 1)      errorbar <- rep(errorbar,2)
  if (length(do.inv.trafo) == 1)  do.inv.trafo <- rep(do.inv.trafo,2)

  main <- list()
  sderr <- list(c(NULL),c(NULL))
  lab <- list()
  for (i in 1:2) {
    M <- switch(design[i],
                template = getMain(sgi, screen=screen[i], channel=channel[i], type=type[i], design=design[i], QueryNr = QueryNr[i], do.inv.trafo = do.inv.trafo[i]),
                query = getMain(sgi, screen=screen[i], channel=channel[i], type=type[i], design=design[i], TemplatePlate = TemplatePlate[i], do.inv.trafo = do.inv.trafo[i]) )
    lab[[i]] <- attr(M,"axislabel")
    main[[i]] <- M
    if (do.inv.trafo[i]) {
      errorbar[i] = "None"
    }
    if (errorbar[i] != "None") {
      sderr[[i]] <- switch(design[i],
                           template = getMain(sgi, screen=screen[i], channel=channel[i], type=errorbar[i], design=design[i], QueryNr = QueryNr[i], do.inv.trafo = FALSE),
                           query = getMain(sgi, screen=screen[i], channel=channel[i], type=errorbar[i], design=design[i], TemplatePlate = TemplatePlate[i], do.inv.trafo = FALSE) )
    }
  }

  if (compare.targets) {
    TID1 <- sgi@targets[sgi@reagents[names(main[[1]]),"TID"],"Symbol"]
    TID2 <- sgi@targets[sgi@reagents[names(main[[2]]),"TID"],"Symbol"]
    RID1 <- names(main[[1]])
    RID2 <- names(main[[2]])
    A <-  c()
    B <-  c()
    for (i in 1:length(TID1)) {
      I <- which(TID2 == TID1[i])
      if (length(I) > 0) {
        I <- I[RID2[I] != RID1[i]]
        if (length(I) > 0) {
          A <- c(A,rep(i,length(I)))
          B <- c(B,I)
        }
      }
    }
    main[[1]] <- main[[1]][A]
    main[[2]] <- main[[2]][B]
    if (length(sderr[[1]]) > 0) { sderr[[1]] <- sderr[[1]][A] }
    if (length(sderr[[2]]) > 0) { sderr[[2]] <- sderr[[2]][B] }    
  } else {
    m <- match(names(main[[1]]),names(main[[2]]))
    main[[1]] <- main[[1]][!is.na(m)]
    main[[2]] <- main[[2]][m[!is.na(m)]]
    if (length(sderr[[1]]) > 0) { sderr[[1]] <- sderr[[1]][!is.na(m)] }
    if (length(sderr[[2]]) > 0) { sderr[[2]] <- sderr[[2]][m[!is.na(m)]] }
  }

  if (!compare.targets) {
    names(main[[1]]) = sgi@targets[sgi@reagents[names(main[[1]]),"TID"],"Symbol"]
    names(main[[2]]) = names(main[[1]])
  }
  plotScatterErrorbars(main[[1]],main[[2]],sderr[[1]],sderr[[2]], xlab=lab[[1]], ylab=lab[[2]], same.range=same.range, controls=controls, controls.col=controls.col,main="Main Effects")
}

reportMainEffects <- function(sgi, verbose = 0, path = ".", dir = "maineffects", prefix = "maineffects", png.args = list(width = 500, height=500), pdf.args = list(width=7, height=7), plot.args=list(), report = NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (verbose > 0) {
    time1 <- Sys.time()
    message("Plot main effects. \r", appendLF = FALSE)
  }

  dir.create(sprintf("%s/%s", path, dir), showWarnings = FALSE)

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

  controls = sgi@targets$Symbol[sgi@targets$group %in% c("pos","neg")]
  
  zz <- 0
  ZZ <- sgi@S * sgi@C
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      if (verbose > 1) {
        zz <- zz + 1
        message(sprintf("Plot main effects. Nr %d from %d.\r", zz, ZZ), appendLF = FALSE)
      }

      p = openPage(sprintf("%s/%s/index-%s-%s-%s.html",path, dir, prefix, s, c))

      IT <- sort(unique( sgi@queryDesign$TemplatePlate ))
      IQ <- sort(unique( sgi@templateDesign$QueryNr ))
      if (length(IQ) > 1) {
        fn2.png <- c()
        fn2.pdf <- c()
        QN = c()
        for (i in 1:(length(IQ)-1)) {
          for (j in (i+1):length(IQ)) {
            fn.png <- sprintf("%s/%s/%s-template-%s-%s-%d-%d.png", path, dir, prefix, s, c, IQ[i], IQ[j])
            fn.pdf <- sprintf("%s/%s/%s-template-%s-%s-%d-%d.pdf", path, dir, prefix, s, c, IQ[i], IQ[j])
            png(width=png.args$width, height=png.args$height, file=fn.png)
            plotMainEffects(sgi, screen=s, channel=c, design="template", QueryNr=c(IQ[i],IQ[j]), plot.args=plot.args, controls=controls)
            dev.off()
            pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
            plotMainEffects(sgi, screen=s, channel=c, design="template", QueryNr=c(IQ[i],IQ[j]), plot.args=plot.args, controls=controls)
            dev.off()

            fn2.png <- c(fn2.png,sprintf("%s-template-%s-%s-%d-%d.png", prefix, s, c, IQ[i], IQ[j]))
            fn2.pdf <- c(fn2.pdf,sprintf("%s-template-%s-%s-%d-%d.pdf", prefix, s, c, IQ[i], IQ[j]))
            QN <- c(QN, sprintf("%d - %d", IQ[i], IQ[j]))
          }
        }
        hwrite("<h3>Main Effects of Template Reagents</h3>",table=FALSE,br=FALSE,page=p)
        hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
        hwrite("The main effect of the template reagent of different copies on the template plates are compared to each other. If a normalization of the time effect is performed, then the normalized main effects are plotted.",table=FALSE,br=TRUE,page=p)
        Img <- hwriteImage(matrix(fn2.png,nr=length(fn2.png),nc=1),table=FALSE)
        row.names(Img) <- QN
        hwrite(Img,link=fn2.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)
      }

      if (length(IT) > 1) {
        fn2.png <- c()
        fn2.pdf <- c()
        QN = c()
        for (i in 1:(length(IT)-1)) {
          for (j in (i+1):length(IT)) {
            fn.png <- sprintf("%s/%s/%s-query-%s-%s-%d-%d.png", path, dir, prefix, s, c, IT[i], IT[j])
            fn.pdf <- sprintf("%s/%s/%s-query-%s-%s-%d-%d.pdf", path, dir, prefix, s, c, IT[i], IT[j])
            png(width=png.args$width, height=png.args$height, file=fn.png)
            plotMainEffects(sgi, screen=s, channel=c, design="query", TemplatePlate=c(IT[i],IT[j]), plot.args=plot.args, controls=controls)
            dev.off()
            pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
            plotMainEffects(sgi, screen=s, channel=c, design="query", TemplatePlate=c(IT[i],IT[j]), plot.args=plot.args, controls=controls)
            dev.off()

            fn2.png <- c(fn2.png,sprintf("%s-query-%s-%s-%d-%d.png", prefix, s, c, IT[i], IT[j]))
            fn2.pdf <- c(fn2.pdf,sprintf("%s-query-%s-%s-%d-%d.pdf", prefix, s, c, IT[i], IT[j]))
            QN <- c(QN, sprintf("%d - %d", IT[i], IT[j]))
          }
        }
        hwrite("<h3>Main Effects of Query Reagents</h3>",table=FALSE,br=FALSE,page=p)
        hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
        hwrite("The main effect of the query reagent on different template plates are compared to each other. If a normalization of the time effect is performed, then the normalized main effects are plotted.",table=FALSE,br=TRUE,page=p)
        Img <- hwriteImage(matrix(fn2.png,nr=length(fn2.png),nc=1),table=FALSE)
        row.names(Img) <- QN
        hwrite(Img,link=fn2.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)
      }

      fn2.png <- c()
      fn2.pdf <- c()
      QN = c()
      for (i in 1:length(IT)) {
        for (j in 1:length(IQ)) {
          fn.png <- sprintf("%s/%s/%s-templatequery-%s-%s-%d-%d.png", path, dir, prefix, s, c, IT[i], IQ[j])
          fn.pdf <- sprintf("%s/%s/%s-templatequery-%s-%s-%d-%d.pdf", path, dir, prefix, s, c, IT[i], IQ[j])
          png(width=png.args$width, height=png.args$height, file=fn.png)
          plotMainEffects(sgi, screen=s, channel=c, design=c("template","query"), TemplatePlate=IT[i],QueryNr=IQ[j], plot.args=plot.args, controls=controls)
          dev.off()
          pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
          plotMainEffects(sgi, screen=s, channel=c, design=c("template","query"), TemplatePlate=IT[i],QueryNr=IQ[j], plot.args=plot.args, controls=controls)
          dev.off()
          
          fn2.png <- c(fn2.png,sprintf("%s-templatequery-%s-%s-%d-%d.png", prefix, s, c, IT[i], IQ[j]))
          fn2.pdf <- c(fn2.pdf,sprintf("%s-templatequery-%s-%s-%d-%d.pdf", prefix, s, c, IT[i], IQ[j]))
          QN <- c(QN, sprintf("%d - %d", IT[i], IQ[j]))
        }
      }
      hwrite("<h3>Main Effects of Template compared to Query Reagents</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
      hwrite("The main effect of the template reagents are compared to the main effect of the query reagents. If a normalization of the time effect is performed, then the normalized main effects are plotted.",table=FALSE,br=TRUE,page=p)
      Img <- hwriteImage(matrix(fn2.png,nr=length(fn2.png),nc=1),table=FALSE)
      row.names(Img) <- QN
      hwrite(Img,link=fn2.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)
      
      fn2.png <- c()
      fn2.pdf <- c()
      QN = c()
      for (i in 1:length(IQ)) {
        fn.png <- sprintf("%s/%s/%s-inddesigntemplate-%s-%s-%d.png", path, dir, prefix, s, c, IQ[i])
        fn.pdf <- sprintf("%s/%s/%s-inddesigntemplate-%s-%s-%d.pdf", path, dir, prefix, s, c, IQ[i])
        png(width=png.args$width, height=png.args$height, file=fn.png)
        plotMainEffects(sgi, screen=s, channel=c, design="template", QueryNr=IQ[i], plot.args=plot.args, controls=controls, compare.targets=TRUE)
        dev.off()
        pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
        plotMainEffects(sgi, screen=s, channel=c, design="template", QueryNr=IQ[i], plot.args=plot.args, controls=controls, compare.targets=TRUE)
        dev.off()
        fn2.png <- c(fn2.png,sprintf("%s-inddesigntemplate-%s-%s-%d.png", prefix, s, c, IQ[i]))
        fn2.pdf <- c(fn2.pdf,sprintf("%s-inddesigntemplate-%s-%s-%d.pdf", prefix, s, c, IQ[i]))
        QN <- c(QN, sprintf("%d", IQ[i]))
      }
      hwrite("<h3>Template Main Effects of Independent Reagents</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
      hwrite("The main effect of the template reagents are compared to each other. Each dot compares different reagents targeting the same target. If a normalization of the time effect is performed, then the normalized main effects are plotted.",table=FALSE,br=TRUE,page=p)
      Img <- hwriteImage(matrix(fn2.png,nr=length(fn2.png),nc=1),table=FALSE)
      row.names(Img) <- QN
      hwrite(Img,link=fn2.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

      fn2.png <- c()
      fn2.pdf <- c()
      QN = c()
      for (i in 1:length(IT)) {
        fn.png <- sprintf("%s/%s/%s-inddesignquery-%s-%s-%d.png", path, dir, prefix, s, c, IT[i])
        fn.pdf <- sprintf("%s/%s/%s-inddesignquery-%s-%s-%d.pdf", path, dir, prefix, s, c, IT[i])
        png(width=png.args$width, height=png.args$height, file=fn.png)
        plotMainEffects(sgi, screen=s, channel=c, design="query", TemplatePlate=IT[i], plot.args=plot.args, controls=controls, compare.targets=TRUE)
        dev.off()
        pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
        plotMainEffects(sgi, screen=s, channel=c, design="query", TemplatePlate=IT[i], plot.args=plot.args, controls=controls, compare.targets=TRUE)
        dev.off()
        fn2.png <- c(fn2.png,sprintf("%s-inddesignquery-%s-%s-%d.png", prefix, s, c, IT[i]))
        fn2.pdf <- c(fn2.pdf,sprintf("%s-inddesignquery-%s-%s-%d.pdf", prefix, s, c, IT[i]))
        QN <- c(QN, sprintf("%d", IT[i]))
      }
      hwrite("<h3>Query Main Effects of Independent Reagents</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
      hwrite("The main effect of the query reagents are compared to each other. Each dot compares different reagents targeting the same target. If a normalization of the time effect is performed, then the normalized main effects are plotted.",table=FALSE,br=TRUE,page=p)
      Img <- hwriteImage(matrix(fn2.png,nr=length(fn2.png),nc=1),table=FALSE)
      row.names(Img) <- QN
      hwrite(Img,link=fn2.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

      fn.png <- sprintf("%s/%s/%s-timeeffect-%s-%s-%d.png", path, dir, prefix, s, c, IT[i])
      fn.pdf <- sprintf("%s/%s/%s-timeeffect-%s-%s-%d.pdf", path, dir, prefix, s, c, IT[i])
      png(width=png.args$width, height=png.args$height, file=fn.png)
      plot(sgi@mainQuery[,s,c],pch=20, main="Query Time Effect Normalization", xlab="time", ylab=sprintf("query main effect (ratio), unnormalized, (%s)", getScale(sgi, channel=c)))
      lines(sgi@mainTimeEffect[,s,c],col="red",lwd=3)
      dev.off()
      pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
      plot(sgi@mainQuery[,s,c],pch=20, main="Query Time Effect Normalization", xlab="time", ylab=sprintf("query main effect (ratio), unnormalized, (%s)", getScale(sgi, channel=c)))
      lines(sgi@mainTimeEffect[,s,c],col="red",lwd=3)
      dev.off()
      fn.png <- sprintf("%s-timeeffect-%s-%s-%d.png", prefix, s, c, IT[i])
      fn.pdf <- sprintf("%s-timeeffect-%s-%s-%d.pdf", prefix, s, c, IT[i])
      hwrite("<h3>Query Time Effect</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
      hwrite("The main effect of the query reagents are plotted in time order.  If a normalization of the time effect is performed, then the red line shows the estimated time effect.",table=FALSE,br=TRUE,page=p)
      Img <- hwriteImage(matrix(fn.png,nr=1,nc=1),table=FALSE)
##       row.names(Img) <- QN
      hwrite(Img,link=fn.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

      fn.png <- sprintf("%s/%s/%s-spatialeffectrow-%s-%s-%d.png", path, dir, prefix, s, c, IT[i])
      fn.pdf <- sprintf("%s/%s/%s-spatialeffectrow-%s-%s-%d.pdf", path, dir, prefix, s, c, IT[i])
      row = match(substr(sgi@templateDesign$Well,1,1),LETTERS)
      png(width=png.args$width, height=png.args$height, file=fn.png)
      plot(row,sgi@mainTemplate[,s,c],pch=20, main="Template Row Effect Normalization", xlab="row", ylab=sprintf("template main effect (ratio), unnormalized, (%s)", getScale(sgi, channel=c)))
      points(sgi@mainSpatialEffectRow[1,,s,c],col="red",pch=19)
      dev.off()
      pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
      plot(row,sgi@mainTemplate[,s,c],pch=20, main="Template Row Effect Normalization", xlab="row", ylab=sprintf("template main effect (ratio), unnormalized, (%s)", getScale(sgi, channel=c)))
      points(sgi@mainSpatialEffectRow[1,,s,c],col="red",pch=19)
      dev.off()
      fn.png <- sprintf("%s-spatialeffectrow-%s-%s-%d.png", prefix, s, c, IT[i])
      fn.pdf <- sprintf("%s-spatialeffectrow-%s-%s-%d.pdf", prefix, s, c, IT[i])
      hwrite("<h3>Template Row Effect</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
      hwrite("The main effect of the template reagents are plotted in row order.  If a normalization of the spatial effect is performed, then the red points show the estimated row effect.",table=FALSE,br=TRUE,page=p)
      Img <- hwriteImage(matrix(fn.png,nr=1,nc=1),table=FALSE)
##       row.names(Img) <- QN
      hwrite(Img,link=fn.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

      fn.png <- sprintf("%s/%s/%s-spatialeffectcol-%s-%s-%d.png", path, dir, prefix, s, c, IT[i])
      fn.pdf <- sprintf("%s/%s/%s-spatialeffectcol-%s-%s-%d.pdf", path, dir, prefix, s, c, IT[i])
      col = as.integer(substr(sgi@templateDesign$Well,2,3))
      png(width=png.args$width, height=png.args$height, file=fn.png)
      plot(col,sgi@mainTemplate[,s,c],pch=20, main="Template Column Effect Normalization", xlab="column", ylab=sprintf("template main effect (ratio), unnormalized, (%s)", getScale(sgi, channel=c)))
      points(sgi@mainSpatialEffectCol[1,,s,c],col="red",pch=19)
      dev.off()
      pdf(width=pdf.args$width, height=pdf.args$height, file=fn.pdf)
      plot(col,sgi@mainTemplate[,s,c],pch=20, main="Template Column Effect Normalization", xlab="column", ylab=sprintf("template main effect (ratio), unnormalized, (%s)", getScale(sgi, channel=c)))
      points(sgi@mainSpatialEffectCol[1,,s,c],col="red",pch=19)
      dev.off()
      fn.png <- sprintf("%s-spatialeffectcol-%s-%s-%d.png", prefix, s, c, IT[i])
      fn.pdf <- sprintf("%s-spatialeffectcol-%s-%s-%d.pdf", prefix, s, c, IT[i])
      hwrite("<h3>Template Column Effect</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
      hwrite("The main effect of the template reagents are plotted in column order.  If a normalization of the spatial effect is performed, then the red points show the estimated column effect.",table=FALSE,br=TRUE,page=p)
      Img <- hwriteImage(matrix(fn.png,nr=1,nc=1),table=FALSE)
##       row.names(Img) <- QN
      hwrite(Img,link=fn.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)
      
      closePage(p,splash=FALSE)
    }
  }
  if (verbose > 0) {
    time2 = Sys.time()
    message(sprintf("Plot main effects. Done. Time: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }
  invisible(NULL)
}



