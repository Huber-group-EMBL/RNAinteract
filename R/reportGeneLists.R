
reportGeneLists <- function(sgi, verbose = 0, path = ".", dir = "hitlist", prefix = "hitlist", report = NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (verbose > 0) {
    time1 <- Sys.time()
    message("Write gene lists. \r", appendLF = FALSE)
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
        message(sprintf("Write gene lists. Nr %d from %d.\r", zz, ZZ), appendLF = FALSE)
      }

      q.value <- getData(sgi, screen = s, channel = c, type="q.value", format="targetMatrix")
      p.value <- getData(sgi, screen = s, channel = c, type="p.value", format="targetMatrix")
      gene1 <- matrix(row.names(p.value),nrow=dim(p.value)[1], ncol=dim(p.value)[2])
      gene2 <- matrix(colnames(p.value),nrow=dim(p.value)[2], ncol=dim(p.value)[1])
      gene2 <- t(gene2)
      ## maint <- getMain(sgi, design="template", screen=s, channel=c, summary="target")
      ## mainq <- getMain(sgi, design="query", screen=s, channel=c, summary="target")
      ## main = invtransform(sgi, 0.5 * (maint + mainq), channel=c)
      ## main1 <- matrix(main,nrow=dim(p.value)[1], ncol=dim(p.value)[2])
      ## main2 <- matrix(main,nrow=dim(p.value)[2], ncol=dim(p.value)[1])
      ## main2 <- t(main2)
      main1 <- getData(sgi, screen = s, channel = c, type="main", design="template", format="targetMatrix", do.inv.trafo=TRUE)
      main2 <- getData(sgi, screen = s, channel = c, type="main", design="query", format="targetMatrix", do.inv.trafo=TRUE)
      NI <- getData(sgi, screen = s, channel = c, type="ni.model", format="targetMatrix", do.inv.trafo=TRUE)
      D <- getData(sgi, screen = s, channel = c, type="data", format="targetMatrix", do.trafo=FALSE)
      PI <- getData(sgi, screen = s, channel = c, type="pi", format="targetMatrix", do.inv.trafo=TRUE)

      g <- c(gene2,gene1)
      m1 <- match(gene1,g)
      m2 <- match(gene2,g)
      m3 <- m1
      m4 <- m2
      m3[m2 < m1] <- m2[m2 < m1]
      m4[m2 < m1] <- m1[m2 < m1]
      m <- m3 * max(m4) + m4
      dup = !duplicated(m)
      dup[m3 == m4] = FALSE
      UT = matrix(dup,nrow=nrow(p.value), ncol=ncol(p.value))
      GL <- data.frame(gene1 = gene1[UT],gene2 = gene2[UT],q.value = q.value[UT], p.value = p.value[UT],
                       main1 = main1[UT], main2 = main2[UT],
                       neg = getMainNeg(sgi, screen=s, channel=c, do.inv.trafo=TRUE),
                       NI = NI[UT], Measured = D[UT], pi = PI[UT])

      ## main.log = 0.5 * (maint + mainq)
      ## main1.log <- matrix(main.log,nrow=dim(p.value)[1], ncol=dim(p.value)[2])
      ## main2.log <- matrix(main.log,nrow=dim(p.value)[2], ncol=dim(p.value)[1])
      ## main2.log <- t(main2.log)
      main1.log <- getData(sgi, screen = s, channel = c, type="main", design="template", format="targetMatrix", do.inv.trafo=FALSE)
      main2.log <- getData(sgi, screen = s, channel = c, type="main", design="query", format="targetMatrix", do.inv.trafo=FALSE)
      NI.log <- getData(sgi, screen = s, channel = c, type="ni.model", format="targetMatrix", do.inv.trafo=FALSE)
      D.log <- getData(sgi, screen = s, channel = c, type="data", format="targetMatrix", do.trafo=TRUE)
      PI.log <- getData(sgi, screen = s, channel = c, type="pi", format="targetMatrix", do.inv.trafo=FALSE)

      GL.log <- data.frame(gene1 = gene1[UT],gene2 = gene2[UT],q.value = q.value[UT], p.value = p.value[UT],
                       main1 = main1.log[UT], main2 = main2.log[UT],
                       neg = getMainNeg(sgi, screen=s, channel=c, do.inv.trafo=FALSE),
                       NI = NI.log[UT], Measured = D.log[UT], pi = PI.log[UT])

      ## isctrl1 <- matrix(names(maint) %in% controls,nrow=dim(p.value)[1], ncol=dim(p.value)[2])
      ## isctrl2 <- matrix(names(maint) %in% controls,nrow=dim(p.value)[2], ncol=dim(p.value)[1])
      ## isctrl2 <- t(isctrl2)
      ## isctrl <- isctrl1[UT] | isctrl2[UT]
      isctrl <- (GL$gene1 %in% controls) | (GL$gene2 %in% controls)
      
      qv <- GL$p.value
      qv[isctrl] <- qv[isctrl] + 2
      I <- order(qv)
      GL <- GL[I,]
      GL.log <- GL.log[I,]
      GL.text <- GL
      for (i in 3:10) {
        GL.text[[i]] <- sprintf("%0.3f",GL.text[[i]])
      }

      p = openPage(sprintf("%s/%s/index-%s-%s-%s.html",path, dir, prefix, s, c))

      write.table(GL, file = sprintf("%s/%s/%s-%s-%s.txt",path, dir, prefix, s, c), row.names = FALSE, sep="\t", quote=FALSE)
      write.table(GL.log, file = sprintf("%s/%s/%s-trafo-%s-%s.txt",path, dir, prefix, s, c), row.names = FALSE, sep="\t", quote=FALSE)
      write.table(PI.log, file = sprintf("%s/%s/%s-matrix-pi-%s-%s.txt",path, dir, prefix, s, c), row.names = TRUE, sep="\t", quote=FALSE)
      write.table(q.value, file = sprintf("%s/%s/%s-matrix-qvalue-%s-%s.txt",path, dir, prefix, s, c), row.names = TRUE, sep="\t", quote=FALSE)

      hwrite("<h3>Interaction Hit List</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<b>(screen=%s; channel=%s)</b>", s, c),table=FALSE,br=TRUE,page=p)
      hwrite("Download as text-file:",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<a href=%s-%s-%s.txt>[list]</a> ",prefix, s, c),table=FALSE,br=FALSE,page=p)
      if (getScale(sgi,c) != "linear-scale") {
        hwrite(sprintf("<a href=%s-trafo-%s-%s.txt>[list (%s)]</a> ",prefix, s, c, getScale(sgi,channel=c)),table=FALSE,br=FALSE,page=p)
      }
      hwrite(sprintf("<a href=%s-matrix-pi-%s-%s.txt>[matrix PI-score]</a> ",prefix, s, c),table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("<a href=%s-matrix-qvalue-%s-%s.txt>[matrix q-value]</a> ",prefix, s, c),table=FALSE,br=TRUE,page=p)
      row.names(GL.text) <- 1:nrow(GL.text)
      hwrite(GL.text,style='text-align:right',row.bgcolor=rep(c("#CCCCCC","#EEEEEE"),length.out=nrow(GL.text)), cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

      closePage(p,splash=FALSE)
    }
  }
  if (verbose > 0) {
    time2 = Sys.time()
    message(sprintf("Write gene lists. Done. Time elapsed: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }
  invisible(NULL)
}



