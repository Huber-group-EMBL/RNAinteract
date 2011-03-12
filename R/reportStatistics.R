
zPrime <- function(x.mean,x.sd, y.mean, y.sd) {
  zprime = 1 - (3*x.sd + 3*y.sd) / abs(x.mean - y.mean)
  return(zprime)
}

reportStatistics <- function(sgi, verbose = 0, path = ".", dir = "stats", prefix = "stat", report = NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (verbose > 0) {
    time1 <- Sys.time()
    message("Write statistics. \r", appendLF = FALSE)
  }
  dir.create(sprintf("%s/%s", path, dir), showWarnings = FALSE)

  if (!is.null(report)) {
    hwrite(sprintf("<b>statistics</b>",prefix), br=TRUE, page=report)
    hwrite("[Z'-score]", link = sprintf("%s/index-zprime-%s.html",dir, prefix), target = "showframe", br=TRUE, page=report)
    hwrite("[Correlation]", link = sprintf("%s/index-cor-%s.html",dir, prefix), target = "showframe", br=TRUE, page=report)
    hwrite("[nrInteractions]", link = sprintf("%s/index-nr-%s.html",dir, prefix), target = "showframe", br=TRUE, page=report)
    hwrite("<br>", br=FALSE, page=report)
  }

  ##########################
  ## Z-prime score
  ##########################
  controls = sgi@targets$Symbol[sgi@targets$group %in% c("pos","neg")]

  neg <- sgi@targets$Symbol[sgi@targets$group == "neg"]
  pos <- sgi@targets$Symbol[sgi@targets$group == "pos"]
  rown <- matrix("", nr=sgi@S * sgi@C * length(neg), nc = 3)
  colnames(rown) <- c("screen", "channel", "neg")
  fn.data.png <- matrix("", nr=sgi@S * sgi@C * length(neg), nc = length(pos))
  fn.data.pdf <- matrix("", nr=sgi@S * sgi@C * length(neg), nc = length(pos))
  fn.density.png <- matrix("", nr=sgi@S * sgi@C * length(neg), nc = length(pos))
  fn.density.pdf <- matrix("", nr=sgi@S * sgi@C * length(neg), nc = length(pos))
  Z <- matrix(NA, nr=sgi@S * sgi@C * length(neg), nc = length(pos))
  colnames(Z) <- pos
  colnames(fn.data.png) <- pos
  colnames(fn.data.pdf) <- pos
  colnames(fn.density.png) <- pos
  colnames(fn.density.pdf) <- pos
  zz <- 0
  Mean = getMain(sgi, type = "main", summary="target",drop=FALSE)
  SD = getMain(sgi, type = "mainsd", summary="target",drop=FALSE)
  D = getData(sgi, normalized=TRUE,drop=FALSE)
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      I = which(sgi@templateDesign$RID[sgi@data2mainTemplate] %in% c(neg,pos))
      r = range(D[I,,c], finite=TRUE)
      r[2] = r[2] + 0.15*(r[2] - r[1])
      for (n in neg) {
        zz <- zz + 1
        for (p in pos) {
          Z[zz,p] <- zPrime(Mean[n,s,c], SD[n,s,c], Mean[p,s,c], SD[p,s,c])

          I = which(sgi@templateDesign$RID[sgi@data2mainTemplate] == p)
          J = which(sgi@templateDesign$RID[sgi@data2mainTemplate] == n)
          fn.data.png[zz,p] = sprintf("%s_data_%s_%s_%s_%s.png",prefix,n,p,s,c)
          fn.data.pdf[zz,p] = sprintf("%s_data_%s_%s_%s_%s.png",prefix,n,p,s,c)
          png(file=sprintf("%s/%s/%s_data_%s_%s_%s_%s.png",path,dir,prefix,n,p,s,c))
          plot(D[J,s,c],ylim=r,col="red",pch=20, main=sprintf("screen %s (%s)", s, c), xlab="plate", ylab=attr(D,"axislab")[c]);
          points(D[I,s,c],col="blue",pch=20)
          legend("topright",c(sprintf("pos. control: %s",p),sprintf("neg. control: %s",n)), fill=c("blue","red"),inset=0.02)
          dev.off()
          pdf(file=sprintf("%s/%s/%s_data_%s_%s_%s_%s.pdf",path,dir,prefix,n,p,s,c))
          plot(D[J,s,c],ylim=r,col="red",pch=20, main=sprintf("screen %s (%s)", s, c), xlab="plate", ylab=attr(D,"axislab")[c]);
          points(D[I,s,c],col="blue",pch=20)
          legend("topright",c(sprintf("pos. control: %s",p),sprintf("neg. control: %s",n)), fill=c("blue","red"),inset=0.02)
          dev.off()
          fn.density.png[zz,p] = sprintf("%s_density_%s_%s_%s_%s.png",prefix,n,p,s,c)
          fn.density.pdf[zz,p] = sprintf("%s_density_%s_%s_%s_%s.pdf",prefix,n,p,s,c)
          png(file=sprintf("%s/%s/%s_density_%s_%s_%s_%s.png",path,dir,prefix,n,p,s,c))
          multidensity(list(pos=D[J,s,c], neg = D[I,s,c]),
                       xlab=attr(D,"axislab")[c], main = sprintf("screen %s (%s)", s, c),lwd=3,
                       legend = list(x="topright",legend = c(sprintf("pos. control: %s",p),sprintf("neg. control: %s",n)), fill=c("blue","red"),inset=0.02))
          dev.off()
          pdf(file=sprintf("%s/%s/%s_density_%s_%s_%s_%s.pdf",path,dir,prefix,n,p,s,c))
          multidensity(list(pos=D[J,s,c], neg = D[I,s,c]),
                       xlab=attr(D,"axislab")[c], main = sprintf("screen %s (%s)", s, c),lwd=3,
                       legend = list(x="topright",legend = c(sprintf("pos. control: %s",p),sprintf("neg. control: %s",n)), fill=c("blue","red"),inset=0.02))
          dev.off()
        }
        rown[zz,1] <- s
        rown[zz,2] <- c
        rown[zz,3] <- n
      }
    }
  }

  p = openPage(sprintf("%s/%s/index-zprime-%s.html",path, dir, prefix))

  write.table(cbind(rown,Z), file = sprintf("%s/%s/stat-zprime-%s.txt",path, dir, prefix), row.names = FALSE, sep="\t", quote=FALSE)

  Z.text <- Z
  Z.text[] <- sprintf("%0.3f",Z)

  hwrite("<h3>Z-prime score for controls</h3>",table=FALSE,br=FALSE,page=p)
  hwrite(sprintf("Download as text-file: <a href=stat-zprime-%s.txt>[Zprime]</a><br>",prefix),table=FALSE,br=FALSE,page=p)
  hwrite(cbind(rown, Z.text),style='text-align:right',col.bgcolor=c(rep("#CCCCCC",3),rep("#EEEEEE",length(pos))), cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

  hwrite("<h3>Controls</h3>",table=FALSE,br=FALSE,page=p)
  Img.data <- hwriteImage(fn.data.png, table=FALSE)
  colnames(Img.data) <- pos
  hwrite(cbind(rown, Img.data),link=cbind(matrix("",nr=nrow(fn.data.pdf),nc=3),fn.data.pdf),style='text-align:right',col.bgcolor=c(rep("#CCCCCC",3),rep("#EEEEEE",length(pos))), cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

  hwrite("<h3>Density of Controls</h3>",table=FALSE,br=FALSE,page=p)
  Img.density <- hwriteImage(fn.density.png, table=FALSE)
  colnames(Img.density)
  hwrite(cbind(rown, Img.density),link=cbind(matrix("",nr=nrow(fn.density.pdf),nc=3),fn.density.pdf), style='text-align:right',col.bgcolor=c(rep("#CCCCCC",3),rep("#EEEEEE",length(pos))), cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

  closePage(p,splash=FALSE)

  ##########################
  ## correlations
  ##########################
  p = openPage(sprintf("%s/%s/index-cor-%s.html",path, dir, prefix))

  C1 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  fn.png.1 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  fn.pdf.1 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  row.names(C1) <- row.names(fn.png.1) <- row.names(fn.pdf.1) <- getChannelNames(sgi)
  colnames(C1)  <- colnames(fn.png.1)  <- colnames(fn.pdf.1)  <- getScreenNames(sgi)
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      res <- getReplicateData(sgi, screen=s, channel=c, type="data", normalized = TRUE)
      C1[c,s] <- cor(res$x,res$y,use="complete.obs")
      main = sprintf("within-screen-replicates screen %s (%s)", s, c)
      png(width=400,height=400, file=sprintf("%s/%s/%s-repscatter-data-%s-%s.png", path, dir, prefix, s, c))
      smoothScatter(res$x, res$y, pch=20, main=main, xlab=res$lab,ylab=res$lab)
      dev.off()
      pdf(file=sprintf("%s/%s/%s-repscatter-data-%s-%s.pdf", path, dir, prefix, s, c))
      smoothScatter(res$x, res$y, pch=20, main=main, xlab=res$lab,ylab=res$lab)
      dev.off()
      fn.png.1[c,s] <- sprintf("%s-repscatter-data-%s-%s.png", prefix, s, c)
      fn.pdf.1[c,s] <- sprintf("%s-repscatter-data-%s-%s.pdf", prefix, s, c)
    }
  }

  C2 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  fn.png.2 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  fn.pdf.2 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  row.names(C2) <- row.names(fn.png.2) <- row.names(fn.pdf.2) <- getChannelNames(sgi)
  colnames(C2)  <- colnames(fn.png.2)  <- colnames(fn.pdf.2)  <- getScreenNames(sgi)
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      res <- getIndDesignData(sgi, screen=s, channel=c, type="data", normalized = TRUE)
      C2[c,s] <- cor(res$x,res$y,use="complete.obs")      
      main = sprintf("ind. reagents, screen %s (%s)", s, c)
      png(width=400,height=400, file=sprintf("%s/%s/%s-indscatter-data-%s-%s.png", path, dir, prefix, s, c))
      smoothScatter(res$x, res$y, pch=20, main=main, xlab=res$lab,ylab=res$lab)
      dev.off()
      pdf(file=sprintf("%s/%s/%s-indscatter-data-%s-%s.pdf", path, dir, prefix, s, c))
      smoothScatter(res$x, res$y, pch=20, main=main, xlab=res$lab,ylab=res$lab)
      dev.off()
      fn.png.2[c,s] <- sprintf("%s-indscatter-data-%s-%s.png", prefix, s, c)
      fn.pdf.2[c,s] <- sprintf("%s-indscatter-data-%s-%s.pdf", prefix, s, c)
    }
  }

  if (sgi@S > 1) {
    C3 <- matrix(NA, nr=sgi@C, nc=sgi@S * (sgi@S-1) / 2)
    fn.png.3 <- matrix(NA, nr=sgi@C, nc=sgi@S * (sgi@S-1) / 2)
    fn.pdf.3 <- matrix(NA, nr=sgi@C, nc=sgi@S * (sgi@S-1) / 2)
    row.names(C3) <- row.names(fn.png.3) <- row.names(fn.pdf.3) <- getChannelNames(sgi)
    colnames(C3)  <- colnames(fn.png.3)  <- colnames(fn.pdf.3) <- rep("", ncol(C3))
    D <- getData(sgi, normalized = TRUE,drop=FALSE)
    for (c in getChannelNames(sgi)) {
      zz <- 0
      for (i1 in 1:(sgi@S-1)) {
        for (i2 in (i1+1):sgi@S) {
          s1 <- getScreenNames(sgi)[i1]
          s2 <- getScreenNames(sgi)[i2]
          zz <- zz + 1
          C3[c,zz] <- cor(D[,s1,c],D[,s2,c],use="complete.obs")
          colnames(C3)[zz] = sprintf("%s - %s", s1, s2)
          main = sprintf("between screen replicate, screen %s - %s (%s)", s1, s2, c)
          png(width=400,height=400, file=sprintf("%s/%s/%s-betweenscatter-data-%s-%s-%s.png", path, dir, prefix, s1, s2, c))
          smoothScatter(D[,s1,c],D[,s2,c], pch=20, main=main, xlab=res$lab,ylab=res$lab)
          dev.off()
          pdf(file=sprintf("%s/%s/%s-betweenscatter-data-%s-%s-%s.pdf", path, dir, prefix, s1, s2, c))
          smoothScatter(D[,s1,c],D[,s2,c], pch=20, main=main, xlab=res$lab,ylab=res$lab)
          dev.off()
          fn.png.3[c,zz] <- sprintf("%s-betweenscatter-data-%s-%s-%s.png", prefix, s1, s2, c)
          fn.pdf.3[c,zz] <- sprintf("%s-betweenscatter-data-%s-%s-%s.pdf", prefix, s1, s2, c)
        }
      }
    }
  }
  C1.text <- C1
  C1.text[] <- sprintf("%0.3f", C1)
  C2.text <- C2
  C2.text[] <- sprintf("%0.3f", C2)
  if (sgi@S > 1) {
    C3.text <- C3
    C3.text[] <- sprintf("%0.3f", C3)
  }
  M1 <- hwrite(C1.text, cellpadding=3, cellspacing=0,border=1)
  M2 <- hwrite(C2.text, cellpadding=3, cellspacing=0,border=1)
  if (sgi@S > 1) {
    M3 <- hwrite(C3.text, cellpadding=3, cellspacing=0,border=1)
  }

  hwrite("<h3>Correlation of readout input data</h3>",table=FALSE,br=FALSE,page=p)
  if (sgi@S > 1) {
    C <- cbind(C1,C2,C3)
  } else {
    C <- cbind(C1,C2)
  }
  write.table(C, file = sprintf("%s/%s/stat-cor-input-%s.txt",path, dir, prefix), row.names = FALSE, sep="\t", quote=FALSE)
  hwrite(sprintf("Download as text-file: <a href=stat-cor-input-%s.txt>[CorInput]</a><br>",prefix),table=FALSE,br=FALSE,page=p)
  if (sgi@S > 1) {
    M <- matrix(c(M1, M2, M3),nr=1, nc=3)
  } else {
    M <- matrix(c(M1, M2,""),nr=1, nc=3)
  }
  colnames(M) <- c("within-screen tech. rep.", "within-screen ind. designs", "between screen")
  hwrite(M, cellpadding=3, cellspacing=0, br=TRUE, border=1, page=p)

  hwrite("<b>Within screen (technical) replicates</b>",table=FALSE,br=TRUE,page=p)
  Img <- hwriteImage(fn.png.1,table=FALSE)
  row.names(Img) <- row.names(fn.png.1)
  colnames(Img) <- colnames(fn.png.1)
  hwrite(Img, link=fn.pdf.1, cellpadding=3, cellspacing=0, br=TRUE,border=1, page=p)
  
  hwrite("<b>Within screen independent designs</b>",table=FALSE,br=TRUE,page=p)
  Img <- hwriteImage(fn.png.2,table=FALSE)
  row.names(Img) <- row.names(fn.png.2)
  colnames(Img) <- colnames(fn.png.2)
  hwrite(Img, link=fn.pdf.2, cellpadding=3, cellspacing=0,br=TRUE, border=1, page=p)

  if (sgi@S > 1) {
    hwrite("<b>Between screen (biological) replicates</b>",table=FALSE,br=TRUE,page=p)
    Img <- hwriteImage(fn.png.3,table=FALSE)
    row.names(Img) <- row.names(fn.png.3)
    colnames(Img) <- colnames(fn.png.3)
    hwrite(Img, link=fn.pdf.3, cellpadding=3, cellspacing=0,br=TRUE, border=1, page=p)
  }

  C1 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  fn.png.1 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  fn.pdf.1 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  row.names(C1) <- row.names(fn.png.1) <- row.names(fn.pdf.1) <- getChannelNames(sgi)
  colnames(C1)  <- colnames(fn.png.1)  <- colnames(fn.pdf.1)  <- getScreenNames(sgi)
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      res <- getReplicateData(sgi, screen=s, channel=c, type="pi")
      C1[c,s] <- cor(res$x,res$y,use="complete.obs")
      main = sprintf("within-screen-replicates screen %s (%s)", s, c)
      png(width=400,height=400, file=sprintf("%s/%s/%s-repscatter-pi-%s-%s.png", path, dir, prefix, s, c))
      smoothScatter(res$x, res$y, pch=20, main=main, xlab=res$lab,ylab=res$lab)
      dev.off()
      pdf(file=sprintf("%s/%s/%s-repscatter-pi-%s-%s.pdf", path, dir, prefix, s, c))
      smoothScatter(res$x, res$y, pch=20, main=main, xlab=res$lab,ylab=res$lab)
      dev.off()
      fn.png.1[c,s] <- sprintf("%s-repscatter-pi-%s-%s.png", prefix, s, c)
      fn.pdf.1[c,s] <- sprintf("%s-repscatter-pi-%s-%s.pdf", prefix, s, c)
    }
  }

  C2 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  fn.png.2 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  fn.pdf.2 <- matrix(NA, nr=sgi@C, nc=sgi@S)
  row.names(C2) <- row.names(fn.png.2) <- row.names(fn.pdf.2) <- getChannelNames(sgi)
  colnames(C2)  <- colnames(fn.png.2)  <- colnames(fn.pdf.2)  <- getScreenNames(sgi)
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      res <- getIndDesignData(sgi, screen=s, channel=c, type="pi")
      C2[c,s] <- cor(res$x,res$y,use="complete.obs")      
      main = sprintf("ind. reagents, screen %s (%s)", s, c)
      png(width=400,height=400, file=sprintf("%s/%s/%s-indscatter-pi-%s-%s.png", path, dir, prefix, s, c))
      smoothScatter(res$x, res$y, pch=20, main=main, xlab=res$lab,ylab=res$lab)
      dev.off()
      pdf(file=sprintf("%s/%s/%s-indscatter-pi-%s-%s.pdf", path, dir, prefix, s, c))
      smoothScatter(res$x, res$y, pch=20, main=main, xlab=res$lab,ylab=res$lab)
      dev.off()
      fn.png.2[c,s] <- sprintf("%s-indscatter-pi-%s-%s.png", prefix, s, c)
      fn.pdf.2[c,s] <- sprintf("%s-indscatter-pi-%s-%s.pdf", prefix, s, c)
    }
  }

  C3 <- matrix(NA, nr=sgi@C, nc=sgi@S * (sgi@S-1) / 2)
  fn.png.3 <- matrix(NA, nr=sgi@C, nc=sgi@S * (sgi@S-1) / 2)
  fn.pdf.3 <- matrix(NA, nr=sgi@C, nc=sgi@S * (sgi@S-1) / 2)
  row.names(C3) <- row.names(fn.png.3) <- row.names(fn.pdf.3) <- getChannelNames(sgi)
  colnames(C3)  <- colnames(fn.png.3)  <- colnames(fn.pdf.3) <- rep("", ncol(C3))
  D <- getData(sgi, type="pi", format="plain",drop=FALSE)
  if (sgi@S > 1) {
    for (c in getChannelNames(sgi)) {
      zz <- 0
      for (i1 in 1:(sgi@S-1)) {
        for (i2 in (i1+1):sgi@S) {
          s1 <- getScreenNames(sgi)[i1]
          s2 <- getScreenNames(sgi)[i2]
          zz <- zz + 1
          C3[c,zz] <- cor(D[,s1,c],D[,s2,c],use="complete.obs")
          colnames(C3)[zz] = sprintf("%s - %s", s1, s2)
          main = sprintf("between screen replicate, screen %s - %s (%s)", s1, s2, c)
          png(width=400,height=400, file=sprintf("%s/%s/%s-betweenscatter-pi-%s-%s-%s.png", path, dir, prefix, s1, s2, c))
          smoothScatter(D[,s1,c],D[,s2,c], pch=20, main=main, xlab=res$lab,ylab=res$lab)
          dev.off()
          pdf(file=sprintf("%s/%s/%s-betweenscatter-pi-%s-%s-%s.pdf", path, dir, prefix, s1, s2, c))
          smoothScatter(D[,s1,c],D[,s2,c], pch=20, main=main, xlab=res$lab,ylab=res$lab)
          dev.off()
          fn.png.3[c,zz] <- sprintf("%s-betweenscatter-pi-%s-%s-%s.png", prefix, s1, s2, c)
          fn.pdf.3[c,zz] <- sprintf("%s-betweenscatter-pi-%s-%s-%s.pdf", prefix, s1, s2, c)
        }
      }
    }
  }

  C1.text <- C1
  C1.text[] <- sprintf("%0.3f", C1)
  C2.text <- C2
  C2.text[] <- sprintf("%0.3f", C2)
  if (sgi@S > 1) {
    C3.text <- C3
    C3.text[] <- sprintf("%0.3f", C3)
  }
  M1 <- hwrite(C1.text, cellpadding=3, cellspacing=0,border=1)
  M2 <- hwrite(C2.text, cellpadding=3, cellspacing=0,border=1)
  if (sgi@S > 1) {
    M3 <- hwrite(C3.text, cellpadding=3, cellspacing=0,border=1)
  }

  hwrite("<h3>Correlation of pairwise interactions</h3>",table=FALSE,br=FALSE,page=p)
  C <- cbind(C1,C2,C3)
  write.table(C, file = sprintf("%s/%s/stat-cor-pi-%s.txt",path, dir, prefix), row.names = FALSE, sep="\t", quote=FALSE)
  hwrite(sprintf("Download as text-file: <a href=stat-cor-pi-%s.txt>[CorPI]</a><br>",prefix),table=FALSE,br=FALSE,page=p)
  if (sgi@S > 1) {
    M <- matrix(c(M1, M2, M3),nr=1, nc=3)
  } else {
    M <- matrix(c(M1, M2,""),nr=1, nc=3)
  }
  colnames(M) <- c("within-screen tech. rep.", "within-screen ind. designs", "between screen")
  hwrite(M, cellpadding=3, cellspacing=0,border=1, br=TRUE, page=p)

  hwrite("<b>Within screen (technical) replicates</b>",table=FALSE,br=TRUE,page=p)
  Img <- hwriteImage(fn.png.1,table=FALSE)
  row.names(Img) <- row.names(fn.png.1)
  colnames(Img) <- colnames(fn.png.1)
  hwrite(Img, link=fn.pdf.1, cellpadding=3, cellspacing=0, br=TRUE,border=1, page=p)
  
  hwrite("<b>Within screen independent designs</b>",table=FALSE,br=TRUE,page=p)
  Img <- hwriteImage(fn.png.2,table=FALSE)
  row.names(Img) <- row.names(fn.png.2)
  colnames(Img) <- colnames(fn.png.2)
  hwrite(Img, link=fn.pdf.2, cellpadding=3, cellspacing=0,br=TRUE, border=1, page=p)

  if (sgi@S > 1) {
    hwrite("<b>Between screen (biological) replicates</b>",table=FALSE,br=TRUE,page=p)
    Img <- hwriteImage(fn.png.3,table=FALSE)
    row.names(Img) <- row.names(fn.png.3)
    colnames(Img) <- colnames(fn.png.3)
    hwrite(Img, link=fn.pdf.3, cellpadding=3, cellspacing=0,br=TRUE, border=1, page=p)
  }
  closePage(p,splash=FALSE)

  ##########################
  ## nr of interactions
  ##########################

  N <- matrix(NA, nr=sgi@C, nc=sgi@S)
  row.names(N) <- getChannelNames(sgi)
  colnames(N) <- getScreenNames(sgi)
  q.value <- getData(sgi, type = "q.value", format="targetMatrix",drop=FALSE)
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      N[c,s] = sum((q.value[,,s,c])[upper.tri(q.value[,,s,c])] <= 0.05, na.rm=TRUE)
    }
  }

  p = openPage(sprintf("%s/%s/index-nr-%s.html",path, dir, prefix))

  write.table(N, file = sprintf("%s/%s/%s-nr.txt",path, dir, prefix), sep="\t", quote=FALSE)

  hwrite("<h3>Number of Interactions</h3>",table=FALSE,br=FALSE,page=p)
  hwrite("The number of interactions with a q-value smaller or equal 5%<br>",table=FALSE,br=FALSE,page=p)
  hwrite(sprintf("Download as text-file: <a href=%s-nr.txt>[NrInteractions]</a><br>",prefix),table=FALSE,br=FALSE,page=p)
  hwrite(N,style='text-align:right', cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

  closePage(p,splash=FALSE)

  if (verbose > 0) {
    time2 = Sys.time()
    message(sprintf("Write statistics. Done. Time: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }
  invisible(NULL)
}



