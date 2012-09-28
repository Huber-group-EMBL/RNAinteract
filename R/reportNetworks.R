
reportNetworks <- function(sgi, verbose = 0, path = ".", dir = "networks", prefix = "networks", Networks, qv = 0.05, withoutgroups = c("pos","neg"), report = NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (verbose > 0) {
    time1 <- Sys.time()
    message("Write networks. \r", appendLF = FALSE)
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

  q.value <- getData(sgi, type = "q.value", format="targetMatrix", withoutgroups=withoutgroups,drop=FALSE)
  NW <- array(0,dim=c(dim(q.value)[1],dim(q.value)[2],ncol(Networks)-2))
  dimnames(NW) <- list(dimnames(q.value)[[1]],dimnames(q.value)[[2]],colnames(Networks)[3:ncol(Networks)])
  m1 <- match(Networks$gene1, dimnames(q.value)[[1]])
  m2 <- match(Networks$gene2, dimnames(q.value)[[2]])
  for (i in 1:nrow(Networks)) {
    NW[m1[i],m2[i],] <- as.numeric(Networks[i,3:ncol(Networks)])
    NW[m2[i],m1[i],] <- as.numeric(Networks[i,3:ncol(Networks)])
  }

  res = list()
  zz <- 0
  ZZ <- sgi@S * sgi@C
  for (s in getScreenNames(sgi)) {
    res[[s]] = list()
    for (c in getChannelNames(sgi)) {
      res[[s]][[c]] = list()
      K <- dim(NW)[3]
      pv <- matrix(NA,nrow=K,ncol=2)
      colnames(pv) <- c("p-value","odds-ratio")
      row.names(pv) <- dimnames(NW)[[3]]
      fn.png <- matrix("",nrow=K,ncol=1)
      fn.pdf <- matrix("",nrow=K,ncol=1)
      C <- list()
      for (k in 1:K) {
        res[[s]][[c]][[k]] = list()
        name <- dimnames(NW)[[3]][k]
        D = data.frame(G=ifelse((q.value[,,s,c])[upper.tri(q.value[,,s,c])] <= qv,1,0), N=(NW[,,k])[upper.tri(NW[,,k])])
        D[[1]] = factor(D[[1]],levels=c(0,1))
        D[[2]] = factor(D[[2]],levels=c(0,1))
        colnames(D) = c(c,name)
        t = table(D)
        C[[k]] <- t
        ftest <- fisher.test(t,alternative="greater")
        pv[k,1] = ftest$p.value
        pv[k,2] = ftest$estimate
        res[[s]][[c]][[k]]$T = t
        res[[s]][[c]][[k]]$p.value = ftest$p.value
        res[[s]][[c]][[k]]$estimate = ftest$estimate

        D = data.frame(g=factor(c(0,0,1,1),levels=c(0,1)),n=factor(c(0,1,0,1),levels=c(0,1)),v=c(0,t[1,2],t[2,1],t[2,2]))
        fn.png[k,1] <- sprintf("%s/%s/%s-ballon_%s_%s_%s.png",path, dir, prefix,s,c,name)
        png(width=150,height=150,filename=fn.png[k,1])
        par(mar=c(0.1,0.1,0.1,0.1))
        balloonplot(D$n,D$g,D$v,xlab=name,ylab=c,main="",cum.margins=FALSE)
        dev.off()
        fn.pdf[k,1] <- sprintf("%s/%s/%s-ballon_%s_%s_%s.pdf",path, dir, prefix,s,c,name)
        pdf(width=3,height=3,file=fn.pdf[k,1])
        par(mar=c(0.1,0.1,0.1,0.1))
        balloonplot(D$n,D$g,D$v,xlab=name,ylab=c,main="",cum.margins=FALSE)
        dev.off()
        fn.png[k,1] <- sprintf("%s-ballon_%s_%s_%s.png",prefix,s,c,name)
        fn.pdf[k,1] <- sprintf("%s-ballon_%s_%s_%s.pdf",prefix,s,c,name)
      }

      p = openPage(sprintf("%s/%s/index-%s-%s-%s.html",path, dir, prefix, s, c))
      write.table(pv, file = sprintf("%s/%s/%s-pv-%s-%s.txt",path, dir, prefix, s,c), sep="\t", quote=FALSE)      
      hwrite("<h3>p-value</h3>",table=FALSE,br=FALSE,page=p)
      hwrite(sprintf("Download as text-file: <a href=%s-pv-%s-%s.txt>[p-value]</a><br>",prefix, s,c),table=FALSE,br=FALSE,page=p)
      hwrite(pv,style='text-align:right', cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)
      hwrite("<h3>Confusion matrices</h3>",table=FALSE,br=FALSE,page=p)
      Img <- hwriteImage(fn.png,table=FALSE)
      hwrite(Img, link=fn.pdf,style='text-align:right', cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)
      closePage(p,splash=FALSE)
    }
  }


  if (verbose > 0) {
    time2 = Sys.time()
    message(sprintf("Write networks. Done. Time: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }
  invisible(res)
}



