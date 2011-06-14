
computePValues <- function (sgi, method = "pooled.ttest", mixTemplateQuery = TRUE, verbose = 0) {
  stopifnot( is( sgi, "RNAinteract" ) )
  if (verbose > 0) {
    time1 <- Sys.time()
    message("Estimate genetic interactions. \r", appendLF = FALSE)
  }

  templateTargets <- dimnames(sgi@p.value)[[1]]
  queryTargets <- dimnames(sgi@p.value)[[2]]
  ID <- array(1:prod(dim(sgi@p.value)[1:2]),dim = dim(sgi@p.value)[1:2],dimnames = dimnames(sgi@p.value)[1:2])
  if (mixTemplateQuery) {
    join <- queryTargets[which(queryTargets %in% templateTargets)]
    if (length(join) > 0) {
      m1 = match(join,queryTargets)
      ID2 = ID[join,join]
      ID2[lower.tri(ID2)] <- t(ID2)[lower.tri(ID2)]
      ID[join,join] = ID2
      ID[1:length(ID)] = match(ID[1:length(ID)], unique(ID[1:length(ID)]) )
    }
  }
  m1 <- match(sgi@targets[sgi@reagents[sgi@templateDesign$RID[sgi@data2mainTemplate],"TID"], "Symbol"], dimnames(ID)[[1]])
  m2 <- match(sgi@targets[sgi@reagents[sgi@queryDesign$RID[sgi@data2mainQuery], "TID"], "Symbol"], dimnames(ID)[[2]])
  id <- ID[(m2-1)*dim(ID)[1]+m1]

  PI <- getData(sgi, type="pi", format = "plain",drop=FALSE)
  maxid = max(ID)
  maxn = max(table(id))

  if (method %in% c("HotellingT2")) {
    p.value = rep(NA,maxid)
    names(p.value) = 1:maxid
    for (s in 1:sgi@S) {
      for (i in 1:maxid) {
        I = which(id == i)
        n = length(I)
        if (n >= sgi@C) {
          X = matrix(NA,nr=n,nc=sgi@C)
          z=0
          for (c in 1:sgi@C) {
            X[,c] = PI[I,s,c]
          }
          try( {
            p.value[i] = HotellingsT2(X)$p.value
          }, silent = TRUE )
        }
      }
      for (c in 1:sgi@C) {
        sgi@p.value[,,s,c] <- p.value[ID]
      }
      I = which(is.finite(p.value))
      qjob = qvalue(p.value[I])
      q.value = rep(NA,length(p.value))
      q.value[I] = qjob$qvalues
      for (c in 1:sgi@C) {
        sgi@q.value[,,s,c] <- q.value[ID]
      }
    }
  } else {
    for (s in 1:sgi@S) {
      for (c in 1:sgi@C) {
        X = matrix(NA,nr=maxid,nc=maxn)
        pi = PI[,s,c]
        id2 = id
        for (i in 1:maxn) {
          m = match(1:maxid,id2)
          X[,i] = pi[m]
          id2[m] = NA
        }
        n <- apply(is.finite(X),1,sum)
        if (method %in% c("pooled.ttest","ttest")) {
          m <- apply(X,1,mean,na.rm=TRUE)
          stddev <- apply(X,1,sd,na.rm=TRUE)
          if (method == "pooled.ttest") {
            S <- median(stddev,na.rm=TRUE)
            stddev[stddev < S] <- S
          }
          t <- sqrt(n) * m / stddev
          p.value <- 2*pt(-abs(t),df=n-1)
        }
        if (method == "limma") {
          fit <- lmFit(X)
          fit <- eBayes(fit)
          p.value <- fit$p.value
        }
        sgi@p.value[,,s,c] <- p.value[ID]
        I = which(is.finite(p.value))
        qjob = qvalue(p.value[I])
        q.value = rep(NA,length(p.value))
        q.value[I] = qjob$qvalues
        sgi@q.value[,,s,c] <- q.value[ID]
      }
    }
  }
  if (verbose > 0) {
    time2 = Sys.time()
    message(sprintf("Estimate genetic interactions. Done. Time: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }
  return(sgi)
}

