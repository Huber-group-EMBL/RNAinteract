
estimateMainEffect <- function(sgi, use.query = NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

  D <- getData(sgi, type="data", format="plain",drop=FALSE)
  I <- which((sgi@data2mainTemplate > 0) & (sgi@data2mainQuery > 0))
  D <- D[I,,,drop=FALSE]
  d2t <- sgi@data2mainTemplate[I]
  d2q <- sgi@data2mainQuery[I]
  ft <- factor(d2t, levels = 1:sgi@NT)
  fq <- factor(d2q, levels = 1:sgi@NQ)
  if (is.null(use.query)) {
    selected <- 1:length(ft)
  } else {
    selected <- which(fq %in% which(sgi@queryDesign$RID %in% use.query))
  }
  for (s in 1:sgi@S) {
    for (c in 1:sgi@C) {
      sgi@mainTemplate[,s,c] <- 0.5 * tapply(D[selected,s,c],ft[selected],median,na.rm=TRUE)
      sgi@mainQuery[,s,c]    <- 0.5 * tapply(D[,s,c],fq,median,na.rm=TRUE)
      for (z in 1:20) {
        sgi@mainTemplate[,s,c] <- tapply(D[selected,s,c] - sgi@mainQuery[d2q[selected],s,c],   ft[selected],median,na.rm=TRUE)
        sgi@mainQuery[,s,c]    <- tapply(D[,s,c] - sgi@mainTemplate[d2t,s,c],fq,median,na.rm=TRUE)
      }

      # Standard deviation and standard error of main effects
      sgi@mainSdTemplate[,s,c] <- tapply(D[selected,s,c] - sgi@mainTemplate[d2t[selected],s,c] - sgi@mainQuery[d2q[selected],s,c],ft[selected],mad,na.rm=TRUE)
      sgi@mainSdQuery[,s,c]    <- tapply(D[,s,c] - sgi@mainTemplate[d2t,s,c] - sgi@mainQuery[d2q,s,c],fq,mad,na.rm=TRUE)
      sgi@mainSderrTemplate[,s,c] <- sgi@mainSdTemplate[,s,c] / sqrt(tapply(rep(1,dim(D[selected,,,drop=FALSE])[1]),ft[selected],sum,na.rm=TRUE))
      sgi@mainSderrQuery[,s,c]    <- sgi@mainSdQuery[,s,c]    / sqrt(tapply(rep(1,dim(D)[1]),fq,sum,na.rm=TRUE))

      # main effect of negative controls
      I <- which(sgi@targets$group == "neg")
      IR <- which(sgi@reagents$TID %in% sgi@targets$TID[I])
      IT <- which(sgi@templateDesign$RID %in% sgi@reagents$RID[IR])
      IQ <- which(sgi@queryDesign$RID %in% sgi@reagents$RID[IR])
      if (length(IT) > 0) {
        sgi@mainNegTemplate[s,c] <- mean(sgi@mainTemplate[IT,s,c], na.rm=TRUE)
      } else {
        sgi@mainNegTemplate[s,c] <- median(sgi@mainTemplate[,s,c], na.rm=TRUE)
        warning("No negative controls in template design! Wild type effect can not be estimated! Median template main effect is used instead.")
      }
      if (length(IQ) > 0) {
        sgi@mainNegQuery[s,c] <- mean(sgi@mainQuery[IQ,s,c], na.rm=TRUE)
      } else {
        sgi@mainNegQuery[s,c] <- median(sgi@mainQuery[,s,c], na.rm=TRUE)
        warning("No negative controls in query design! Wild type effect can not be estimated! Median query main effect is used instead.")
      }
      sgi@mainNeg[s,c] <- sgi@mainNegTemplate[s,c] + sgi@mainNegQuery[s,c]
    }
  }

  a <- array(sgi@mainNegTemplate, dim=c(sgi@S,sgi@C,sgi@NT))
  sgi@mainTemplate[] <- sgi@mainTemplate - aperm(a,c(3,1,2))

  a <- array(sgi@mainNegQuery, dim=c(sgi@S,sgi@C,sgi@NQ))
  sgi@mainQuery[] <- sgi@mainQuery - aperm(a,c(3,1,2))

  return(sgi)
}

computePI <- function(sgi) {
  D <- getData(sgi, type="data", format="plain",drop=FALSE)
  I <- which((sgi@data2mainTemplate > 0) & (sgi@data2mainQuery > 0))
  d2t <- sgi@data2mainTemplate[I]
  d2q <- sgi@data2mainQuery[I]

  sgi@ni.model[] <- NA
  sgi@ni.model[I,,] <- sgi@mainTemplate[d2t,,,drop=FALSE] + sgi@mainQuery[d2q,,,drop=FALSE] + array(rep(sgi@mainNeg,each=length(d2q)),dim=c(length(d2q),sgi@S,sgi@C))

  sgi@pi[I,,] <- D - sgi@ni.model[I,,,drop=FALSE]

  return(sgi)
}

normalizeMainEffectQuery <- function(sgi, batch = NULL, time = NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (is.null(batch)) {
    batch = rep(1,sgi@NQ)
  }
  if (is.null(time)) {
    time <- 1:sgi@NQ
  }

  for (s in 1:sgi@S) {
    for (c in 1:sgi@C) {

      I = unique(batch)
      te <- rep(0,sgi@NQ)
      for (i in 1:length(I)) {
        J <- which(batch == I[i])
        x1 <- time[J]
        ##           x2 <- x1^2
        y <- sgi@mainQuery[J,s,c]
        ##           model <- rlm(y ~ x1+x2,k=0.3)
        ##           model <- rlm(y ~ x1,k=0.3)
        ##           te[J] = model$fitted.values
        x2 = x1
        y2 = y
        L = ceiling(length(x2)/2)
        for (l in 1:L) {
          model = lm(y2~x2)
          MaxRes = which.max(abs(model$res))
          x2=x2[-MaxRes]
          y2=y2[-MaxRes]
        }
        model = lm(y2~x2)
        te[J] <- predict(model,newdata=data.frame(x2=x1))
      }
      sgi@mainTimeEffect[,s,c] <- te
    }
  }

  return(sgi)
}

normalizeMainEffectTemplate <- function(sgi, screen = NULL, channel = NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

  use.screens = getScreenNames(sgi)
  if (!is.null(screen)) {
    use.screens = screen[screen %in% use.screens]
    if (sum(!(screen %in% use.screens)) > 0) {
      warning(sprintf("The following screen names are unknown: %s.", paste(screen[!(screen %in% use.screens)])))
    }
  }
  use.channels = getChannelNames(sgi)
  if (!is.null(channel)) {
    use.channels = channel[channel %in% use.channels]
    if (sum(!(channel %in% use.channels)) > 0) {
      warning(sprintf("The following channel names are unknown: %s.", paste(channel[!(channel %in% use.channels)],collapse=", ")))
    }
  }
  S <- length(use.screens)
  C <- length(use.channels)
  if (S < 1) { stop("At least one screen has to be selected.") }
  if (C < 1) { stop("At least one channel has to be selected.") }
  
  for (s in use.screens) {
    for (c in use.channels) {

      TP = unique(sgi@templateDesign$TemplatePlate)
      Mat = matrix(NA,nr=sgi@pdim[1],nc=sgi@pdim[2])
      row = match(substr(sgi@templateDesign$Well,1,1),LETTERS)
      col=as.integer(substr(sgi@templateDesign$Well,2,3))
      for (i in 1:length(TP)) {
        Mat[] = NA
        Ind = which(sgi@templateDesign$TemplatePlate == TP[i])
        for (j in Ind) {
          Mat[row[j],col[j]] = sgi@mainTemplate[j,s,c]
        }
        MP = medpolish(Mat,na.rm=TRUE)
        for (j in Ind) {
          sgi@mainSpatialEffect[j,s,c] <- MP$row[row[j]] + MP$col[col[j]]
        }
        sgi@mainSpatialEffectRow[i,,s,c] <- MP$row
        sgi@mainSpatialEffectCol[i,,s,c] <- MP$col
      }
    }
  }

  return(sgi)
}

