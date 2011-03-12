
summarizeScreens <- function(sgi, screens, newscreenname="mean") {
  stopifnot( is( sgi, "RNAinteract" ) )

  sgi2 <- sgi
  sgi2@data <- sgi2@data[,1,,drop=FALSE]
  sgi2@screenNames <- newscreenname
  sgi2@S <- 1L
  sgi2@data <- array(0, dim=c(sgi2@F,1,sgi2@C),dimnames = list(dimnames(sgi@data)[[1]],sgi2@screenNames,sgi2@channelNames))
  sgi2@mainTemplate      <- array(0, dim=c(sgi2@NT,1,sgi2@C),dimnames = list(1:sgi2@NT,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainQuery         <- array(0, dim=c(sgi2@NQ,1,sgi2@C),dimnames = list(1:sgi2@NQ,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainSderrTemplate <- array(0, dim=c(sgi2@NT,1,sgi2@C),dimnames = list(1:sgi2@NT,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainSderrQuery    <- array(0, dim=c(sgi2@NQ,1,sgi2@C),dimnames = list(1:sgi2@NQ,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainSdTemplate    <- array(0, dim=c(sgi2@NT,1,sgi2@C),dimnames = list(1:sgi2@NT,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainSdQuery       <- array(0, dim=c(sgi2@NQ,1,sgi2@C),dimnames = list(1:sgi2@NQ,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainTimeEffect    <- array(0, dim=c(sgi2@NQ,1,sgi2@C),dimnames = list(1:sgi2@NQ,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainSpatialEffect <- array(0, dim=c(sgi2@NT,sgi2@S,sgi2@C),dimnames = list(1:sgi2@NT,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainSpatialEffectRow <- array(0, dim=c(length(unique(sgi2@templateDesign$TemplatePlate)),sgi2@pdim[1],sgi2@S,sgi2@C),dimnames = list(NULL,NULL,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainSpatialEffectCol <- array(0, dim=c(length(unique(sgi2@templateDesign$TemplatePlate)),sgi2@pdim[2],sgi2@S,sgi2@C),dimnames = list(NULL,NULL,sgi2@screenNames,sgi2@channelNames))
  sgi2@mainNeg           <- array(0, dim=c(1,sgi2@C),dimnames = list(sgi2@screenNames,sgi2@channelNames))
  sgi2@mainNegTemplate   <- array(0, dim=c(1,sgi2@C),dimnames = list(sgi2@screenNames,sgi2@channelNames))
  sgi2@mainNegQuery      <- array(0, dim=c(1,sgi2@C),dimnames = list(sgi2@screenNames,sgi2@channelNames))
  sgi2@ni.model          <- array(NA, dim=c(sgi2@F,1,sgi2@C),dimnames = list(1:sgi2@F,sgi2@screenNames,sgi2@channelNames))
  sgi2@pi                <- array(NA, dim=c(sgi2@F,1,sgi2@C),dimnames = list(1:sgi2@F,sgi2@screenNames,sgi2@channelNames))
  sgi2@plateeffect       <- array(NA, dim=c(sgi2@F,1,sgi2@C),dimnames = list(1:sgi2@F,sgi2@screenNames,sgi2@channelNames))
  sgi2@p.value           <- sgi@p.value[,,1,,drop=FALSE]
  dimnames(sgi2@p.value)[[3]] <- newscreenname
  sgi2@p.value[] <- NA
  sgi2@q.value           <- sgi@q.value[,,1,,drop=FALSE]
  dimnames(sgi2@q.value)[[3]] <- newscreenname
  sgi2@q.value[] <- NA

  sgi2@data[,1,]              <- apply(sgi@data[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@mainTemplate[,1,]      <- apply(sgi@mainTemplate[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@mainQuery[,1,]         <- apply(sgi@mainQuery[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@mainSderrTemplate[,1,] <- apply(sgi@mainSderrTemplate[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@mainSderrQuery[,1,]    <- apply(sgi@mainSderrQuery[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@mainSdTemplate[,1,]    <- apply(sgi@mainSdTemplate[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@mainSdQuery[,1,]       <- apply(sgi@mainSdQuery[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@mainTimeEffect[,1,]    <- apply(sgi@mainTimeEffect[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@mainSpatialEffect[,1,] <- apply(sgi@mainSpatialEffect[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@mainSpatialEffectRow[,,1,] <- apply(sgi@mainSpatialEffectRow[,,screens,,drop=FALSE],c(1,2,4),mean,na.rm=TRUE)
  sgi2@mainSpatialEffectCol[,,1,] <- apply(sgi@mainSpatialEffectCol[,,screens,,drop=FALSE],c(1,2,4),mean,na.rm=TRUE)
  sgi2@mainNeg[1,]            <- apply(sgi@mainNeg[screens,],c(2),mean,na.rm=TRUE)
  sgi2@mainNegTemplate[1,]    <- apply(sgi@mainNegTemplate[screens,],c(2),mean,na.rm=TRUE)
  sgi2@mainNegQuery[1,]       <- apply(sgi@mainNegQuery[screens,],c(2),mean,na.rm=TRUE)
  sgi2@ni.model[,1,]          <- apply(sgi@ni.model[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@pi[,1,]                <- apply(sgi@pi[,screens,],c(1,3),mean,na.rm=TRUE)
  sgi2@plateeffect[,1,]       <- apply(sgi@plateeffect[,screens,],c(1,3),mean,na.rm=TRUE)
  if (any(!is.na(sgi@p.value[,,screens,]))) {
    sgi2@p.value[,,1,]          <- apply(sgi@p.value[,,screens,],c(1,2,4),min,na.rm=TRUE)
  }
  if (any(!is.na(sgi@q.value[,,screens,]))) {
    sgi2@q.value[,,1,]          <- apply(sgi@q.value[,,screens,],c(1,2,4),min,na.rm=TRUE)
  }

  return(sgi2)

}


bindscreens <- function(sgi1, sgi2) {
  stopifnot( is( sgi1, "RNAinteract" ) & is( sgi2, "RNAinteract" ) )

  sgi1@data                <- abind(sgi1@data,sgi2@data,along=2)
  sgi1@screenNames         <- c(sgi1@screenNames,sgi2@screenNames)
  sgi1@S                   <- sgi1@S + sgi2@S
  sgi1@mainTemplate        <- abind(sgi1@mainTemplate,sgi2@mainTemplate, along=2)
  sgi1@mainQuery           <- abind(sgi1@mainQuery,sgi2@mainQuery, along=2)
  sgi1@mainSderrTemplate   <- abind(sgi1@mainSderrTemplate,sgi2@mainSderrTemplate, along=2)
  sgi1@mainSderrQuery      <- abind(sgi1@mainSderrQuery,sgi2@mainSderrQuery, along=2)
  sgi1@mainSdTemplate      <- abind(sgi1@mainSdTemplate,sgi2@mainSdTemplate, along=2)
  sgi1@mainSdQuery         <- abind(sgi1@mainSdQuery,sgi2@mainSdQuery, along=2)
  sgi1@mainTimeEffect      <- abind(sgi1@mainTimeEffect,sgi2@mainTimeEffect, along=2)
  sgi1@mainSpatialEffect   <- abind(sgi1@mainSpatialEffect,sgi2@mainSpatialEffect, along=2)
  sgi1@mainSpatialEffectRow<- abind(sgi1@mainSpatialEffectRow,sgi2@mainSpatialEffectRow, along=3)
  sgi1@mainSpatialEffectCol<- abind(sgi1@mainSpatialEffectCol,sgi2@mainSpatialEffectCol, along=3)
  sgi1@mainNeg             <- abind(sgi1@mainNeg,sgi2@mainNeg, along=1)
  sgi1@mainNegTemplate     <- abind(sgi1@mainNegTemplate,sgi2@mainNegTemplate, along=1)
  sgi1@mainNegQuery        <- abind(sgi1@mainNegQuery,sgi2@mainNegQuery, along=1)
  sgi1@ni.model            <- abind(sgi1@ni.model,sgi2@ni.model, along=2)
  sgi1@pi                  <- abind(sgi1@pi,sgi2@pi, along=2)
  sgi1@plateeffect         <- abind(sgi1@plateeffect,sgi2@plateeffect, along=2)
  sgi1@p.value             <- abind(sgi1@p.value,sgi2@p.value, along=3)
  sgi1@q.value             <- abind(sgi1@q.value,sgi2@q.value, along=3)

  return(sgi1)
}

sgisubset <- function(sgi, screen=getScreenNames(sgi), channel=getChannelNames(sgi)) {
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

  sgi@data                <- sgi@data[,use.screens,use.channels,drop=FALSE]
  sgi@screenNames         <- use.screens
  sgi@S                   <- length(use.screens)
  sgi@channelNames        <- use.channels
  sgi@C                   <- length(use.channels)
##                     transformation    = "character",
  sgi@mainTemplate        <- sgi@mainTemplate[,use.screens,use.channels,drop=FALSE]
  sgi@mainQuery           <- sgi@mainQuery[,use.screens,use.channels,drop=FALSE]
  sgi@mainSderrTemplate   <- sgi@mainSderrTemplate[,use.screens,use.channels,drop=FALSE]
  sgi@mainSderrQuery      <- sgi@mainSderrQuery[,use.screens,use.channels,drop=FALSE]
  sgi@mainSdTemplate      <- sgi@mainSdTemplate[,use.screens,use.channels,drop=FALSE]
  sgi@mainSdQuery         <- sgi@mainSdQuery[,use.screens,use.channels,drop=FALSE]
  sgi@mainTimeEffect      <- sgi@mainTimeEffect[,use.screens,use.channels,drop=FALSE]
  sgi@mainSpatialEffect   <- sgi@mainSpatialEffect[,use.screens,use.channels,drop=FALSE]
  sgi@mainSpatialEffectRow<- sgi@mainSpatialEffectRow[,,use.screens,use.channels,drop=FALSE]
  sgi@mainSpatialEffectCol<- sgi@mainSpatialEffectCol[,,use.screens,use.channels,drop=FALSE]
  sgi@mainNeg             <- sgi@mainNeg[use.screens,use.channels,drop=FALSE]
  sgi@mainNegTemplate     <- sgi@mainNegTemplate[use.screens,use.channels,drop=FALSE]
  sgi@mainNegQuery        <- sgi@mainNegQuery[use.screens,use.channels,drop=FALSE]
  sgi@ni.model            <- sgi@ni.model[,use.screens,use.channels,drop=FALSE]
  sgi@pi                  <- sgi@pi[,use.screens,use.channels,drop=FALSE]
  sgi@plateeffect         <- sgi@plateeffect[,use.screens,use.channels,drop=FALSE]
  sgi@p.value             <- sgi@p.value[,,use.screens,use.channels,drop=FALSE]
  sgi@q.value             <- sgi@q.value[,,use.screens,use.channels,drop=FALSE]

  return(sgi)
}

sgisubsetQueryDesign <- function(sgi, query.targets = NULL, query.reagents = NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

# Does not reorder reagents and targets
  use.queryRID <- sgi@reagents$RID
  if (!is.null(query.reagents)) {
    if (!is.list(query.reagents)) {
      query.reagents = list(RID = query.reagents)
    }
  }
  if (!is.null(query.targets)) {
    use.TID = rep(FALSE,nrow(sgi@targets))
    if (!is.list(query.targets)) {
      query.targets = list(TID = query.targets)
    }
    if (is.list(query.targets)) {
      for (i in 1:length(query.targets)) {
        if ((i == 1) & (is.null(names(query.targets)[1]))) {
          names(query.targets)[1] = "TID"
        }
        if (names(query.targets)[i] %in% colnames(sgi@targets)) {
          useit = rep(FALSE,nrow(sgi@targets))
          if (names(query.targets)[i] == "TID") {
            if (is.logical(query.targets[[i]])) {
              useit[] = query.targets[[i]]
            } else if (is.numeric(query.targets[[i]])) {
              useit[query.targets[[i]]] = TRUE
            } else if (is.character(query.targets[[i]])) {
              useit[sgi@targets$TID %in% query.targets[[i]]] = TRUE
            } else {
              warning(sprintf("Index of TID is of mode %s. Cannot be used for subsetting.",mode(query.targets[[i]])))
            }
          } else {
            useit[] = sgi@targets[[names(query.targets)[i]]] %in% query.targets[[i]]
          }
          use.TID = use.TID | useit
        } else {
          warning(sprintf("column name %s not found in sgi@targets. Not used for subsetting.",names(query.targets)[i]))
        }
      }
    }
    if (is.null(query.reagents)) {
      query.reagents = list(TID=sgi@targets$TID[use.TID])
    } else {
      query.reagents = c(query.reagents,list(TID=sgi@targets$TID[use.TID]))
    }
  }
  if (!is.null(query.reagents)) {
    use.queryRID = rep(FALSE,nrow(sgi@reagents))
    if (is.list(query.reagents)) {
      for (i in 1:length(query.reagents)) {
        if ((i == 1) & (is.null(names(query.reagents)[1]))) {
          names(query.reagents)[1] = "RID"
        }
        if (names(query.reagents)[i] %in% colnames(sgi@reagents)) {
          useit = rep(FALSE,nrow(sgi@reagents))
          if (names(query.reagents)[i] == "RID") {
            if (is.logical(query.reagents[[i]])) {
              useit[] = query.reagents[[i]]
            } else if (is.numeric(query.reagents[[i]])) {
              useit[query.reagents[[i]]] = TRUE
            } else if (is.character(query.reagents[[i]])) {
              useit[sgi@reagents$RID %in% query.reagents[[i]]] = TRUE
            } else {
              warning(sprintf("Index of RID is of mode %s. Cannot be used for subsetting.",mode(query.reagents[[i]])))
            }
          } else {
            useit[] = sgi@reagents[[names(query.reagents)[i]]] %in% query.reagents[[i]]
          }
          use.queryRID = use.queryRID | useit
        } else {
          warning(sprintf("column name %s not found in sgi@reagents. Not used for subsetting.",names(query.reagents)[i]))
        }
      }
    }
    use.queryRID <- sgi@reagents$RID[use.queryRID]
    message(sprintf("%d query reagents covering %d target gene selected.", length(use.queryRID), length(unique(sgi@reagents[use.queryRID,"TID"]))))
  }

  use.data = which(sgi@reagents$RID[sgi@data2mainQuery] %in% use.queryRID)
  sgi@data                <- sgi@data[use.data,,,drop=FALSE]
  sgi@well                <- sgi@well[use.data]
  sgi@plate               <- sgi@plate[use.data]
  sgi@F                   <- length(use.data)
  sgi@mainQuery           <- sgi@mainQuery[sgi@queryDesign$RID %in% use.queryRID,,,drop=FALSE]
  sgi@mainSderrQuery      <- sgi@mainSderrQuery[sgi@queryDesign$RID %in% use.queryRID,,,drop=FALSE]
  sgi@mainSdQuery         <- sgi@mainSdQuery[sgi@queryDesign$RID %in% use.queryRID,,,drop=FALSE]
  sgi@mainTimeEffect      <- sgi@mainTimeEffect[sgi@queryDesign$RID %in% use.queryRID,,,drop=FALSE]
  sgi@data2mainTemplate   <- sgi@data2mainTemplate[use.data]
  sgi@data2mainQuery      <- match(sgi@data2mainQuery[use.data], which(sgi@queryDesign$RID %in% use.queryRID))
  sgi@queryDesign         <- sgi@queryDesign[sgi@queryDesign$RID %in% use.queryRID,]
  sgi@NQ                  <- nrow(sgi@queryDesign)
  sgi@ni.model            <- sgi@ni.model[use.data,,,drop=FALSE]
  sgi@pi                  <- sgi@pi[use.data,,,drop=FALSE]
  sgi@plateeffect         <- sgi@plateeffect[use.data,,,drop=FALSE]
  templateSymbol          <- sgi@targets$Symbol[sgi@targets$TID %in% sgi@reagents[sgi@templateDesign$RID,"TID"]]
  querySymbol             <- sgi@targets$Symbol[sgi@targets$TID %in% sgi@reagents[sgi@queryDesign$RID,"TID"]]
  sgi@p.value             <- array(NA, dim=c(length(templateSymbol),length(querySymbol),sgi@S,sgi@C),
                                   dimnames = list(template=templateSymbol,query=querySymbol,screen=sgi@screenNames,channel=sgi@channelNames))
  sgi@q.value             <- array(NA, dim=c(length(templateSymbol),length(querySymbol),sgi@S,sgi@C),
                                   dimnames = list(template=templateSymbol,query=querySymbol,screen=sgi@screenNames,channel=sgi@channelNames))

  sgi@mainTemplate[] = 0
  sgi@mainQuery[] = 0
  sgi@mainSderrTemplate[] = 0
  sgi@mainSderrQuery[] = 0
  sgi@mainSdTemplate[] = 0
  sgi@mainSdQuery[] = 0
  sgi@mainTimeEffect[] = 0
  sgi@mainSpatialEffect[] = 0
  sgi@mainSpatialEffectRow[] = 0
  sgi@mainSpatialEffectCol[] = 0
  sgi@mainNeg[] = 0
  sgi@mainNegTemplate[] = 0
  sgi@mainNegQuery[] = 0
  sgi@ni.model[] = NA
  sgi@pi[] = NA
  sgi@plateeffect[] = NA
  sgi@p.value[] = NA
  sgi@q.value[] = NA
  
  return(sgi)
}




