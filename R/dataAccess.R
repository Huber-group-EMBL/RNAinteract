
getScreenNames <- function(sgi) {
  stopifnot( is( sgi, "RNAinteract" ) )

  return(sgi@screenNames)
}

getChannelNames <- function(sgi) {
  stopifnot( is( sgi, "RNAinteract" ) )

  return(sgi@channelNames)
}

getData <- function(sgi,
                    type="data",
                    format = "plain",
                    design = "template",
                    mixTemplateQuery = TRUE,
                    screen=NULL, channel=NULL,
                    do.trafo=TRUE, do.inv.trafo = FALSE, normalized = FALSE,
                    withoutgroups=c(), drop=TRUE) {
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

  if ((S == 1) && (C == 1) && drop) {
    drop = TRUE
  } else {
    drop = FALSE
  }

  if (type %in% c("p.value","q.value")) {
##     if (!(format %in% c("targetMatrix","targetTQMatrix"))) {
##       stop("format has to be targetMatrix or targetTQMatrix for type p.value or q.value")
##     }
    res <- switch(type,
                  p.value = sgi@p.value[,,use.screens,use.channels,drop=FALSE],
                  q.value = sgi@q.value[,,use.screens,use.channels,drop=FALSE])
    lab <- switch(type,
                  p.value = rep("p-value",C),
                  q.value = rep("q-value",C))
##     if (format == "targetTQMatrix") {
##       namesr <- dimnames(sgi@p.value)[[1]]   # unique(sgi@targets[sgi@reagents[sgi@templateDesign$RID,"TID"],"Symbol"])
##       namesc <- dimnames(sgi@p.value)[[2]]   # unique(sgi@targets[sgi@reagents[sgi@queryDesign$RID,"TID"],"Symbol"])
##       R <- match(namesr, dimnames(res)[[1]])
##       C <- match(namesc, dimnames(res)[[2]])
##       res <- res[R,C,,,drop=FALSE]      
##     }
    if (length(withoutgroups) > 0) {
      gr = sgi@targets[dimnames(res)[[1]],"group"]
      gc = sgi@targets[dimnames(res)[[2]],"group"]
      res <- res[!(gr %in% withoutgroups),!(gc %in% withoutgroups),,,drop=FALSE]
    }
    if (drop) {
      res = res[,,1,1]
    }
    attr(res,"axislabel") <- lab
    names(attr(res,"axislabel")) <- use.channels
    return(res)
  }

  # Get a data array
  D <- switch(type,
              data      = sgi@data[,use.screens,use.channels,drop=FALSE],
              pi        = sgi@pi[,use.screens,use.channels,drop=FALSE],
              plateeffect = sgi@plateeffect[,use.screens,use.channels,drop=FALSE],
              ni.model  = sgi@ni.model[,use.screens,use.channels,drop=FALSE],
              main      = switch(design,
                template = sgi@mainTemplate[sgi@data2mainTemplate,use.screens,use.channels,drop=FALSE],
                query = sgi@mainQuery[sgi@data2mainQuery,use.screens,use.channels,drop=FALSE]),
              mainsderr = switch(design,
                template = sgi@mainSderrTemplate[sgi@data2mainTemplate,use.screens,use.channels,drop=FALSE],
                query = sgi@mainSderrQuery[sgi@data2mainQuery,use.screens,use.channels,drop=FALSE]),
              mainsd    = switch(design,
                template = sgi@mainSdTemplate[sgi@data2mainTemplate,use.screens,use.channels,drop=FALSE],
                query = sgi@mainSdQuery[sgi@data2mainQuery,use.screens,use.channels,drop=FALSE]),
              maintime  = switch(design,
                template = NULL,
                query = sgi@mainTimeEffect[sgi@data2mainQuery,use.screens,use.channels,drop=FALSE]),
              mainspatial  = switch(design,
                template = sgi@mainSpatialEffect[sgi@data2mainTemplate,use.screens,use.channels,drop=FALSE]),
                query = NULL,
              p.value   = NULL,  # Will be assigned in format "templateMatrix"
              q.value   = NULL,  # Will be assigned in format "templateMatrix"
              stop("getData: type unknown")
              )

  if ((type == "main") & (normalized)) {
    if (design == "template") {
      D = D - sgi@mainSpatialEffect[sgi@data2mainTemplate,use.screens,use.channels,drop=FALSE]
    }
    if (design == "query") {
      D = D - sgi@mainTimeEffect[sgi@data2mainQuery,use.screens,use.channels,drop=FALSE]
    }
  }
  
  # transform the data and memorize the scale of the data
  scale <- rep("linear scale", C)
  names(scale) <- use.channels
  if ((type == "data") & do.trafo) {
    for (c in use.channels) {
      fct <- getTransformation(sgi, channel = c)
      D[,,c] <- fct(D[,,c])
      scale[c] <- getScale(sgi, channel = c)
    }
  }

  # normalize the data (plate effects and query time effect)
  if ((type == "data") & normalized) {
    for (s in use.screens) {
      for (c in use.channels) {
        if (!do.trafo) {
          x <- transform(sgi, D[,s,c], channel = c)
          x <- x - sgi@mainTimeEffect[sgi@data2mainQuery,s,c] - sgi@mainSpatialEffect[sgi@data2mainTemplate,s,c]
          d[,s,c] <- invtransform(sgi, D[,s,c], channel = c)
        } else {
          D[,s,c] <- D[,s,c] - sgi@mainTimeEffect[sgi@data2mainQuery,s,c] - sgi@mainSpatialEffect[sgi@data2mainTemplate,s,c]
        }
      }
    }
  }

  # inverse transform the data and memorize the scale of the data
  if (!(type %in% c("data","p.value","q.value"))) {
    if (do.inv.trafo) {
      for (ch in use.channels) {
        fct = getInvTransformation(sgi, channel = ch)
        D[,,ch] = fct(D[,,ch])
      }
    } else {
      for (ch in use.channels) {
        scale[ch] <- getScale(sgi, channel = ch)
      }
    }
  }

  # reformating data
  if (format == "plain") {
    res <- D[,use.screens,use.channels,drop=drop]
  }
  if (format == "platelist") {
    P <- max(sgi@plate)
    platelist <- list()
    for (s in use.screens) {
      platelist[[s]] <- list()
      for (c in use.channels) {
        platelist[[s]][[c]] <- list()
        for (p in 1:P) {
##           platelist[[s]][[c]][[p]] <- matrix(d[sgi@plate == p,s,c],nr=sgi@pdim[1],nc=sgi@pdim[2])
          platelist[[s]][[c]][[p]] <- D[sgi@plate == p,s,c]
##           dim(platelist[[s]][[c]][[p]]) <- sgi@pdim
        }
      }
    }
    if (drop) {
      platelist <- platelist[[1]][[1]]
    }
    res <- platelist
  }
  if (format %in% c("targetMatrix")) {
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

    ## PI <- getData(sgi, type="pi", format = "plain",drop=FALSE)
    maxid = max(ID)
    maxn = max(table(id))

    res <- array(NA,dim=c(dim(ID),S,C),
                 dimnames = list(template=dimnames(ID)[[1]], query=dimnames(ID)[[2]],
                   screen = use.screens, channel = use.channels))
    for (s in 1:S) {
      for (c in 1:C) {
        X = matrix(NA,nr=maxid,nc=maxn)
        ## pi = PI[,s,c]
        id2 = id
        for (i in 1:maxn) {
          m = match(1:maxid,id2)
          X[,i] = D[m,s,c]
          id2[m] = NA
        }
        d <- apply(X,1,mean,na.rm=TRUE)
        res[,,s,c] <- d[ID]
      }
    }
    if (length(withoutgroups) > 0) {
      gr = sgi@targets$group[match(dimnames(ID)[[1]], sgi@targets$Symbol)]
      gc = sgi@targets$group[match(dimnames(ID)[[2]], sgi@targets$Symbol)]
      res <- res[!(gr %in% withoutgroups),!(gc %in% withoutgroups),,,drop=FALSE]
    }
    if (drop) {
      res <- res[,,,]
    }
  }
##   if (format %in% c("reagentTQMatrix","reagentMatrix","targetTQMatrix","targetMatrix")) {
##     if (format == "reagentTQMatrix") {
##       namesr <- unique(sgi@templateDesign$RID[sgi@data2mainTemplate])
##       namesc <- unique(sgi@queryDesign$RID[sgi@data2mainQuery])
##       dr <- sgi@templateDesign$RID[sgi@data2mainTemplate]
##       dc <- sgi@queryDesign$RID[sgi@data2mainQuery]
##       mr <- match(dr, namesr)
##       mc <- match(dc, namesc)
##     }
##     if (format == "targetTQMatrix") {
##       namesr <- dimnames(sgi@p.value)[[1]]   
##       namesc <- dimnames(sgi@p.value)[[2]]   
## ##       namesr <- unique(sgi@targets[sgi@reagents[sgi@templateDesign$RID,"TID"],"Symbol"])
## ##       namesc <- unique(sgi@targets[sgi@reagents[sgi@queryDesign$RID,"TID"],"Symbol"])
##       dr <- sgi@targets[sgi@reagents[sgi@templateDesign$RID[sgi@data2mainTemplate],"TID"],"Symbol"]
##       dc <- sgi@targets[sgi@reagents[sgi@queryDesign$RID[sgi@data2mainQuery],"TID"],"Symbol"]
##       mr <- match(dr, namesr)
##       mc <- match(dc, namesc)
##     }
##     if (format == "reagentMatrix") {
##       namesr <- sgi@reagents$RID
##       namesc <- sgi@reagents$RID
##       dr <- sgi@templateDesign$RID[sgi@data2mainTemplate]
##       dc <- sgi@queryDesign$RID[sgi@data2mainQuery]
##       mr <- match(dr, namesr)
##       mc <- match(dc, namesc)
##       I = which(mr < mc)
##       if (length(I) > 0) {
##         tmp <- mr[I]
##         mr[I] <- mc[I]
##         mc[I] <- tmp
##       }
##     }
##     if (format == "targetMatrix") {
## ##       namesr <- sgi@targets$Symbol
## ##       namesc <- sgi@targets$Symbol
##       namesr <- dimnames(sgi@p.value)[[1]]   
##       namesc <- dimnames(sgi@p.value)[[2]]   
##       dr <- sgi@targets[sgi@reagents[sgi@templateDesign$RID[sgi@data2mainTemplate],"TID"],"Symbol"]
##       dc <- sgi@targets[sgi@reagents[sgi@queryDesign$RID[sgi@data2mainQuery],"TID"],"Symbol"]
##       mr <- match(dr, namesr)
##       mc <- match(dc, namesc)
##       I = which(mr < mc)
##       if (length(I) > 0) {
##         tmp <- mr[I]
##         mr[I] <- mc[I]
##         mc[I] <- tmp
##       }
##     }
##     nr <- length(namesr)
##     nc <- length(namesc)

##     d <- dim(D)
##     IDXT <- array(mr,dim=d)
##     IDXQ <- array(mc-1,dim=d)
##     IDXS <- aperm(array((1:d[2])-1,dim=d[c(2,1,3)]), c(2,1,3))
##     IDXC <- aperm(array((1:d[3])-1,dim=d[c(3,1,2)]), c(2,3,1))
##     IDX <- IDXT + nr * IDXQ + nr * nc * IDXS + nr * nc * d[2] * IDXC

##     A <- array(NA, dim=c(nr, nc, S, C))
##     A[] <- tapply(D,factor(IDX, levels=1:prod(c(nr,nc,S,C))), mean, na.rm=TRUE)
##     if (format %in% c("reagentTQMatrix","targetTQMatrix")) {
##       dimnames(A) <- list(template = namesr, query = namesc, screens = use.screens, channels = use.channels)
##     } else {
##       dimnames(A) <- list(targets = namesr, targets = namesc, screens = use.screens, channels = use.channels)
##     }
##     if (format %in% c("reagentMatrix","targetMatrix")) {
##       for (s in use.screens) {
##         for (ch in use.channels) {
##           B <- A[,,s,ch]
##           B[upper.tri(B)] <- t(B)[upper.tri(B)]
##           A[,,s,ch] <- B
##         }
##       }
##     }
##     if (length(withoutgroups) > 0) {
##       if (format %in% c("reagentMatrix","reagentTQMatrix")) {
##         gr = sgi@targets$group[match(sgi@reagents[namesr,"TID"], sgi@targets$Symbol)]
##         gc = sgi@targets$group[match(sgi@reagents[namesc,"TID"], sgi@targets$Symbol)]
##         A <- A[!(gr %in% withoutgroups),!(gc %in% withoutgroups),,,drop=FALSE]
##       } else {
##         gr = sgi@targets$group[match(namesr, sgi@targets$Symbol)]
##         gc = sgi@targets$group[match(namesc, sgi@targets$Symbol)]
##         A <- A[!(gr %in% withoutgroups),!(gc %in% withoutgroups),,,drop=FALSE]
##       }
##     }
##     if (drop) {
##       A = A[,,1,1]
##     }
##     res <- A
##   }
  L <- switch(type,
              data      = "input data",
              pi        = "pairwise interaction",
              plateeffect = "plateeffect",
              ni.model  = "non-interacting model value",
              mainraw   = switch(design,
                template = "raw main effect template reagent",
                query = "raw main effect query reagent"),
              main      = switch(design,
                template = "main effect (ratio), template reagent",
                query = "main effect (ratio), query reagent"),
              mainabs   = switch(design,
                template = "main effect (absolute), template reagent",
                query = "main effect (absolute), query reagent"),
              mainsderr = switch(design,
                template = "main effect stderr, template reagent",
                query = "main effect stderr, query reagent"),
              mainsd    = switch(design,
                template = "main effect stddev, query reagent",
                query = "main effect stddev, query reagent"),
              maintime  = switch(design,
                template = NULL,
                query = "main effect stderr, query reagent")
              )
  attr(res,"axislabel") <- sprintf("%s (%s)", L, scale)
  names(attr(res,"axislabel")) <- use.channels
  return(res)
}

## getCellHTSObject <- function(sgi, type="data") {
##   stopifnot( is( sgi, "RNAinteract" ) )

##   use.screens <- getScreenNames(sgi)
##   use.channels <- getChannelNames(sgi)
  
##   D <- switch(type,
##               data      = sgi@data[,use.screens,use.channels,drop=FALSE],
##               pi        = sgi@pi[,use.screens,use.channels,drop=FALSE],
##               plateeffect = sgi@plateeffect[,use.screens,use.channels,drop=FALSE],
##               ni.model  = sgi@ni.model[,use.screens,use.channels,drop=FALSE],
##               main      = switch(design,
##                 template = sgi@mainTemplate[sgi@data2mainTemplate,use.screens,use.channels,drop=FALSE],
##                 query = sgi@mainQuery[sgi@data2mainQuery,use.screens,use.channels,drop=FALSE]),
##               mainsderr = switch(design,
##                 template = sgi@mainSderrTemplate[sgi@data2mainTemplate,use.screens,use.channels,drop=FALSE],
##                 query = sgi@mainSderrQuery[sgi@data2mainQuery,use.screens,use.channels,drop=FALSE]),
##               mainsd    = switch(design,
##                 template = sgi@mainSdTemplate[sgi@data2mainTemplate,use.screens,use.channels,drop=FALSE],
##                 query = sgi@mainSdQuery[sgi@data2mainQuery,use.screens,use.channels,drop=FALSE]),
##               maintime  = switch(design,
##                 template = NULL,
##                 query = sgi@mainTimeEffect[sgi@data2mainQuery,use.screens,use.channels,drop=FALSE]),
##               mainspatial  = switch(design,
##                 template = sgi@mainSpatialEffect[sgi@data2mainTemplate,use.screens,use.channels,drop=FALSE]),
##                 query = NULL,
##               stop("getData: type unknown")
##               )


##   # create assayData
##   dat <- list()
##   for (ch in 1:sgi@C) {
##     dat[[dimnames(D)[[3]][ch] ]] <- D[,,ch]
##   }
##   adata <- do.call(assayDataNew, c(storage.mode="lockedEnvironment", dat))

##   # create phenoData
##   pdata <- new("AnnotatedDataFrame",
##                data=data.frame(replicate=seq_len(sgi@S),
##                                assay=rep("assay name", sgi@S),
##                                stringsAsFactors=FALSE),
##                varMetadata=data.frame(labelDescription=c("Replicate number",
##                                         "Biological assay"),
##                channel=factor(rep("_ALL_", 2L),levels=c(getChannelNames(sgi), "_ALL_")),
##                row.names=c("replicate", "assay"),
##                stringsAsFactors=FALSE))

##   # create featureData
##   fdata <- new("AnnotatedDataFrame", 
##                  data <- data.frame(plate=sgi@plate,
##                                     well=sgi@well,
##                                     controlStatus=factor(rep("unknown", sgi@F)),
##                                     stringsAsFactors=FALSE),
##                  varMetadata=data.frame(labelDescription=c("Plate number", "Well ID",
##                                                            "Well annotation"),
##                                         row.names=c("plate", "well", "controlStatus"),
##                                         stringsAsFactors=FALSE))

##   batch <- NULL
##   res <- new("cellHTS", 
##              assayData=adata,
##              phenoData=pdata,
##              featureData=fdata)

##   return(res)
## }

getMain <- function(sgi, type = "main", design = "template", summary = "none", QueryNr = NULL, TemplatePlate = NULL, do.inv.trafo = FALSE, format = "plain", withoutgroups=c(), screen=NULL, channel=NULL, normalized = TRUE, drop=TRUE) {
  stopifnot( is( sgi, "RNAinteract" ) )

  use.screens = getScreenNames(sgi)
  if (!is.null(screen)) {
    use.screens = screen[screen %in% use.screens]
    if (sum(!(screen %in% use.screens)) > 0) {
      message(sprintf("The following screen names are unknown: %s", screen[!(screen %in% use.screens)]))
    }
  }
  use.channels = getChannelNames(sgi)
  if (!is.null(channel)) {
    use.channels = channel[channel %in% use.channels]
    if (sum(!(channel %in% use.channels)) > 0) {
      message(sprintf("The following channel names are unknown: %s", channel[!(channel %in% use.channels)]))
    }
  }
  S <- length(use.screens)
  C <- length(use.channels)
  if ((S == 1) && (C == 1) && drop) {
    drop = TRUE
  } else {
    drop = FALSE
  }

  M <- switch(type,
              main        = switch(design,
                                   template = sgi@mainTemplate[,use.screens,use.channels,drop=FALSE],
                                   query = sgi@mainQuery[,use.screens,use.channels,drop=FALSE]),
              mainsderr   = switch(design,
                                   template = sgi@mainSderrTemplate[,use.screens,use.channels,drop=FALSE],
                                   query = sgi@mainSderrQuery[,use.screens,use.channels,drop=FALSE]),
              mainsd      = switch(design,
                                   template = sgi@mainSdTemplate[,use.screens,use.channels,drop=FALSE],
                                   query = sgi@mainSdQuery[,use.screens,use.channels,drop=FALSE]),
              maintime    = switch(design,
                                   template = NULL,
                                   query = sgi@mainTimeEffect[,use.screens,use.channels,drop=FALSE]),
              mainspatial = switch(design,
                                   template = sgi@mainSpatialEffect[,use.screens,use.channels,drop=FALSE],
                                   query = NULL)
              )

  if ((type == "main") & normalized) {
    M <- M - switch(design,
                    template = sgi@mainSpatialEffect[,use.screens,use.channels,drop=FALSE],
                    query = sgi@mainTimeEffect[,use.screens,use.channels,drop=FALSE])
  }

  dimnames(M) <- switch(design,
                        template = list(sgi@templateDesign$RID,use.screens, use.channels),
                        query = list(sgi@queryDesign$RID,use.screens,use.channels) )

  if (summary == "reagent") {
    if (design == "template") {
      m <- match(sgi@templateDesign$RID, sgi@reagents$RID)
    } else {
      m <- match(sgi@queryDesign$RID, sgi@reagents$RID)
    }
    Mtmp <- array(NA, dim=c(nrow(sgi@reagents),S,C), dimnames=list(sgi@reagents$RID,use.screens,use.channels))
    f <- factor(m,levels = 1:nrow(sgi@reagents))
    for (s in use.screens) {
      for (c in use.channels) {
        Mtmp[,s,c] <- tapply(M[,s,c],f,mean)
      }
    }
    M <- Mtmp
  }
  
  if (summary == "target") {
    if (design == "template") {
      m1 <- match(sgi@templateDesign$RID, sgi@reagents$RID)
      m <- match(sgi@reagents$TID[m1],sgi@targets$TID)
    } else {
      m1 <- match(sgi@queryDesign$RID, sgi@reagents$RID)
      m <- match(sgi@reagents$TID[m1],sgi@targets$TID)
    }
    Mtmp <- array(NA, dim=c(nrow(sgi@targets),S,C), dimnames=list(sgi@targets$Symbol,use.screens,use.channels))
    f <- factor(m,levels = 1:nrow(sgi@targets))
    for (s in use.screens) {
      for (c in use.channels) {
        Mtmp[,s,c] <- tapply(M[,s,c],f,mean)
      }
    }
    M <- Mtmp
  }
  
  scale <- rep("linear scale", C)
  names(scale) <- use.channels
  if (do.inv.trafo) {
    for (c in use.channels) {
      fct = getInvTransformation(sgi, channel = c)
      M[,,c] = fct(M[,,c])
    }
  } else {
    for (ch in use.channels) {
      scale[ch] <- getScale(sgi, channel = ch)
    }
  }


  if (!is.null(QueryNr) | !is.null(TemplatePlate)) {
    IQ <- rep(TRUE, dim(M)[1])
    IT <- rep(TRUE, dim(M)[1])
    if (!is.null(QueryNr)) {
      IQ = switch(design,
        template = (sgi@templateDesign$QueryNr %in% QueryNr),
        query = (sgi@queryDesign$QueryNr %in% QueryNr) )
    }
    if (!is.null(TemplatePlate)) {
      IT <- switch(design,
                   template = (sgi@templateDesign$TemplatePlate %in% TemplatePlate),
                   query = (sgi@queryDesign$TemplatePlate %in% TemplatePlate) )
    }
    M <- M[IT & IQ,,,drop=FALSE]
  }

  if (length(withoutgroups) > 0) {
    if (format == "target") {
      g = sgi@targets$group[match(dimnames(M)[[1]], sgi@targets$Symbol)]
      M <- M[!(g %in% withoutgroups),,,drop=FALSE]
    } else {
      g = sgi@targets$group[match(sgi@reagents[dimnames(M)[[1]],"TID"], sgi@targets$Symbol)]
      M <- M[!(g %in% withoutgroups),,,drop=FALSE]
    }
  }

  if (drop) {
    M <- M[,1,1]
  }

  L <- switch(type,
              mainraw   = switch(design, template = "raw main effect, template",        query = "raw main effect query"),
              main      = switch(design, template = "main effect (ratio), template",    query = "main effect (ratio), query"),
              mainabs   = switch(design, template = "main effect (absolute), template", query = "main effect (absolute), query"),
              mainsderr = switch(design, template = "main effect stderr, template",     query = "main effect stderr, query"),
              mainsd    = switch(design, template = "main effect stddev, template",     query = "main effect stddev, query"),
              maintime  = switch(design, template = NULL,                               query = "main effect stderr, query")
              )

  sel <- ""
  if (!is.null(QueryNr)) { sel <- sprintf("%s, of query-nr %s", sel, paste(QueryNr,collapse=",")) }
  if (!is.null(TemplatePlate)) { sel <- sprintf("%s, on temp.-plate %s", sel, paste(TemplatePlate,collapse=",")) }
  attr(M,"axislabel") <- sprintf("%s%s (%s)", L, sel, scale)
  names(attr(M,"axislabel")) <- use.channels

  return(M)

}


getMainNeg <- function(sgi, type = "all", do.inv.trafo = FALSE, screen = NULL, channel = NULL, drop=TRUE) {

  use.screens = getScreenNames(sgi)
  if (!is.null(screen)) {
    use.screens = screen[screen %in% use.screens]
    if (sum(!(screen %in% use.screens)) > 0) {
      message(sprintf("The following screen names are unknown: %s", screen[!(screen %in% use.screens)]))
    }
  }
  use.channels = getChannelNames(sgi)
  if (!is.null(channel)) {
    use.channels = channel[channel %in% use.channels]
    if (sum(!(channel %in% use.channels)) > 0) {
      message(sprintf("The following channel names are unknown: %s", channel[!(channel %in% use.channels)]))
    }
  }
  S <- length(use.screens)
  C <- length(use.channels)
  if ((S == 1) && (C == 1) && drop) {
    drop = TRUE
  } else {
    drop = FALSE
  }

  nc <- switch(type,
               template = sgi@mainNegTemplate[use.screens,use.channels,drop=FALSE],
               query = sgi@mainNegQuery[use.screens,use.channels,drop=FALSE],
               sgi@mainNeg[use.screens,use.channels,drop=FALSE] )

  if (do.inv.trafo) {
    for (s in use.screens) {
      for (c in use.channels) {
        nc[s,c] <- invtransform(sgi,nc[s,c],channel=c)
      }
    }
  }

  if (drop) {
    nc <- nc[1,1]
  }
  
  return(nc)
}

getReplicateData <- function(sgi, screen, channel, type="data", design = "template", do.trafo = TRUE, do.inv.trafo = FALSE, normalized = FALSE) {
  stopifnot( is( sgi, "RNAinteract" ) )

  RID1 <- sgi@templateDesign$RID[sgi@data2mainTemplate]
  RID2 <- sgi@queryDesign$RID[sgi@data2mainQuery]
  I <- RID1 > RID2
  Rtmp <- RID2[I]
  RID2[I] <- RID1[I]
  RID1[I] <- Rtmp

  IDX <- sprintf("%s____%s",RID1,RID2)

  I <- order(IDX)
  IDX <- IDX[I]

  n <- length(IDX)
  J <- which(!is.na(IDX[1:(n-1)]) & !is.na(IDX[2:n]) & (IDX[1:(n-1)] == IDX[2:n]))
  A <- (I[1:n])[J]
  B <- (I[2:(n-1)])[J]

  D <- getData(sgi, type = type, design = design, do.trafo = do.trafo, do.inv.trafo = do.inv.trafo, format = "plain", screen = screen, channel = channel, normalized = normalized)
  lab <- attr(D, "axislabel")[channel]

  res = list(x = D[A], y = D[B], lab = lab)
  return(res)
}


getIndDesignData <- function(sgi, screen, channel, type="data", design = "template", do.trafo = TRUE, do.inv.trafo = FALSE, normalized = FALSE) {
  stopifnot( is( sgi, "RNAinteract" ) )

  RID1 <- sgi@templateDesign$RID[sgi@data2mainTemplate]
  RID2 <- sgi@queryDesign$RID[sgi@data2mainQuery]

  design <- rep(1,nrow(sgi@reagents))
  design[duplicated(sgi@reagents$TID)] <- 2

  TID1 <- sgi@reagents$TID[match(RID1,sgi@reagents$RID)]
  TID2 <- sgi@reagents$TID[match(RID2,sgi@reagents$RID)]
  D1 <- design[match(RID1,sgi@reagents$RID)]
  D2 <- design[match(RID2,sgi@reagents$RID)]

  I <- TID1 > TID2
  Ttmp <- TID2[I]
  TID2[I] <- TID1[I]
  TID1[I] <- Ttmp
  dtmp <- D2[I]
  D2[I] <- D1[I]
  D1[I] <- dtmp

  V1 <- which((D1 == 1) & (D2 == 1))
  V2 <- which((D1 == 2) & (D2 == 2))
  IDX <- sprintf("%s____%s",TID1,TID2)
  I1 <- V1[which(!duplicated(IDX[V1]))]
  I2 <- V2[which(!duplicated(IDX[V2]))]
  m <- match(IDX[I1],IDX[I2])
  A <- I1[which(is.finite(m))]
  B <- I2[m[which(is.finite(m))]]

  D <- getData(sgi, type = type, design = design, do.trafo = do.trafo, do.inv.trafo = do.inv.trafo, format = "plain", screen = screen, channel = channel, normalized = normalized)
  lab <- attr(D, "axislabel")[channel]

  res <- list(x = D[A], y = D[B], lab = lab)
  return(res)
}

