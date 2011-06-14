
createRNAinteract <- function(data,
                      well,
                      plate,
                      pdim,
                      Reagents,
                      Targets,
                      TemplateDesign,
                      QueryDesign,
                      Transformation = NULL) {
  sgi <- new("RNAinteract")

  sgi@data <- data
  sgi@screenNames <- dimnames(data)[[2]]
  sgi@channelNames <- dimnames(data)[[3]]
  sgi@well <- well
  sgi@plate <- plate
  sgi@pdim <- pdim
  sgi@NT <- nrow(TemplateDesign)
  sgi@NQ <- nrow(QueryDesign)
  sgi@C  <- length(dimnames(data)[[3]])
  sgi@S  <- length(dimnames(data)[[2]])
  sgi@F  <- nrow(data)

  A <- regexpr("[0123456789]",TemplateDesign$Well)
  B <- nchar(TemplateDesign$Well)
  Row <- substr(TemplateDesign$Well,1,A-1)
  Col <- as.integer(substr(TemplateDesign$Well,A,B))
  TemplateDesign$Well = sprintf("%s%0.2d",Row,Col)

  sgi@reagents          <- Reagents
  row.names(sgi@reagents) = sgi@reagents$RID
  sgi@targets           <- Targets
  row.names(sgi@targets) = sgi@targets$TID
  sgi@templateDesign    <- TemplateDesign
  sgi@queryDesign       <- QueryDesign
  sgi@transformation    <- rep("log2",sgi@C)
  names(sgi@transformation) <- sgi@channelNames
  if (!is.null(Transformation)) {
    sgi@transformation[names(Transformation)] <- Transformation
  }

  sgi@mainTemplate      <- array(0, dim=c(sgi@NT,sgi@S,sgi@C),dimnames = list(1:sgi@NT,sgi@screenNames,sgi@channelNames))
  sgi@mainQuery         <- array(0, dim=c(sgi@NQ,sgi@S,sgi@C),dimnames = list(1:sgi@NQ,sgi@screenNames,sgi@channelNames))
  sgi@mainSderrTemplate <- array(0, dim=c(sgi@NT,sgi@S,sgi@C),dimnames = list(1:sgi@NT,sgi@screenNames,sgi@channelNames))
  sgi@mainSderrQuery    <- array(0, dim=c(sgi@NQ,sgi@S,sgi@C),dimnames = list(1:sgi@NQ,sgi@screenNames,sgi@channelNames))
  sgi@mainSdTemplate    <- array(0, dim=c(sgi@NT,sgi@S,sgi@C),dimnames = list(1:sgi@NT,sgi@screenNames,sgi@channelNames))
  sgi@mainSdQuery       <- array(0, dim=c(sgi@NQ,sgi@S,sgi@C),dimnames = list(1:sgi@NQ,sgi@screenNames,sgi@channelNames))
  sgi@mainTimeEffect    <- array(0, dim=c(sgi@NQ,sgi@S,sgi@C),dimnames = list(1:sgi@NQ,sgi@screenNames,sgi@channelNames))
  sgi@mainSpatialEffect <- array(0, dim=c(sgi@NT,sgi@S,sgi@C),dimnames = list(1:sgi@NT,sgi@screenNames,sgi@channelNames))
  sgi@mainSpatialEffectRow <- array(0, dim=c(length(unique(TemplateDesign$TemplatePlate)),sgi@pdim[1],sgi@S,sgi@C),dimnames = list(NULL,NULL,sgi@screenNames,sgi@channelNames))
  sgi@mainSpatialEffectCol <- array(0, dim=c(length(unique(TemplateDesign$TemplatePlate)),sgi@pdim[2],sgi@S,sgi@C),dimnames = list(NULL,NULL,sgi@screenNames,sgi@channelNames))
  sgi@mainNeg           <- array(0, dim=c(sgi@S,sgi@C),dimnames = list(sgi@screenNames,sgi@channelNames))
  sgi@mainNegTemplate   <- array(0, dim=c(sgi@S,sgi@C),dimnames = list(sgi@screenNames,sgi@channelNames))
  sgi@mainNegQuery      <- array(0, dim=c(sgi@S,sgi@C),dimnames = list(sgi@screenNames,sgi@channelNames))
  sgi@data2mainTemplate <- rep(NA_integer_,sgi@F)
  sgi@data2mainQuery    <- rep(NA_integer_,sgi@F)
  sgi@ni.model          <- array(NA, dim=c(sgi@F,sgi@S,sgi@C),dimnames = list(1:sgi@F,sgi@screenNames,sgi@channelNames))
  sgi@pi                <- array(NA, dim=c(sgi@F,sgi@S,sgi@C),dimnames = list(1:sgi@F,sgi@screenNames,sgi@channelNames))
  sgi@plateeffect       <- array(NA, dim=c(sgi@F,sgi@S,sgi@C),dimnames = list(1:sgi@F,sgi@screenNames,sgi@channelNames))
  templateSymbol        <- sgi@targets$Symbol[sgi@targets$TID %in% sgi@reagents[sgi@templateDesign$RID,"TID"]]
  querySymbol           <- sgi@targets$Symbol[sgi@targets$TID %in% sgi@reagents[sgi@queryDesign$RID,"TID"]]
  sgi@p.value           <- array(NA, dim=c(length(templateSymbol),length(querySymbol),sgi@S,sgi@C),
                                 dimnames = list(template=templateSymbol,query=querySymbol,screen=sgi@screenNames,channel=sgi@channelNames))
  sgi@q.value           <- array(NA, dim=c(length(templateSymbol),length(querySymbol),sgi@S,sgi@C),
                                 dimnames = list(template=templateSymbol,query=querySymbol,screen=sgi@screenNames,channel=sgi@channelNames))
##   sgi@p.value           <- array(NA, dim=c(nrow(sgi@targets),nrow(sgi@targets),sgi@S,sgi@C),dimnames = list(sgi@targets$Symbol,sgi@targets$Symbol,sgi@screenNames,sgi@channelNames))
##   sgi@q.value           <- array(NA, dim=c(nrow(sgi@targets),nrow(sgi@targets),sgi@S,sgi@C),dimnames = list(sgi@targets$Symbol,sgi@targets$Symbol,sgi@screenNames,sgi@channelNames))

  for (i in 1:sgi@NQ) {
    IT <- which((TemplateDesign$TemplatePlate == QueryDesign$TemplatePlate[i]) & (TemplateDesign$QueryNr == QueryDesign$QueryNr[i]))
    IS <- which((plate == QueryDesign$Plate[i]) & (well %in% TemplateDesign$Well[IT]))
    if (length(IS) > 0) {
      sgi@data2mainQuery[IS] <- as.integer(i)
    }

    m <- match(well[IS], TemplateDesign$Well[IT])
    IS <- IS[!is.na(m)]
    m <- m[!is.na(m)]
    sgi@data2mainTemplate[IS] <- IT[m]

  }

  return(sgi)
}

createRNAinteractFromFiles <- function(name="anonymous",
                               filePlatelist = "Platelist.txt",
                               fileReagents = "Reagents.txt",
                               fileTargets = "Targets.txt",
                               fileTemplateDesign = "TemplateDesign.txt",
                               fileQueryDesign = "QueryDesign.txt",
                               path = ".",
                               pdim = NULL,
                               Transformation = "log2") {
  chts <- createCellHTSFromFiles(name=name, filePlatelist = filePlatelist, path = path, pdim = pdim)

  Reagents       <- read.table(sprintf("%s/%s",path,fileReagents),
                               comment.char="", quote="",
                               header=TRUE, sep="\t",stringsAsFactors=FALSE)
  Targets        <- read.table(sprintf("%s/%s",path,fileTargets),
                               comment.char="", quote="",
                               header=TRUE, sep="\t",stringsAsFactors=FALSE)
  TemplateDesign <- read.table(sprintf("%s/%s",path,fileTemplateDesign),
                               comment.char="", quote="",
                               header=TRUE, sep="\t",stringsAsFactors=FALSE)
  QueryDesign    <- read.table(sprintf("%s/%s",path,fileQueryDesign),
                               comment.char="", quote="",
                               header=TRUE, sep="\t",stringsAsFactors=FALSE)

  return(createRNAinteract(Data(chts), well(chts), plate(chts), pdim(chts), Reagents, Targets, TemplateDesign, QueryDesign, Transformation = Transformation))
}




