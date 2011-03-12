
createCellHTSFromFiles <- function(filePlatelist = "Platelist.txt",
                                   name="anonymous",
                                   path=".",pdim = NULL) {
  platelist <- read.table(sprintf("%s/%s",path,filePlatelist), sep="\t", header=TRUE,stringsAsFactors=FALSE)
  fn <- unique(platelist$Filename)

  data <- cbind(Filename=fn[1],read.table(sprintf("%s/%s",path,fn[1]), sep="\t", header=TRUE,stringsAsFactors=FALSE))
  if (length(fn) > 1) {
    for (i in 2:(length(fn))) {
      data <- rbind(data,cbind(Filename=fn[i],read.table(sprintf("%s/%s",path,fn[i]), sep="\t", header=TRUE,stringsAsFactors=FALSE)))
    }
  }
  if (any(duplicated(data[,c("Platebarcode","Well")]))) {
    stop("The columns Platebarcode and Well have to uniquely identify an entry in the input data files")
  }
  m <- match(data$Platebarcode,platelist$Platebarcode)
  if (any(is.na(m))) {
    I <- which(is.na(m))
    tmp <- data[I,c("Filename","Platebarcode")]
    tmp <- unique(tmp)
    message("The following bar codes are not found in platelist.\n")
    message("Filename\tPlatebarcode")
    message(paste(sprintf("%s\t%s",tmp$Filename,df$Platebarcode),collapse="\n"))
  }
  data <- cbind(Plate=platelist$Plate[m],Screen=platelist$Replicate[m], data,stringsAsFactors=FALSE)

  # Estimate plate dimension
  A <- regexpr("[0123456789]",data$Well)
  B <- nchar(data$Well)

  Row <- match(substr(data$Well,1,A-1),LETTERS)
  Col <- as.integer(substr(data$Well,A,B))
  if (is.null(pdim)) {
    mr <- max(Row)
    mc <- max(Col)

    if (!all(is.finite(c(mr,mc)))) {
      stop("The plate dimension can not be estimated. Please specify pdim when calling readMultiWellPlates")
    } else {
      pdim <- c(nrow=4L,ncol=6L)
      if ((mr > 4) | (mc > 6)) {
        pdim <- c(nrow=8L,ncol=12L)
      }
      if ((mr > 8) | (mc > 12)) {
        pdim <- c(nrow=16L,ncol=24L)
      }
    }
    message(sprintf("It looks like %d-well plates are used. nrow=%d ncol=%d",prod(pdim), pdim[1], pdim[2]))
    message("If this is not correct, specify the pdim argument and re-run.")
  } else {
    pdim = as.integer(pdim)
  }

  nscreen = max(platelist$Replicate)
  nwell <- prod(pdim)
  nplate <- max(data$Plate)
  id <- (data$Plate - 1) * nwell + (Row - 1) * pdim[2] + Col

  # create assayData
  I <- which(!(names(data) %in% c("Plate","Screen", "Filename", "Platebarcode", "Well")) & sapply(data,is.numeric))
  nchannel <- length(I)
  dat <- list()
  for (i in 1:nchannel) {
    dat[[names(data)[I[i]] ]] <- matrix(NA, nrow=nplate*nwell, ncol=nscreen)
    for (j in 1:nscreen) {
      J <- which(data$Screen == j)
      if (length(J) > 0) {
        dat[[names(data)[I[i]] ]][id[J],j] <- data[J,I[i]]
      } else {
        message(sprintf("No data found for channel %s and screen %d",names(data)[I[i]], j))
      }
    }
  }
  adata <- do.call(assayDataNew, c(storage.mode="lockedEnvironment", dat))

  # create phenoData
  pdata <- new( getClassDef( "AnnotatedDataFrame", package="Biobase"),
               data=data.frame(replicate=seq_len(nscreen),
                               assay=rep(name, nscreen),
                               stringsAsFactors=FALSE),
               varMetadata=data.frame(labelDescription=c("Replicate number",
                                        "Biological assay"),
               channel=factor(rep("_ALL_", 2L),
               levels=c(names(dat), "_ALL_")),
               row.names=c("replicate", "assay"),
               stringsAsFactors=FALSE))

  # create featureData
  fplate <- rep(1:nplate,each=nwell)
  fwell <- rep(sprintf("%s%0.2d",rep(LETTERS[1:pdim[1]],each=pdim[2]),1:pdim[2]),times = nplate)
  fdata <- new( getClassDef( "AnnotatedDataFrame", package="Biobase"), 
                 data <- data.frame(plate=fplate,
                                    well=fwell,
                                    controlStatus=factor(rep("unknown", nwell*nplate)),
                                    stringsAsFactors=FALSE),
                 varMetadata=data.frame(labelDescription=c("Plate number", "Well ID",
                                                           "Well annotation"),
                                        row.names=c("plate", "well", "controlStatus"),
                                        stringsAsFactors=FALSE))

  intensityFiles <- NULL
  batch <- NULL
  
  res <- new(  "cellHTS", 
             assayData=adata,
             phenoData=pdata,
             featureData=fdata,
             plateList=platelist)
}

