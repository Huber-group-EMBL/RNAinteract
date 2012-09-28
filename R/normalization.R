
normalizePlateEffect <- function(sgi, type="Bscore", maxit=20,verbose=0) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (verbose>0) {
    time1 <- Sys.time()
    message("Normalize plate effects. \r", appendLF = FALSE)
  }

  Plate <- sgi@plate + 96*(sgi@queryDesign$QueryNr[sgi@data2mainQuery]-1)
  Well <- sgi@well
  wellnames <- sprintf("%s%0.2d",rep(LETTERS[1:sgi@pdim[1]],times=sgi@pdim[2]),rep(1:sgi@pdim[2],each=sgi@pdim[1]))
  names(wellnames) <- wellnames
  wellnames[1:192+192] <- wellnames[1:192]
  Well <- wellnames[Well]
  pdim = c(16,12)
  
  D <- getData(sgi,type="pi",format="plain")
  plate <- matrix(NA,nrow=pdim[1],ncol=pdim[2])
  tmp <- rep(0,length(plate))
  names(tmp) <- sprintf("%s%0.2d",rep(LETTERS[1:pdim[1]],times=pdim[2]),rep(1:pdim[2],each=pdim[1]))
  row.names(plate) <- LETTERS[1:pdim[1]]
  colnames(plate) <- 1:pdim[2]

  X <- as.matrix(expand.grid(1:pdim[1],1:pdim[2]))
  X[,1] <- X[,1] - mean(X[,1])
  X[,2] <- X[,2] - mean(X[,2])
  X[,1] <- X[,1] / sd(X[,1])
  X[,2] <- X[,2] / sd(X[,2])
  X2 <- X^2
  X3 <- X^3
  X4 <- X^4
  X5 <- X^5
  X6 <- X[,1]*X[,2]
  X7 <- X[,1]^2*X[,2]
  X8 <- X[,1]*X[,2]^2
##   X2D <- cbind(X,X2,X3,X4,X5,X6,X7,X8)
  X2D <- cbind(X,X6)
  for (i in ncol(X2D)) {
    X2D[,i] <- X2D[,i] - mean(X2D[,i])
    X2D[,i] <- X2D[,i] / sd(X2D[,i])
  }

  Xr <- as.matrix(expand.grid(1:pdim[1],1:pdim[2]))[,1]
  Xr <- Xr %% 2
  X[,1] <- 0
  X[Xr==0,1] = X[Xr==0,2]
  X[Xr==0,2] = 0
  X2[,1] <- 0
  X2[Xr==0,1] = X2[Xr==0,2]
  X2[Xr==0,2] = 0
  X3[,1] <- 0
  X3[Xr==0,1] = X3[Xr==0,2]
  X3[Xr==0,2] = 0
  X4[,1] <- 0
  X4[Xr==0,1] = X4[Xr==0,2]
  X4[Xr==0,2] = 0
  X5[,1] <- 0
  X5[Xr==0,1] = X5[Xr==0,2]
  X5[Xr==0,2] = 0
##   XMD <- cbind(X,X2,X3,X4,X5)
  XMD <- cbind(X,X2,X3)

  for (i in ncol(XMD)) {
    XMD[,i] <- XMD[,i] - mean(XMD[,i])
    XMD[,i] <- XMD[,i] / sd(XMD[,i])
  }

  P <- max(Plate)
  
  plateeffect <- D
  plateeffect[] <- 0


  posn <- 1:prod(pdim)
  row <- 1 + (posn - 1)%%pdim[1]
  col <- 1 + (posn - 1)%/%pdim[1]
  centre <- 0.5 + c(pdim[1]/2, pdim[2]/2)
  xpos <- col - centre[2]
  ypos <- centre[1] - row

  
  if (type == "Bscore") {
    ZZ <- sgi@S * sgi@C * P
  } else {
    ZZ <- sgi@S * sgi@C * P * maxit
  }
  zz <- 0
  for (s in getScreenNames(sgi)) {
    for (c in getChannelNames(sgi)) {
      if (type == "Bscore") {
        for (p in 1:P) {
          zz <- zz + 1
          if (verbose > 1) {
            message(sprintf("Normalization plate %d from %d.\r",zz, ZZ), appendLF=FALSE)
          }
          I <- Plate == p
          tmp[Well[I]] <- D[I,s,c]
          plate[] <- tmp
          M = medpolish(plate,na.rm=TRUE)
          PE <- matrix(M$row,nrow=pdim[1],ncol=pdim[2]) + t(matrix(M$col,nrow=pdim[2],ncol=pdim[1])) + M$overall
          plate <- plate - PE
          tmp[] <- plate[]
          D[I,s,c] <- tmp[Well[I]]
          tmp[] <- PE[]
          plateeffect[I,s,c] <- tmp[Well[I]]
        }
      }

      if (type == "spatial") {
        for (p in 1:P) {
          zz <- zz + 1
          if (verbose > 1) {
            message(sprintf("Normalization plate %d from %d.\r",zz, ZZ), appendLF=FALSE)
          }
          I <- Plate == p
          tmp[Well[I]] <- D[I,s,c]
          df = data.frame(y = tmp, xpos = xpos, ypos = ypos)
          m = locfit(y ~ lp(xpos, ypos), data = df, lfproc = locfit.robust)
          predx = predict(m, newdata = df)
          tmp[] <- df$y - predx
          D[I,s,c] <- tmp[Well[I]]
          tmp[] <- predx
          plateeffect[I,s,c] <- tmp[Well[I]]
        }
      }
      
    }
  }
  sgi@pi <- D
  sgi@plateeffect <- plateeffect

  if (verbose>0) {
    time2 = Sys.time()
    message(sprintf("Normalize plate effects. Done. Time: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }

  return(sgi)
}







