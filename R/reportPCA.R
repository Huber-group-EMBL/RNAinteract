
embedPCA <- function(sgi, screen, channel, dim=4, embed = "template", withoutgroups=c()) {
  X <- getData(sgi, type="pi", format="targetMatrix", screen=screen, channel=channel, withoutgroups=withoutgroups)
  if (any(!is.finite(X))) {
    if (any(is.finite(X))) {
      M = mean(X[is.finite(X)])
    } else {
      M = 0
    }
    X[!is.finite(X)] = M
  }
  pc <- prcomp(X)
  PC <- X %*% pc$rotation
  PC <- PC[,1:dim]
  return (PC)
}

#################


## linearEmbedding <- function(X, y = NULL, type="PCA", stop=-2, dim=4, symm = TRUE) {
##   model = list()

##   if (any(!is.finite(X))) {
##     if (any(is.finite(X))) {
##       M = mean(X[is.finite(X)])
##     } else {
##       M = 0
##     }
##     X[!is.finite(X)] = M
##   }
##   model$type = type
##   model$X = X
##   if (type == "PCA") {
##     model$main = "Principal Component Analysis"
##     model$pc <- prcomp(X)
##     PC <- X %*% model$pc$rotation
##     model$PC <- PC[,1:dim]
##   }
##   if (type == "LDA") {
##     if (!is.factor(y)) { y = factor(y) }
##     model$main = "Sparse Linear Discriminant Analysis"
##     model$y = y
##     model$sdamodel <- sda(X,y,stop=stop)
##     train <- predict(model$sdamodel, X)
##     model$PC <- train$x
##     model$selectedFeatures <- model$sdamodel$varIndex
##   }

##   if (dim(X)[1] != dim(X)[2]) {
##     symm = FALSE
##   }

##   model$symm = symm
##   model$hc.row = hclust(dist(model$PC))
##   if (symm) {
##     model$hc.col = model$hc.row
##   } else {
##     model$hc.col = hclust(dist(t(X)))
##   }
##   return(model)
## }

## getConfusionMatrix <- function(X, y = NULL, type="LDA", stop=-2) {

##   if (any(!is.finite(X))) {
##     if (any(is.finite(X))) {
##       M = mean(X[is.finite(X)])
##     } else {
##       M = 0
##     }
##     X[!is.finite(X)] = M
##   }

##   y = factor(y)

##   df = NULL
##   if (type == "LDA") {
##     ypred = y
##     ypred[] = NA
##     posterior = matrix(NA,nr=length(ypred),nc=4)
##     row.names(posterior) = row.names(X)
##     for (i in 1:length(y)) {
##       model <- sda(X[-i,,drop=FALSE],y[-i],stop=stop)
##       x <- predict(model, X)$x[,1:2]
##       model2 = lda(x[-i,,drop=FALSE],y[-i])
##       pp = predict(model2,x[i,,drop=FALSE])
##       ypred[i] <- pp$class
##       posterior[i,] <- pp$posterior
##       colnames(posterior) = colnames(pp$posterior)
##     }
##     df = data.frame(true = y, predicted = ypred, row.names=row.names(X))
##     save(posterior,y,file="posterior.rdb")
##   }
  
##   return(df)
## }

## linearEmbeddingNew <- function(model, Xnew, symm = TRUE) {

##   if (any(!is.finite(Xnew))) {
##     if (any(is.finite(model$X))) {
##       M = mean(model$X[is.finite(model$X)])
##     } else {
##       M = 0
##     }
##     Xnew[!is.finite(Xnew)] = M
##   }

##   if (model$type == "PCA") {
##     PC <- Xnew %*% model$pc$rotation
##     model$PCnew <- PC[,1:dim(model$PC)[2]]
##   }
##   if (model$type == "LDA") {
##     test <- predict(model$sdamodel, Xnew)
##     model$PCnew = test$x
##     model$classnew = test$class
##     model$posterior = test$posterior
##   }

##   if (dim(Xnew)[1] != dim(Xnew)[2]) {
##     symm = FALSE
##   }
##   model$symm.new = symm
##   model$hc.new.row = hclust(dist(model$PCnew))
##   if (symm) {
##     model$hc.new.col = model$hc.new.row
##   } else {
##     model$hc.new.col = hclust(dist(t(Xnew)))
##   }
##   return(model)
## }

## plotEmbedding <- function(PC, d=c(1,2), pch=19, col=NULL, main="Principal Component Analysis",
##                           text = rep("", nrow(PC)), legend.args = NULL,
##                           cex.points=1.4, cex.text=1.0) {
##   if (is.null(col)) {
##     col = rep("black", nrow(PC))
##   }
##   plot(PC[,d[1]],PC[,d[2]],pch=19,col=par("bg"),cex=0, main=main, xlab=sprintf("Dimension %d",d[1]), ylab = sprintf("Dimension %d",d[2]), xaxt="n", yaxt="n",bty="n")

##   # plot points
##   if ("white" %in% col) {
##     I = which(col == "white")
##     points(PC[I,d[1]],PC[I,d[2]],pch=21,col="gray50", bg=col[I],cex=cex.points)
##   }
##   if ("gray" %in% col) {
##     I = which(col == "gray")
##     points(PC[I,d[1]],PC[I,d[2]],pch=21,col="gray50", bg=col[I],cex=cex.points)
##   }
##   if ("black" %in% col) {
##     I = which(col == "black")
##     points(PC[I,d[1]],PC[I,d[2]],pch=21,col="gray50", bg=col[I],cex=cex.points)
##   }
##   I = which(!(col %in% c("white","gray","black")))
##   if (length(I) > 0) {
##     points(PC[I,d[1]],PC[I,d[2]],pch=21,col="gray50", bg=col[I],cex=cex.points)
##   }

##   # plot text
##   if ("gray" %in% col) {
##     I = which(col == "gray")
##     text(PC[I,d[1]],PC[I,d[2]],text[I],col=col[I],cex=cex.text)
##   }
##   if ("black" %in% col) {
##     I = which(col == "black")
##     text(PC[I,d[1]],PC[I,d[2]],text[I],col=col[I],cex=cex.text)
##   }
##   I = which(!(col %in% c("white","gray","black")))
##   if (length(I) > 0) {
##     text(PC[I,d[1]],PC[I,d[2]],text[I],col=col[I],cex=cex.text)
##   }

##   # legend
##   if (!is.null(legend.args)) {
##     do.call("legend", legend.args)
##   }
## }

## plotEmbeddingLDA <- function(model, pch=NULL, col=NULL, groupnew=NULL, groupcol=groupcol, main="Linear Discriminant Analysis",
##                              text = rep("", nrow(PC)), legend.args = NULL,
##                              cex.points=1.4, cex.text=1.0) {
##   PC <- model$PCnew
##   d <- c(1,2)

##   if (is.null(col)) {
##     col = rep("black", nrow(PC))
##   }
##   if (is.null(pch)) {
##     pch = rep(1, nrow(PC))
##   }

##   ldamodel <- lda(model$PC[,1:2],model$y)
##   rx <- range(model$PCnew[,1])
##   ry <- range(model$PCnew[,2])
##   tmp <- 0.05*(rx[2] - rx[1])
##   rx <- c(rx[1]-tmp,rx[2]+tmp)
##   tmp <- 0.05*(ry[2] - ry[1])
##   ry <- c(ry[1]-tmp,ry[2]+tmp)
##   Xnew <- matrix(0,nr=160000,nc=2)
##   x <- seq(rx[1],rx[2],length.out=400)
##   y <- seq(ry[1],ry[2],length.out=400)
##   Xnew[,1] <- rep(x, each=400)
##   Xnew[,2] <- rep(rev(y),times=400)
##   pp = predict(ldamodel,Xnew)
##   cl <- as.character(pp$class)
##   C = t(col2rgb(groupcol))/255
##   C = C[match(colnames(pp$posterior), row.names(C)),]
## ##   C[C >= 1] = 1.0
## ##   C[C <= 0.0] = 0.0
##   C = 1-0.2*(1-C)
##   M <- array(0.0, dim=c(400,400,3))
##   for (i in 1:3) {
##     M[,,i] <- pp$posterior %*% C[,i,drop=FALSE]
##   }
##   M[M > 0.9999] = 0.9999
## #  save(model,pch,col,groupnew,groupcol,file="~/tmp.rdb")
##   dx = 0.5*(x[2] - x[1])
##   dy = 0.5*(y[2] - y[1])
##   vp = dataViewport(width=1, height=1,
##               xscale = c(x[1]-dx, x[400]+dx),
##               yscale = c(y[1]-dy, y[400]+dy),
##               clip="off")
##   pushViewport(vp)
##   vpimg = dataViewport(width = 1, height = 1, xscale = c(0.5, 400.5), yscale = c(0.5,400.5), clip = "off")
##   pushViewport(vpimg)
##   grid.raster(image = as.raster(M),
##                 y = 0.5, height = 400.5, hjust=0, 
##                 x = 0.5,  width  = 400.5, vjust=0,
##                 default.units = "native", interpolate=FALSE)
##   popViewport()


## ##  save(ldamodel,pp,file="tmplda.rdb")
## ##   zlim=c(1,length(groupcol))
## ##   cc = rgb(t(col2rgb(groupcol))/255,alpha=0.2)
## ##   save(x,y,M,zlim,cc,PC,col,file="tmplp.rdb")

## ##   M <- matrix(match(cl, names(groupcol)),nr=400,nc=400)
## ##   image(x,y,M,zlim=c(1,length(groupcol)),col=rgb(t(col2rgb(groupcol))/255,alpha=0.2), main=main)

##   # plot points
##   if ("white" %in% col) {
##     I = which(col == "white")
##     grid.points(PC[I,d[1]],PC[I,d[2]],pch=pch[I],gp=gpar(col="gray50", fill=col[I]))
## #    points(PC[I,d[1]],PC[I,d[2]],pch=pch[I],col="gray50", bg=col[I],cex=cex.points)
##   }
##   if ("gray" %in% col) {
##     I = which(col == "gray")
##     grid.points(PC[I,d[1]],PC[I,d[2]],pch=pch[I],gp=gpar(col="gray50", fill=col[I]))
## #    points(PC[I,d[1]],PC[I,d[2]],pch=pch[I],col="gray50", bg=col[I],cex=cex.points)
##   }
##   if ("black" %in% col) {
##     I = which(col == "black")
##     grid.points(PC[I,d[1]],PC[I,d[2]],pch=pch[I],gp=gpar(col="gray50", fill=col[I]))
## #    points(PC[I,d[1]],PC[I,d[2]],pch=pch[I],col="gray50", bg=col[I],cex=cex.points)
##   }
##   I = which(!(col %in% c("white","gray","black")))
##   if (length(I) > 0) {
##     grid.points(unit(PC[I,d[1]],units="native"),unit(PC[I,d[2]],units="native"),pch=pch[I],gp=gpar(col="gray50", fill=col[I]))
## #    points(PC[I,d[1]],PC[I,d[2]],pch=pch[I],col="gray50", bg=col[I],cex=cex.points)
##   }

##   if (!is.null(groupnew)) {
##     cl <- as.character(predict(ldamodel,model$PCnew[,1:2])$class)
##     I <- ((cl != groupnew[row.names(model$PCnew)]) | ((groupnew[row.names(model$PCnew)] != "others") & (pch == 21)))
##     grid.text(row.names(model$PCnew)[I],unit(PC[I,d[1]],units="native"),unit(PC[I,d[2]],units="native"),just=c("left","bottom"),gp=gpar(col="black"))
## #    text(PC[I,d[1]],PC[I,d[2]],row.names(model$PCnew)[I],adj=c(0,0),col="black",cex=cex.text)
##   }

##   popViewport()

##   # legend
## ##   if (!is.null(legend.args)) {
## ##     do.call("legend", legend.args)
## ##   }
## }

## plotHeatmap <- function(X, hc.row = NULL, hc.col = NULL, symm = TRUE, pi.min=NULL, pi.max=NULL, cex=1.0, xlab="features",col.row=NULL, col.col=NULL) {

##   if (any(!is.finite(X))) {
##     if (any(is.finite(X))) {
##       M = mean(X[is.finite(X)])
##     } else {
##       M = 0
##     }
##     X[!is.finite(X)] = M
##   }
##   if (is.null(hc.row)) {
##     hc.row <- hclust(dist(X))
##   }
##   if (is.null(hc.col)) {
##     hc.col <- hclust(dist(t(X)))
##   }
##   dd.row = as.dendrogram(hc.row)
##   dd.col = as.dendrogram(hc.col)
##   ord.row = order.dendrogram(dd.row)
##   ord.col = order.dendrogram(dd.col)

## ##   rbcol = rev(colorRampPalette(brewer.pal(9,"RdBu"))(513))
##   rbcol <- colorRampPalette(c("cornflowerblue","black","yellow"))(513)

##   if (is.null(pi.min)) {
##     pi.min = mad(X,center=0.0)
##   }
##   if (is.null(pi.max)) {
##     pi.max = 3*pi.min
##   }
##   at = c(seq(-pi.max,-pi.min,length.out=257),seq(pi.min,pi.max,length.out=257))
  
##   X = t(X[ord.row, ord.col])
##   X[X > pi.max] = pi.max
##   X[X < -pi.max] = -pi.max

##   if (is.null(col.row)) {
##     add.row = list()
##   } else {
##     add.row <- list(rect=list(fill=col.row))
##   }
##   if (is.null(col.col)) {
##     add.col = list()
##   } else {
##     add.col <- list(rect=list(fill=col.col))
##   }
##   print(levelplot(X,
##                   aspect = "iso",
##                   scales = list(x = list(rot = 90), cex = cex),
##                   col.regions=rbcol,
##                   ylab = expression(pi-score),
##                   xlab = xlab,
##                   at=at,
##                   colorkey = list(space = "left"),
##                   legend =
##                   list(right =
##                        list(fun = dendrogramGrob,
##                             args =
##                             list(x = dd.row, ord = ord.row,
##                                  side = "right",
##                                  size = 4,
##                                  add = add.row)),
##                        top =
##                        list(fun = dendrogramGrob,
##                             args =
##                             list(x = dd.col, ord = ord.col, 
##                                  side = "top",
##                                  size = 4,
##                                  add = add.col)))))
## }

## reportLinearEmbedding <- function(model, X, symm = TRUE, path=".", dir="PCA", prefix = "pca",
##                                   confusion = NULL, width = 480, height = 480,
##                                   col=NULL, pch=NULL, text = row.names(X),
##                                   cex.points = 1.4, cex.text=0.0001, cex.text.lda=1.0, cex.heatmap = 0.3, groupcol = NULL, groupnew = NULL,
##                                   legend.args = list(x="topright",legend=c("RasMapK", "RasMapK-Inh", "JNK","others"), inset=0.02, fill=c("red","orange","blue","white"))) {

##   dir.create(sprintf("%s/%s", path, dir), showWarnings = FALSE)

## ##   I <- order(col)
## ##   X <- X[I,I]
## ##   col <- col[I]
## ##   text <- text[I]
##   model <- linearEmbeddingNew(model, X, symm)
## posterior = model$posterior

##   names <- row.names(X)
## ##   if (is.null(col)) {
## ##     col <- rep("black",nrow(X))
## ##   }
##   col.row = col
##   if (model$symm.new) {
##     col.col = col.row
##   } else {
##     col.col = NULL
##   }

##   md <- dim(model$PCnew)[2]
##   if (md > 3) { md = 3 }
##   z <- ifelse(is.null(model$selectedFeatures), 1, 2)
##   if (!is.null(groupcol)) { z <- z+1 }
##   fn.png <- matrix("", nr=md*(md-1)/2+z, nc=1)
##   fn.pdf <- matrix("", nr=md*(md-1)/2+z, nc=1)
##   rn <- rep("", md*(md-1)/2)
##   zz <- 0
##   if (!is.null(model$selectedFeatures)) {
##     png(width=width,height=height,file=sprintf("%s/%s/%s-heatmap-sel.png",path, dir, prefix))
##     plotHeatmap(X[,model$selectedFeatures], hc.row = model$hc.new.row, col.row=col.row, cex=cex.heatmap, xlab="selected features")
##     dev.off()
##     pdf(file=sprintf("%s/%s/%s-heatmap-sel.pdf",path, dir, prefix))
##     plotHeatmap(X[,model$selectedFeatures], hc.row = model$hc.new.row, col.row=col.row, cex=cex.heatmap, xlab="selected features")
##     dev.off()
##     zz <- zz + 1
##     fn.png[zz,1] <- sprintf("%s-heatmap-sel.png",prefix)
##     fn.pdf[zz,1] <- sprintf("%s-heatmap-sel.pdf",prefix)
##     rn[zz] <- ""
##   }
##   png(width=width,height=height,file=sprintf("%s/%s/%s-heatmap.png",path, dir, prefix))
##   plotHeatmap(X, hc.row = model$hc.new.row, hc.col=model$hc.new.col, col.row=col.row, col.col=col.col, cex=cex.heatmap)
##   dev.off()
##   pdf(file=sprintf("%s/%s/%s-heatmap.pdf",path, dir, prefix))
##   plotHeatmap(X, hc.row = model$hc.new.row, hc.col=model$hc.new.col, col.row=col.row, col.col=col.col, cex=cex.heatmap)
##   dev.off()
##   zz <- zz + 1
##   fn.png[zz,1] <- sprintf("%s-heatmap.png",prefix)
##   fn.pdf[zz,1] <- sprintf("%s-heatmap.pdf",prefix)
##   rn[zz] <- ""
##   if (!is.null(groupcol)) {
##     zz <- zz + 1
##     png(width=width,height=height,file=sprintf("%s/%s/%s-phenomap.png",path, dir, prefix))
##     plotEmbeddingLDA(model, text = text, col=col, pch=pch, main=model$main, groupcol=groupcol,groupnew=groupnew,
##                      legend.args = legend.args, cex.points=cex.points, cex.text=cex.text.lda)
##     dev.off()
##     pdf(file=sprintf("%s/%s/%s-phenomap.pdf",path, dir, prefix))
##     plotEmbeddingLDA(model, text = text, col=col, pch=pch, main=model$main, groupcol=groupcol,groupnew=groupnew,
##                      legend.args = legend.args, cex.points=cex.points, cex.text=cex.text.lda)
##     dev.off()
##     fn.png[zz,1] <- sprintf("%s-phenomap.png",prefix)
##     fn.pdf[zz,1] <- sprintf("%s-phenomap.pdf",prefix)
##     rn[zz] <- "phenomap"
##   }
##   if (dim(model$PCnew)[2] > 1) {
##     for (i in 1:(md-1)) {
##       for (j in (i+1):md) {
##         zz <- zz + 1
##         png(width=width,height=height,file=sprintf("%s/%s/%s-%d-%d.png",path, dir, prefix, i, j))
##         plotEmbedding(model$PCnew, d=c(i,j), text = text, col=col, pch=pch, main=model$main,
##                       legend.args = legend.args, cex.points=cex.points, cex.text=cex.text)
##         dev.off()
##         pdf(file=sprintf("%s/%s/%s-%d-%d.pdf",path, dir, prefix, i, j))
##         plotEmbedding(model$PCnew, d=c(i,j), text = text, col=col, pch=pch, main=model$main,
##                       legend.args = legend.args, cex.points=cex.points, cex.text=cex.text)
##         dev.off()
##         fn.png[zz,1] <- sprintf("%s-%d-%d.png",prefix, i, j)
##         fn.pdf[zz,1] <- sprintf("%s-%d-%d.pdf",prefix, i, j)
##         rn[zz] <- sprintf("%d - %d", i, j)
##       }
##     }
##   }

##   p = openPage(sprintf("%s/%s/index-%s.html",path, dir, prefix))

##   hwrite(sprintf("<h3>%s</h3>",model$main),table=FALSE,br=FALSE,page=p)
##   Img <- hwriteImage(fn.png,table=FALSE)
##   row.names(Img) <- row.names(fn.png)
##   hwrite(Img, link=fn.pdf, cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

##   if (!is.null(confusion)) {
##     hwrite("<h3>Cross Validation of Classifier</h3>",table=FALSE,br=TRUE,page=p)
##     C = as.matrix(table(confusion))
##     class(C) = "matrix"
##     hwrite(C,table=TRUE,br=TRUE,page=p)
##     png(width=500,height=500,file=sprintf("%s/%s/%s-confusion.png",path, dir, prefix))
##     balloonplot(table(confusion))
##     dev.off()
##     pdf(file=sprintf("%s/%s/%s-confusion.pdf",path, dir, prefix))
##     balloonplot(table(confusion))
##     dev.off()
##     hwriteImage(sprintf("%s-confusion.png",prefix),link=sprintf("%s-confusion.pdf",prefix), page=p)
##     I = order(as.integer(confusion[,1])*100+as.integer(confusion[,2]))
##     confusion = confusion[I,]
##     hwrite(as.matrix(confusion),br=TRUE,page=p)
##   }

##   closePage(p,splash=FALSE)

##   invisible(model)
## }

## reportPCA <- function(sgi, verbose = TRUE, path = ".", dir = "PCA", prefix = "pca", page = NULL, groups = NULL, groupcol = NULL, withoutgroups=c("neg", "pos")) {
##   stopifnot( is( sgi, "RNAinteract" ) )

##   if (verbose) {
##     time1 <- Sys.time()
##     cat("Write PCA. \r")
##   }
##   dir.create(sprintf("%s/%s", path, dir), showWarnings = FALSE)

##   res <- list()
##   if (!is.null(page)) {
##     hwrite(sprintf("<b>%s</b>",prefix), br=TRUE, page=page)
##     for (c in getChannelNames(sgi)) {
##       hwrite(paste(c," ",collapse=""), br=FALSE, page=page)
##       for (s in getScreenNames(sgi)) {
##         hwrite(sprintf("[%s]",s), link = sprintf("%s/index-%s-%s-%s.html",dir, prefix, s, c), target = "showframe", br=FALSE, page=page)
##       }
##       hwrite("<br>", br=FALSE, page=page)
##     }
##     hwrite("<br>", br=FALSE, page=page)
##   }

##   zz <- 0
##   ZZ <- sgi@S * sgi@C
##   for (s in getScreenNames(sgi)) {
##     res[[s]] <- list()
##     for (c in getChannelNames(sgi)) {
##       if (verbose) {
##         zz <- zz + 1
##         cat("Write PCA. Nr ", zz, " from ", ZZ,".\r")
##       }

##       D <- getData(sgi, type="pi", format="targetMatrix",screen=s,channel=c, withoutgroups=withoutgroups)
##       col <- rep("white", nrow(D))
##       if (!is.null(groups)) {
##         if (is.null(groupcol)) {
##           groupcol <- rainbow(length(groups))
##         }
##         for (i in 1:length(groups)) {
##           col[row.names(D) %in% groups[[i]]] <- groupcol[i]
##         }
##       }
##       model <- linearEmbedding(X=D, type="PCA", dim=4)
##       res[[s]][[c]] <- reportLinearEmbedding(model, D, width=680, height=680, cex.heatmap=0.3, path=path, dir=dir, prefix=sprintf("%s-%s-%s",prefix,s,c), col=col)
##     }
##   }

##   if (verbose) {
##     time2 = Sys.time()
##     cat("Write PCA. Done. Time elapsed: ",format(difftime(time2, time1, tz = "",units = c("auto"))),".                                      \r\n")
##   }
##   invisible(res)
## }

## reportLDA <- function(sgi, traingroups, verbose = TRUE, path = ".", dir = "LDA", prefix = "lda", page = NULL, groups = NULL, groupcol = NULL, do.cv = FALSE, withoutgroups=c("neg", "pos")) {
##   stopifnot( is( sgi, "RNAinteract" ) )

##   if (verbose) {
##     time1 <- Sys.time()
##     cat("Write LDA. \r")
##   }
##   dir.create(sprintf("%s/%s", path, dir), showWarnings = FALSE)

##   res <- list()
##   if (!is.null(page)) {
##     hwrite(sprintf("<b>%s</b>",prefix), br=TRUE, page=page)
##     for (c in getChannelNames(sgi)) {
##       hwrite(paste(c," ",collapse=""), br=FALSE, page=page)
##       for (s in getScreenNames(sgi)) {
##         hwrite(sprintf("[%s]",s), link = sprintf("%s/index-%s-%s-%s.html",dir, prefix, s, c), target = "showframe", br=FALSE, page=page)
##       }
##       hwrite("<br>", br=FALSE, page=page)
##     }
##     hwrite("<br>", br=FALSE, page=page)
##   }

##   zz <- 0
##   ZZ <- sgi@S * sgi@C
##   for (s in getScreenNames(sgi)) {
## #  for (s in "mean") {
##     res[[s]] <- list()
##     for (c in getChannelNames(sgi)) {
## #    for (c in "nrcells") {
##       if (verbose) {
##         zz <- zz + 1
##         cat("Write LDA. Nr ", zz, " from ", ZZ,".\r")
##       }

##       D <- getData(sgi, type="pi", format="targetMatrix",screen=s,channel=c, withoutgroups=withoutgroups)
##       pch <- rep(21, nrow(D))
##       col <- rep("white", nrow(D))
##       groupnew <- NULL
##       if (!is.null(groups)) {
##         if (is.null(groupcol)) {
##           groupcol <- rainbow(length(groups))
##         }
##         for (i in 1:length(groups)) {
##           col[row.names(D) %in% groups[[i]]] <- groupcol[i]
##         }
##         groupnew <- rep(NA, nrow(D))
##         names(groupnew) <- row.names(D)
##         for (i in 1:length(groups)) {
##           groupnew[row.names(D) %in% groups[[i]]] <- names(groups)[i]
##         }
##       }
##       pch[row.names(D) %in% unlist(traingroups)] <- 24
##       y <- rep(NA, nrow(D))
##       for (i in 1:length(traingroups)) {
##         y[row.names(D) %in% traingroups[[i]]] <- names(traingroups)[i]
##       }
##       Dtrain = D[!is.na(y),]
##       y = y[!is.na(y)]
##       model <- linearEmbedding(X=Dtrain, y=y, type="LDA", stop=-2)
##       confusion <- NULL
##       if (do.cv) {
##         confusion <- getConfusionMatrix(X=Dtrain, y=y, type="LDA", stop=-2)
## ##         file.rename("posterior.rdb",sprintf("posterior_%s_%s.rdb",s,c))
##       }
##       res[[s]][[c]] <- reportLinearEmbedding(model, D, confusion = confusion, width=680, height=680, cex.heatmap=0.3, path=path, dir=dir, prefix=sprintf("%s-%s-%s",prefix,s,c), col=col, pch=pch, groupcol=groupcol,groupnew=groupnew)
##       file.rename("posteriornew.rdb",sprintf("posteriornew_%s_%s.rdb",s,c))
##       res[[s]][[c]]$col <- col
##     }
##   }

##   if (verbose) {
##     time2 = Sys.time()
##     cat("Write LDA. Done. Time elapsed: ",format(difftime(time2, time1, tz = "",units = c("auto"))),".                                      \r\n")
##   }
##   invisible(res)
## }


