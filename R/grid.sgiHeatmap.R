
grid.sgiDendrogram <- function(hc = NULL, side = "top", vp=NULL) {
  dd.col = as.dendrogram(hc)
  ord.col = order.dendrogram(dd.col)
  dg <- dendrogramGrob(x = dd.col, ord = ord.col, side = side, size.add=0.3, size = 2, vp=vp)
  grid.draw(dg)
}


grid.sgiColorkey <- function(X, pi.min=NULL, pi.max=NULL, height=NULL,main="") {
  rbcol <- colorRampPalette(c("cornflowerblue","black","yellow"))(513)
  C = t(col2rgb(rbcol))/255

  if (is.null(pi.min)) {
    pi.min = mad(X,center=0.0)
  }
  if (is.null(pi.max)) {
    pi.max = 3*pi.min
  }
  pi.max2 = pi.max
  pi.max2 = -pretty(c(-pi.max, pi.max),min.n=2)[1]
  pi = abs(X)
  pi = pi[(pi >= pi.min) & (pi <= pi.max)]
  quant = c(quantile(pi, probs = seq(0,1,length.out=256)),max(pi,na.rm=TRUE))
  at = c(-quant,quant)
  zrange = c(-pi.max2,pi.max2)
  ## yPos = pretty(zrange)
  ## yPos = yPos[2:(length(yPos)-1)]
  ## yText = paste(yPos)
  steps = seq(-pi.max2, pi.max2, length=256)
  colid = sapply(steps, function(e) { sum(e > at) } )
  colid[colid < 1] = 1
  colid[colid > nrow(C)] = nrow(C)
  colid = rev(colid)
  Img = array(0, dim=c(length(steps),1,3))
  Img[,1,1] = C[colid,1]
  Img[,1,2] = C[colid,2]
  Img[,1,3] = C[colid,3]
  if (is.null(height)) { height = unit(1, "npc") }
  pushViewport(viewport(height=height,layout=grid.layout(nrow = 1, ncol = 4,heights=unit(20,"lines"),
                                        widths = unit.c(unit(1.3*1.8,"strheight","cell number"), unit(0.2,"strwidth","1.555"), unit(0.6,"strwidth","1.555"), unit(0.2,"strwidth","1.555") ))))
  pushViewport(viewport(layout.pos.col=3))
  grid.raster(Img,
              interpolate = FALSE,
              x = unit(0.5,"npc"), width  = unit(1, "npc"), 
              y = unit(0.5, "npc"), height = unit(1,"npc")-unit(3.2,"lines"),
              just=c("center", "center"))
  grid.rect(x = unit(0.5,"npc"), width  = unit(1, "npc"), 
              y = unit(0.5, "npc"), height = unit(1,"npc")-unit(3.2,"lines"),
              just=c("center", "center"))
  grid.text(label = -pi.max2,
            x = unit(0.5,"npc"),
            y = unit(0, "npc"),
            default.units = "native", just=c("center", "bottom"),gp=gpar(cex=1.8))
  grid.text(label = pi.max2,
            x = unit(0.5,"npc"),
            y = unit(1, "npc"),
            default.units = "native", just=c("center", "top"),gp=gpar(cex=1.8))
  popViewport()
  pushViewport(viewport(layout.pos.col=1))
  grid.text(label = main,
            x = unit(0,"lines"),
            y = unit(0.5,"npc"),
            rot=90,
            gp=gpar(cex=1.8),
            default.units = "native", just=c("center", "top"))
  popViewport()
  popViewport()
  h = convertHeight(unit(1,"npc"),"cm")
  h
}

grid.sgiHeatmap.2 <- function(X, hc.row = NULL, hc.col = NULL, symm = TRUE, pi.min=NULL, pi.max=NULL, cex=1.0, xlab="features", col.row=NULL, col.col=NULL) {

  if (is.null(hc.row)) {
    hc.row <- hclust(dist(X))
  }
  if (is.null(hc.col)) {
    hc.col <- hclust(dist(t(X)))
  }
  dd.row = as.dendrogram(hc.row)
  dd.col = as.dendrogram(hc.col)
  ord.row = order.dendrogram(dd.row)
  ord.col = order.dendrogram(dd.col)

  rbcol <- colorRampPalette(c("cornflowerblue","black","yellow"))(513)

  if (is.null(pi.min)) {
    pi.min = mad(X,center=0.0)
  }
  if (is.null(pi.max)) {
    pi.max = 3*pi.min
  }
  pi.max2 = pi.max
  pi.max2 = -pretty(c(-pi.max, pi.max),min.n=2)[1]
  pi = abs(X)
  pi = pi[(pi >= pi.min) & (pi <= pi.max)]
  quant = c(quantile(pi, probs = seq(0,1,length.out=256)),pi.max2)
  at = c(-quant,quant)
  zrange = c(-pi.max2,pi.max2)
  ## yPos = pretty(zrange)
  ## yPos = yPos[2:(length(yPos)-1)]
  ## yText = paste(yPos)
  steps = seq(-pi.max2, pi.max2, length=256)
  colid = sapply(steps, function(e) { sum(e > at) } )
  colid[colid < 1] = 1
  colid[colid > nrow(C)] = nrow(C)
  colid = rev(colid)
  
  ## pi = abs(upperTriangle(X))
  ## pi = X
  ## pi = pi[(pi >= pi.min)]
  ## quant = quantile(pi, probs = seq(0,1,length.out=257))
  ## at = c(-quant,quant)
  ## X = t(X[ord.row, ord.col])
  X = t(X[ord.row, ord.col])
  X[X > pi.max] = pi.max
  X[X < -pi.max] = -pi.max
  if (is.null(col.row)) {
    add.row = list()
  } else {
    add.row <- list(rect=list(fill=col.row))
  }
  if (is.null(col.col)) {
    add.col = list()
  } else {
    add.col <- list(rect=list(fill=col.col))
  }
  
  ## diag(X) = 0.0
  lp = levelplot(X,
                  aspect = "iso",
                  scales = list(draw=FALSE,alternating=0),
                  panel = function(...) {
                    panel.levelplot(...)
                    panel.axis(side = "right", at=1:ncol(X), labels = colnames(X),text.cex=cex,outside = TRUE)
                    panel.axis(side = "bottom", at=1:nrow(X), labels = row.names(X),text.cex=cex,outside = TRUE)
                  },
                  par.settings = list(clip=list(panel="off",strip="off")),
                  col.regions=rbcol,
                  ylab = NULL,
                  xlab = NULL,
                  at=at,
                  colorkey=NULL,
    )
  return(lp)
}

grid.sgiHeatmap <- function(PI, pi.max=NULL, main=expression(paste(pi,"-score")), hc.row = NULL, hc.col = NULL) {
  vp=viewport(name="hm", layout = grid.layout(3,4,
                           widths = unit.c(unit(1.3*1.8,"strheight","cell number"), unit(1,"strwidth","1.555"),unit(1,"null"), unit(1,"strwidth","abcdefghi")),
                           heights = unit.c(unit(0.1,"lines"), unit(1,"null"),unit(0.1,"lines")),
                           respect=TRUE)
    )
  pushViewport(vp)
  pushViewport(viewport(name="dhm", layout.pos.row=2, layout.pos.col=3,clip=FALSE))
  print(grid.sgiHeatmap.2(PI, hc.row = hc.row, hc.col = hc.col, pi.max=pi.max), newpage = FALSE)
  popViewport()
  pushViewport(viewport(name="dhm", layout.pos.row=2, layout.pos.col=2,clip=FALSE))
  grid.sgiColorkey(PI, main=main, pi.max=pi.max)
  popViewport()
  popViewport()
}

ordertree <- function(hc, level) {
  res <- c()
  if (level < 0) {
    res <- -level
  } else {
    res <- c(ordertree(hc, hc$merge[level,1]), ordertree(hc, hc$merge[level,2]))
  }
  res
}

swaptree <- function(hc, level) {
  tmp <- hc$merge[level,1]
  hc$merge[level,1] <- hc$merge[level,2]
  hc$merge[level,2] <- tmp
  hc$order <- ordertree(hc, nrow(hc$merge))
  hc
}



