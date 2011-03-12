
getDoublePerturbationLabels <- function(mainEffect, dpEffect, mainTarget, text) {
  I = which((nchar(text) > 0) & (dpEffect >= mainEffect+mainTarget))
  if (length(I) > 0) {
    x = mainEffect[I]
    y = dpEffect[I]
    label = text[I]
    name = sprintf("labeltop%d",1:length(I))
    toplabels = data.frame(x=x,y=y,label=label,name=name,stringsAsFactors=FALSE)
    I = order(y-x)
    toplabels = toplabels[I,]
  } else {
    toplabels = data.frame()
  }
  I = which((nchar(text) > 0) & (dpEffect < mainEffect+mainTarget))
  if (length(I) > 0) {
    x = mainEffect[I]
    y = dpEffect[I]
    label = text[I]
    name = sprintf("labelbottom%d",1:length(I))
    bottomlabels = data.frame(x=x,y=y,label=label,name=name,stringsAsFactors=FALSE)
    I = order(x-y)
    bottomlabels = bottomlabels[I,]
  } else {
    bottomlabels = data.frame()
  }
  list(top=toplabels,bottom=bottomlabels)
}

makeDoublePerturbationLabels <- function(labels, gpText=NULL) {
  if (nrow(labels$top) > 0) {
    glabeltop = textGrob(labels$top$label,
      unit(labels$top$x,"native"),
      unit(labels$top$y,"native"),
      name="texttop",
      gp=gpar(col="red"),
      vp=vpPath("dplayout","dpvpclip"),
      just=c("right","bottom")
      )
  } else { glabeltop = NULL }
  if (nrow(labels$bottom) > 0) {
    glabelbottom = textGrob(labels$bottom$label,
      unit(labels$bottom$x,"native"),
      unit(labels$bottom$y,"native"),
      name="textbottom",
      gp=gpar(col="blue",gpText),
      vp=vpPath("dplayout","dpvpclip"),
      just=c("left","top")
      )
  } else { glabelbottom = NULL }
  gList(glabeltop, glabelbottom)
}

makeDoublePerturbationPoints <- function(mainEffect, dpEffect, pch=21, size=unit(1, "char"), fill=NULL, gpPoints=NULL) {
  if (!is.null(fill)) {
    if ("col" %in% names(fill)) {
      gpPoints = gpar(fill = fill$col, gpPoints)
    } else {
      if ("values" %in% names(fill)) {
        if ("at" %in% names(fill)) {
          at = fill$at
        } else {
          range = range(fill$values, finite=TRUE)
          at = seq(range[1], range[2], length.out=ifelse("colramp" %in% names(fill), length(fill$colramp)+1, 257))
        }
        if ("colramp" %in% names(fill)) {
          if (length(fill$colramp)+1 != length(fill$at)) {
            warning("length of fill$colramp has to be length of at - 1!")
          }
          colramp = fill$colramp
        } else {
          if ("colors" %in% names(fill)) {
            colors = fill$colors
          } else {
            colors = rainbow(32)
          }
          colramp = colorRampPalette(colors)(length(at)-1)
        }
        at = at[1:(length(at)-1)]
        colid = sapply(fill$values, function(e) { sum(e > at) } )
        colid[colid < 1] = 1
        colid[colid > length(colramp)] = length(colramp)
        col = colramp[colid]
        gpPoints = gpar(fill=col, gpPoints)
      } else {
        if (!("fill" %in% names(gpPoints))) {
          gpPoints = gpar(fill="lightblue", gpPoints)
        }
      }
    }
  } else {
    if (!("fill" %in% names(gpPoints))) {
      gpPoints = gpar(fill="lightblue", gpPoints)
    }
  }

  pointsGrob( mainEffect, dpEffect,
              pch=pch,
              size=size,
              name="points",
              gp=gpPoints,
              vp=vpPath("dplayout","dpvpclip") )
}

makeDoublePerturbationWTlines <- function(mainEffect, dpEffect, gpWTLines=NULL) {
  gList(
        linesGrob( c(0,1),
                  unit(c(0,0),"native"),
                  name="wtLineX",
                  gp=gpWTLines,
                  vp=vpPath("dplayout","dpvpclip")),
        linesGrob( unit(c(0,0),"native"),
                  c(0,1),
                  name="wtLineY",
                  gp=gpWTLines,
                  vp=vpPath("dplayout","dpvpclip"))
        )
}

makeDoublePerturbationMainEffects <- function(mainTarget, range, gpMain = NULL) {
  gList(
        linesGrob(c(0,1),
                  unit(c(mainTarget,mainTarget),"native"),
                  name="mainEffect",
                  gp=gpMain,
                  vp=vpPath("dplayout","dpvpclip")),
        linesGrob(unit(c(range[1],range[2]),"native"),
                  unit(c(range[1],range[2]),"native"),
                  name="mainEffectTarget",
                  gp=gpMain,
                  vp=vpPath("dplayout","dpvpclip"))
        )
}

makeDoublePerturbationNI <- function(mainTarget, range, gpNI = NULL) {
  linesGrob(unit(c(range[1],range[2]),"native"),
            unit(c(range[1]+mainTarget,range[2]+mainTarget),"native"),
            name="niLine",
            gp=gpNI,
            vp=vpPath("dplayout","dpvpclip"))
}

makeDoublePerturbationMain <- function(main, gpMain = NULL) {
  textGrob(main, unit(0.5,"npc"),
           unit(0.5,"npc"),
           name="title",
           gp=gpar(gpMain, fontface="bold"),
           vp=vpPath("dplayout","main"))
}

makeDoublePerturbationAxis <- function(mainEffect, dpEffect, xlab=NULL, ylab=NULL, gpAxis=NULL, axisOnOrigin = FALSE, drawBox = TRUE) {
  if (axisOnOrigin) {
    vpx = vpPath("dplayout","dpvp", "dpvpxaxis")
    vpy = vpPath("dplayout","dpvp", "dpvpyaxis")
  } else {
    vpx = vpPath("dplayout","dpvp")
    vpy = vpPath("dplayout","dpvp")
  }
  if (drawBox) {
    gl <- gList( rectGrob(name="axisbox",gp=gpAxis,vp=vpPath("dplayout","dpvp")),
                xaxisGrob(name="xaxis",gp=gpAxis,vp=vpx),
                yaxisGrob(name="yaxis",gp=gpAxis,vp=vpy)
                )
  } else {
    gl <- gList( xaxisGrob(name="xaxis",gp=gpAxis,vp=vpx),
                yaxisGrob(name="yaxis",gp=gpAxis,vp=vpy)
                )
  }
  if (!is.null(xlab)) {
    gl <- gList( gl, textGrob(xlab,
                              unit(0.5,"npc"),
                              unit(-3,"lines"),
                              name="xlab",
                              gp=gpAxis,vp=vpPath("dplayout","dpvp")) )
  }
  if (!is.null(ylab)) {
    gl <- gList( gl, textGrob(ylab,
                              unit(-3,"lines"),
                              unit(0.5,"npc"),
                              name="ylab",
                              rot=90,
                              gp=gpar(gpAxis, rot=90),vp=vpPath("dplayout","dpvp")) )
  }
  gl
}

makeDoublePerturbationVP <- function(mainEffect, dpEffect, range, main, xlab, ylab) {
  w = ifelse (is.null(ylab), 3, 4)
  h1 = ifelse (is.null(xlab), 3, 4)
  h2 = ifelse (is.null(main), 1, 2)
  vpTree(viewport(name="dplayout", layout = grid.layout(3,3,
                                      widths = unit(c(w,1,1),c("lines","null","lines")),
                                      heights = unit(c(h2,1,h1),c("lines","null","lines")),
                                      respect=TRUE)
                   ),
         vpList( viewport(name="main", layout.pos.row=1, layout.pos.col=2,
                          clip=FALSE),
                vpTree(viewport(name="dpvp", layout.pos.row=2, layout.pos.col=2,
                                xscale=range,
                                yscale=range,
                                clip=FALSE),
                       vpList( viewport(name="dpvpxaxis", y = unit(0,"native"), height=unit(range[2], "native"),
                                        xscale=range,
                                        yscale=c(0,range[2]),
                                        just="bottom",
                                        clip=FALSE),
                              viewport(name="dpvpyaxis", x = unit(0,"native"), width=unit(range[2], "native"),
                                        xscale=c(0,range[2]),
                                        yscale=range,
                                       just="left",
                                       clip=FALSE))
                        ),
                 viewport(name="dpvpclip", layout.pos.row=2, layout.pos.col=2,
                          xscale=range,
                          yscale=range,
                          clip=TRUE)
                )
          )
}

shiftboxes <- function(df) {
  n.shift = rep(0,nrow(df))
  if (nrow(df) > 1) {
    dw = dh = median(df$h) * 0.5
    df$x = df$x + 1.5*dw
    df$y = df$y + 1.5*dh
    if (dh < 1.0e-300) { dh = 1 }
    for (i in 2:nrow(df)) {
      n=1
      z=0
      while ((n > 0) & (z <= i)) {
        z=z+1
        n=0
        for (j in 1:(i-1)) {
          has.overlap.x = has.overlap.y = FALSE
          if ( ((df$x[j]+df$w[j] >= df$x[i]) & (df$x[j]+df$w[j] <= df$x[i] + df$w[i]))
              | ((df$x[j]        >= df$x[i]) & (df$x[j]         <= df$x[i] + df$w[i]))
              | ((df$x[j]        <= df$x[i]) & (df$x[j]+df$w[j] >= df$x[i] )) ) {
            has.overlap.x = TRUE
            shift.x = df$x[j]+df$w[j] - df$x[i] + dw
          }
          if ( ((df$y[j]+df$h[j] >= df$y[i]) & (df$y[j]+df$h[j] <= df$y[i] + df$h[i]))
              | ((df$y[j]        >= df$y[i]) & (df$y[j]         <= df$y[i] + df$h[i]))
              | ((df$y[j]        <= df$y[i]) & (df$y[j]+df$h[j] >= df$y[i] )) ) {
            has.overlap.y = TRUE
            shift.y = df$y[j]+df$h[j] - df$y[i] + dh
          }
          if (has.overlap.x & has.overlap.y) {
            n=n+1
            n.shift[i] = n.shift[i] + 1
            df$is.shifted[i] = TRUE
            df$x[i] = df$x[i] + min(shift.x,shift.y)
            df$y[i] = df$y[i] + min(shift.x,shift.y)
          }
        }
      }
    }
  }
  return(df)
}

postDrawDetails.doublePerturbation <- function(x, recording) {
  if (x$avoid.overlap) {
    n = downViewport("dpvp")
    if (nrow(x$labels$top) > 0) {
      w=h=rep(0.0,nrow(x$labels$top))
      for (i in 1:nrow(x$labels$top)) {
        w[i] = convertWidth(stringWidth(x$labels$top$label[i]),"native",valueOnly = TRUE)
        h[i] = convertHeight(stringHeight(x$labels$top$label[i]),"native",valueOnly = TRUE)
      }
      df = data.frame(name=x$labels$top$name,label=x$labels$top$label,x=x$labels$top$x,y=x$labels$top$y,w=w,h=h,is.shifted=FALSE,stringsAsFactors=FALSE)
      df$x = -1 * df$x
      df = shiftboxes(df)
      df$x = -1 * df$x
      grid.polyline(unit(c(df$x,x$labels$top$x),"native"),unit(c(df$y,x$labels$top$y),"native"),
                    id=rep(1:nrow(df),2),gp=gpar(col="lightgray"),name="labellinetop")
      grid.text(df$label,unit(df$x,"native"),unit(df$y,"native"),
                name="labeltop", just=c("right","bottom"))
    }
    if (nrow(x$labels$bottom) > 0) {
      w=h=rep(0.0,nrow(x$labels$bottom))
      for (i in 1:nrow(x$labels$bottom)) {
        w[i] = convertWidth(stringWidth(x$labels$bottom$label[i]),"native",valueOnly = TRUE)
        h[i] = convertHeight(stringHeight(x$labels$bottom$label[i]),"native",valueOnly = TRUE)
      }
      df = data.frame(name=x$labels$bottom$name,label=x$labels$bottom$label,x=x$labels$bottom$x,y=x$labels$bottom$y,w=w,h=h,is.shifted=FALSE,stringsAsFactors=FALSE)
      df$y = -1 * df$y
      df = shiftboxes(df)
      df$y = -1 * df$y
      grid.polyline(unit(c(df$x,x$labels$bottom$x),"native"),unit(c(df$y,x$labels$bottom$y),"native"),
                    id=rep(1:nrow(df),2),gp=gpar(col="lightgray"),name="labellinebottom")
      grid.text(df$label,unit(df$x,"native"),unit(df$y,"native"),
                name="labelbottom", just=c("left","top"))
    }
    upViewport(n)
  }
}

doublePerturbationGrob <- function( mainEffect, dpEffect, mainEffectTarget, range=NULL,
                                   main=NULL, xlab=NULL, ylab=NULL,
                                   text=NULL, avoid.overlap=TRUE,
                                   axisOnOrigin = FALSE,
                                   drawBox = TRUE,
                                   pch = 21, size=unit(1, "char"), fill = NULL,
                                   gpMain = gpar(lty="dashed", lwd=3, col="cyan"),
                                   gpNI = gpar(lty="dashed", lwd=3, col="orange"),
                                   gpPoints = gpar(pch=21),
                                   gpText = NULL,
                                   gpAxis = NULL,
                                   gpWTLines=NULL,
                                   name=NULL, gp=NULL, vp=NULL) {
  if (is.null(range)) {
    rangex = range(mainEffect[is.finite(mainEffect)])
    rangey = range(dpEffect[is.finite(mainEffect)])
    range = c(min(rangex[1],rangey[1]),max(rangex[2],rangey[2]))
    dr = diff(range)
    range = range + 0.05*c(-dr,dr)
  }
  labels = getDoublePerturbationLabels(mainEffect, dpEffect, mainEffectTarget, text)
  igt <- gTree(mainEffect=mainEffect, dpEffect=dpEffect, mainEffectTarget,labels=labels, avoid.overlap=avoid.overlap,
               range=range, pch=pch, size=size, fill=fill,
               gpAxis=gpAxis,
               gpWTLines=gpWTLines,
               gpMain=gpMain,
               gpNI=gpNI,
               gpPoints=gpPoints,
               gpText=gpText,
               name=name, gp=gp, vp=vp,
               childrenvp=makeDoublePerturbationVP(mainEffect, dpEffect, range, main, xlab, ylab),
               children=gList(
                 makeDoublePerturbationAxis(mainEffect, dpEffect, xlab=xlab, ylab=ylab, gpAxis=gpAxis, axisOnOrigin = axisOnOrigin, drawBox = drawBox),
                 makeDoublePerturbationWTlines(mainEffect, dpEffect, gpWTLines=gpWTLines),
                 makeDoublePerturbationMainEffects(mainEffectTarget, range, gpMain=gpMain),
                 makeDoublePerturbationNI(mainEffectTarget, range, gpNI=gpNI),
                 makeDoublePerturbationPoints(mainEffect, dpEffect, pch=pch, size=size, fill=fill, gpPoints=gpPoints),
                 makeDoublePerturbationMain(main)
##                  ifelse(avoid.overlap, gList(list()), makeDoublePerturbationLabels(labels, gpText = gpText) )
                 ),
               cl="doublePerturbation")
  igt
}

grid.doublePerturbation <- function(..., draw = TRUE) {
  igt <- doublePerturbationGrob(...)
  if (draw) {
    grid.draw(igt)
  }
  igt
}

