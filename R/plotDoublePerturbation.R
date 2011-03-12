
plotDoublePerturbation <- function( sgi, screen, channel, target, 
                                  withoutgroups=c("neg", "pos"), design = "template",
                                  main = sprintf("Double Perturbation Plot of %s, screen %s (%s)", target, screen, channel),
                                  xlab = "single perturbation level",
                                  ylab = "double perturbation level",
                                  range = NULL,
                                  show.labels = "none", label.par = 0.05, label = NULL, avoid.overlap = TRUE,
                                  col = NULL, fill = list(colors=c("cornflowerblue","black","yellow")),
                                  D=NULL, MT=NULL, MQ=NULL, PV=NULL, QV=NULL, PI=NULL,
                                  ... ) {
  stopifnot( is( sgi, "RNAinteract" ))

  if (is.null(D)) {
    D <- getData(sgi, screen = screen, channel = channel,
                 type="data", format="targetMatrix", normalized = TRUE,
                 mixTemplateQuery = FALSE,
                 withoutgroups=withoutgroups)
    MainNeg <- getMainNeg(sgi, screen = screen, channel = channel)
    D <- D - MainNeg
  }
  if (is.null(MT)) {
    MT <- getData(sgi, screen = screen, channel = channel,
                  type="main", format="targetMatrix", normalized = TRUE,
                  mixTemplateQuery = FALSE,
                  design="template", withoutgroups=withoutgroups)
  }
  if (is.null(MQ)) {
    MQ <- getData(sgi, screen = screen, channel = channel,
                  type="main", format="targetMatrix", normalized = TRUE,
                  mixTemplateQuery = FALSE,
                  design="query", withoutgroups=withoutgroups)
  }
  useit <- which(row.names(D) != target)
  d <- D[useit,target]
  mt <- MT[useit,target]
  mq <- MQ[1,target]

  text = NULL
  if (!is.null(label)) {
    if (is.null(names(label))) {
      text = names(d)
      text[!(text %in% label)] = ""
    } else {
      text = rep("",length(d))
      m = match(names(label), names(d))
      text[m[is.finite(m)]] = label[is.finite(m)]
    }
  } else {
    if (!is.null(show.labels)) {
      text = names(d)
      text = switch( show.labels,
        none = rep("",length(d)),
        all = text,
        p.value = {
          if (is.null(PV)) {
            PV <- getData(sgi, screen = screen, channel = channel,
                          type="p.value", format="targetMatrix",
                          mixTemplateQuery = FALSE,
                          withoutgroups=withoutgroups)
          }
          text[PV[useit,target] > label.par] = ""
          text
        },
        q.value = {
          if (is.null(QV)) {
            QV <- getData(sgi, screen = screen, channel = channel,
                          type="q.value", format="targetMatrix",
                          mixTemplateQuery = FALSE,
                          withoutgroups=withoutgroups)
          }
          text[QV[useit,target] > label.par] = ""
          text
        }
        )
    }
  }
  if (!is.null(col)) {
    if (is.null(names(col))) {
      fill = NULL
    } else {
      fill = list()
      fill$col = rep("lightgray",length(d))
      m = match(names(col), names(d))
      fill$col[m[is.finite(m)]] = col[is.finite(m)]
    }    
  } else {
    if (is.list(fill)) {
      if (!any(c("colors","colramp") %in% names(fill))) {
        stop("list fill has to have either a colors or a colramp entry")
      }
      if (!all(c("values","at") %in% names(fill))) {
        if (is.null(PI)) {
          PI <- getData(sgi, screen = screen, channel = channel,
                        type="pi", format="targetMatrix",
                        mixTemplateQuery = FALSE,
                        withoutgroups=withoutgroups)
        }
        r <- mad(PI, center=0.0, na.rm=TRUE)
        pi <- PI[useit,target]
        if (r > 1.0e-200) { pi <- pi / r}
        r = range(abs(pi),finite=TRUE)
        if (!("values" %in% names(fill))) {
          fill$values = pi
        }
        if (!("at" %in% names(fill))) {
          at = seq(1,3,length.out=127)
          fill$at = c(r[1]-4,-rev(at),at,r[2]+4)
        }
      }
    } else {
      fill = NULL
    }
  }

  if (is.null(range)) {
    range <- range(range(D,finite=TRUE),range(MT,finite=TRUE),range(MQ,finite=TRUE))
    dr <- diff(range)
    range <- range + 0.1*c(-dr,dr)
  }
  scale <- getScale(sgi,channel)

  mainEffectTarget <- mq
  mainEffect <- mt
  dpEffect <- d

  igt <- grid.doublePerturbation(mainEffect=mainEffect,
                                 dpEffect=dpEffect,
                                 mainEffectTarget=mainEffectTarget,
                                 main=main, xlab=xlab, ylab=ylab, text=text,
                                 fill=fill,
                                 range=range, ...)
  igt
}

