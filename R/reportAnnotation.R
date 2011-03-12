
reportAnnotation <- function(sgi, verbose = 0, path = ".", dir = "annotation", prefix = "annotation", report = NULL) {
  stopifnot( is( sgi, "RNAinteract" ) )

  if (verbose > 0) {
    time1 <- Sys.time()
    message("Write annotation. \r", appendLF = FALSE)
  }
  dir.create(sprintf("%s/%s", path, dir), showWarnings = FALSE)

  if (!is.null(report)) {
    hwrite(sprintf("<b>%s</b>",prefix), br=TRUE, page=report)
    hwrite("[targets]", link = sprintf("%s/index-%s-targets.html",dir, prefix), target = "showframe", br=TRUE, page=report)
    hwrite("[reagents]", link = sprintf("%s/index-%s-reagents.html",dir, prefix), target = "showframe", br=FALSE, page=report)
    hwrite("<br>", br=TRUE, page=report)
  }

  p = openPage(sprintf("%s/%s/index-%s-targets.html",path, dir, prefix))
  
  write.table(sgi@targets, file = sprintf("%s/%s/%s-targets.txt",path, dir, prefix), row.names = FALSE, sep="\t", quote=FALSE)

  hwrite("<h3>Targets</h3>",table=FALSE,br=FALSE,page=p)
  hwrite(sprintf("Download as text-file: <a href=%s-targets.txt>[targets]</a> ",prefix),table=FALSE,br=FALSE,page=p)
  hwrite(sgi@targets,style='text-align:right',row.bgcolor=rep(c("#CCCCCC","#EEEEEE"),length.out=nrow(sgi@targets)), cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

  closePage(p,splash=FALSE)

  p = openPage(sprintf("%s/%s/index-%s-reagents.html",path, dir, prefix))

  write.table(sgi@reagents, file = sprintf("%s/%s/%s-reagents.txt",path, dir, prefix), row.names = FALSE, sep="\t", quote=FALSE)

  hwrite("<h3>Reagents</h3>",table=FALSE,br=FALSE,page=p)
  hwrite(sprintf("Download as text-file: <a href=%s-reagents.txt>[reagents]</a> ",prefix),table=FALSE,br=FALSE,page=p)
  hwrite(sgi@reagents,style='text-align:right',row.bgcolor=rep(c("#CCCCCC","#EEEEEE"),length.out=nrow(sgi@reagents)), cellpadding=3, cellspacing=0,border=1,br=TRUE,page=p)

  closePage(p,splash=FALSE)


  if (verbose > 0) {
    time2 = Sys.time()
    message(sprintf("Write annotation. Done. Time: %s.",format(difftime(time2, time1, tz = "",units = c("auto")))))
  }
  invisible(NULL)

}

