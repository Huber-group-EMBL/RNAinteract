
startReport <- function(outputpath) {
  dir.create(outputpath, showWarnings = FALSE)
  con <- file(sprintf("%s/index.html", outputpath), "wb")
  writeLines("<html>", con)
  writeLines("<frameset cols=\"25%,75%\">", con)
  writeLines("  <frame src=\"frame_left.html\">", con)
  writeLines("  <frame src=\"frame_right.html\" name=showframe>", con)
  writeLines("</frameset>", con)
  writeLines("</html>", con)
  close(con)

  report = openPage(sprintf("%s/frame_left.html", outputpath))
  return(report)
}

endReport <- function(report) {
  closePage(report, splash=FALSE)
}


