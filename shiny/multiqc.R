shiny.wgbs.serveMultiqcReports <- function (multiqcPath, datasets) {

  dir.create(multiqcPath, showWarnings = FALSE)

  for (datasetName in names(datasets)) {

    multiqcFile <- datasets[[datasetName]]$fullReport

    targetFile <- paste(multiqcPath, "/", datasetName, ".html", sep = "")

    file.copy(multiqcFile, targetFile)

  }

  addResourcePath("multiqc", multiqcPath)

}

shiny.wgbs.serveQualimapReports <- function (datasets) {

  for (datasetName in names(datasets)) {

    addResourcePath(paste0("qualimap-", datasetName), datasets[[datasetName]]$qualimapDirectory)

  }

}