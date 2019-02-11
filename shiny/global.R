options(repos = c(CRAN = "https://cran.rstudio.com"))

shiny.wgbs.loadOrInstall <- function (packageName, type = "CRAN") {

  isPackageInstalled <- packageName %in% rownames(installed.packages())

  if (!isPackageInstalled) {

    if (type == "CRAN") {

      install.packages(packageName)

    } else if (type == "bioc") {

      source("https://bioconductor.org/biocLite.R")
      biocLite(packageName)

    }

  }

  library(packageName, character.only = TRUE)
}

cranPackages <- c(
  "shiny",
  "shinydashboard",
  "data.table",
  "httr",
  "ggplot2",
  "htmlwidgets",
  "DT",
  "httpuv"
)

biocPackages <- c(
  "biomaRt",
  "rentrez"
)

for (package in cranPackages) {

  shiny.wgbs.loadOrInstall(package)

}

for (package in biocPackages) {

  shiny.wgbs.loadOrInstall(package, type = "bioc")

}

shiny.wgbs.annotatedDMRColumns <- c(
  "chr",
  "start",
  "end",
  "num_cpgs",
  "diff",
  "qValue",
  "tool",
  "length",
  "gene_overlap",
  "gene_name",
  "cgi_overlap",
  "promoter_overlap",
  "promoter_name",
  "num_repeats",
  "mean_cov"
)

shiny.wgbs.annotatedDMRCharacterColumns <- c(
  "chr",
  "tool",
  "gene_name",
  "promoter_name"
)

source("datasets.R")
source("multiqc.R")
source("bamServing.R")