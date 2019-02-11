options(repos = c(CRAN = "https://cran.rstudio.com"))

shiny.wgbs.loadOrInstall <- function (packageName, type = "CRAN") {

  isPackageInstalled <- packageName %in% rownames(installed.packages())

  if (!isPackageInstalled) {

    if (type == "CRAN") {

      install.packages(packageName)

    } else if (type == "bioc") {

      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(packageName, version = "3.8")

    }

  }

  library(packageName, character.only = TRUE)
}

cranPackages <- c(
  "shiny",
  "shinydashboard",
  "data.table",
  "ggplot2",
  "htmlwidgets",
  "DT",
  "httpuv"
)

for (package in cranPackages) {

  shiny.wgbs.loadOrInstall(package)

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