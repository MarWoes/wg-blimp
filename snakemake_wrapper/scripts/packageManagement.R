if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")
}


options(repos = c(CRAN = "https://cran.rstudio.com"))

ensurePackageInstallation <- function (packageName, type = "CRAN") {

  isPackageInstalled <- packageName %in% rownames(installed.packages())

  if (!isPackageInstalled) {

    if (type == "CRAN") {

      install.packages(packageName, quiet = TRUE)

    } else if (type == "bioc") {

      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(packageName)

    }

  }

}

loadOrInstall <- function (packageName, type = "CRAN") {

  ensurePackageInstallation(packageName, type)

  library(packageName, character.only = TRUE)
}

# install dependencies

cranPackages <- c(
  "data.table",
  "KernSmooth",
  "parallel",
  "ggplot2",
  "stringr"
)

biocPackages <- c(
  "BSgenome",
  "bsseq",
  "GenomicRanges",
  "MethylSeekR",
  "rtracklayer"
)

for (package in cranPackages) {
  loadOrInstall(package)
}

for (package in biocPackages) {
  loadOrInstall(package, type = "bioc")
}

# prepare shiny app
ensurePackageInstallation("shiny")

writeLines(capture.output(sessionInfo()), snakemake@output$r_session_info)
