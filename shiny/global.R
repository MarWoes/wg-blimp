options(repos = c(CRAN = "https://cran.rstudio.com"))

library(shiny)
library(shinydashboard)
library(data.table)
library(ggplot2)
library(htmlwidgets)
library(DT)
library(httpuv)

source("datasets.R")
source("multiqc.R")
source("bamServing.R")
