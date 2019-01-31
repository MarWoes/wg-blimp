if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")
}

library(data.table)

wgbs.metilene.transformMethylationValues <- function (inputFile, outputFile) {

  fileContent <- fread(inputFile, data.table = FALSE, stringsAsFactors = FALSE)

  # see https://github.com/dpryan79/MethylDackel#single-cytosine-methylation-metrics-extraction
  methylatedReadsColumnIndex   <- 5
  unmethylatedReadsColumnIndex <- 6

  methylationRatios <- fileContent[,methylatedReadsColumnIndex] / (fileContent[,methylatedReadsColumnIndex] + fileContent[,unmethylatedReadsColumnIndex])

  transformedValues <- data.frame(
    chr = fileContent[,1],
    start = fileContent[,2],
    end = fileContent[,3],
    value = methylationRatios,
    stringsAsFactors = FALSE
  )

  write.table(transformedValues, file = outputFile, col.names = FALSE, quote = FALSE, sep = "\t", row.names = FALSE)
}

if (exists("snakemake")) {
  wgbs.metilene.transformMethylationValues(
    snakemake@input$bedGraph,
    snakemake@output$bedGraph
  )

  # save.image("scripts/bedgraph-snakemake.RData")
}
