if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")
}

library(data.table)

wgbs.methylation.computeMeanMethylation <- function(methylationTable) {

  methylatedValues <- methylationTable$V5
  allValues <- methylationTable$V5 + methylationTable$V6

  return(mean(methylatedValues / allValues))

}

wgbs.methylation.computeMethylationMetrics <- function(bedGraphFiles, methylationMetricsFile, conversionChromosomes = c("MT", "LAMBDA")) {

  methylationValues <- lapply(bedGraphFiles, fread, skip = 1)

  numTrailingCharacters <- nchar("_CpG.bedGraph")

  sampleNames <- basename(bedGraphFiles)
  sampleNames <- substr(sampleNames, 1, nchar(sampleNames) - numTrailingCharacters)

  names(methylationValues) <- sampleNames

  conversionRates <- sapply(methylationValues, function (sampleValues) {

    chromosomeConversionRates <- sapply(conversionChromosomes, function(chr) {

      valuesOnChromosome <- sampleValues[V1 == chr]

      return(1 - wgbs.methylation.computeMeanMethylation(valuesOnChromosome))

    })

    return(chromosomeConversionRates)
  })

  rownames(conversionRates) <- paste("conversion (", conversionChromosomes, ")", sep = "")

  methylationMetrics <- cbind(data.table(sample = names(methylationValues)), t(conversionRates))
  methylationMetrics$methylation <- sapply(methylationValues, wgbs.methylation.computeMeanMethylation)

  write.table(methylationMetrics, methylationMetricsFile, row.names = FALSE)

  return(methylationMetrics)
}

if (exists("snakemake")) {
  wgbs.methylation.computeMethylationMetrics(
    snakemake@input$bed_graphs,
    snakemake@output$methylation_metrics
  )

  # save.image("scripts/conversion-snakemake.Rdata")
}
