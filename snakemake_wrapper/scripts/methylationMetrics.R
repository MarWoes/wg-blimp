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

wgbs.methylation.computeMethylationMetrics <- function(bedGraphFiles, methylationMetricsFile, conversionChromosomes) {

  methylationValues <- lapply(bedGraphFiles, fread, skip = 1)

  numTrailingCharacters <- nchar("_CpG.bedGraph")

  sampleNames <- basename(bedGraphFiles)
  sampleNames <- substr(sampleNames, 1, nchar(sampleNames) - numTrailingCharacters)

  names(methylationValues) <- sampleNames

  methylationMetrics <- data.table(
    sample = names(methylationValues),
    `methylation (total)` = sapply(methylationValues, wgbs.methylation.computeMeanMethylation)
  )

  if (length(conversionChromosomes) > 0) {
  
    conversionRates <- data.table(do.call(rbind, sapply(methylationValues, function (sampleValues) {
  
      chromosomeConversionRates <- sapply(conversionChromosomes, function(chr) {

        valuesOnChromosome <- sampleValues[V1 == chr]
  
        return(wgbs.methylation.computeMeanMethylation(valuesOnChromosome))
  
      })
  
      return(chromosomeConversionRates)
    }, simplify = FALSE)))

    colnames(conversionRates) <- paste("methylation (", conversionChromosomes, ")", sep = "")
    
    methylationMetrics <- cbind(methylationMetrics, conversionRates)
  }

  write.table(methylationMetrics, methylationMetricsFile, row.names = FALSE)

  return(methylationMetrics)
}

if (exists("snakemake")) {
  
  # save.image("conversion-snakemake.Rdata")
  
  wgbs.methylation.computeMethylationMetrics(
    snakemake@input$bed_graphs,
    snakemake@output$methylation_metrics,
    snakemake@params$methylation_rate_on_chromosomes
  )
}
