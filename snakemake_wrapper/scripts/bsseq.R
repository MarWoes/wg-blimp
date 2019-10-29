if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")
}

library(bsseq)

callDmrs <- function (methylDackelBedGraphFiles, sampleNames, group1Samples, group2Samples, threads, min_cpg, min_diff, rdatFile, csvFile, pdfFile) {

  biocParallelParam <- MulticoreParam(workers = threads)
  
  methylationData <- read.bismark(methylDackelBedGraphFiles,
                                  strandCollapse = FALSE,
                                  BPPARAM = biocParallelParam)
  
  sampleNames(methylationData) <- sampleNames

  # remove all NaN CpG loci: this should not happen in a WGBS experiment anyway, so the discarded
  # regions should be considered low-quality
  methylationData <- methylationData[!is.nan(rowSums(getMeth(methylationData, type = "raw")))]

  smoothedData <- BSmooth(methylationData, BPPARAM = biocParallelParam, verbose = TRUE)

  # add some color for later plotting
  pData <- pData(smoothedData)
  pData$col <- c(rep("red", length(group1Samples)), rep("blue", length(group2Samples)))
  pData(smoothedData) <- pData

  # filter out rows containing NA's
  invalidRows <- is.na(rowSums(getMeth(smoothedData)))

  print(paste("Filtering out", sum(invalidRows), "rows containing NA"))

  tstats <- BSmooth.tstat(smoothedData[!invalidRows],
                          group1 = group1Samples,
                          group2 = group2Samples,
                          estimate.var = "same",
                          mc.cores = threads)

  dmrs0 <- dmrFinder(tstats)

  write.table(dmrs0, file = csvFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

  pdf(file = pdfFile, width = 10, height = 5)

  dmrs <- subset(dmrs0, n >= min_cpg & abs(meanDiff) >= min_diff)

  maxNumOfRegions <- ifelse(nrow(dmrs) > 100, 100, nrow(dmrs))
  plotManyRegions(smoothedData[!invalidRows], regions = dmrs[seq_len(maxNumOfRegions),], addRegions = dmrs, extend = 5000)

  dev.off()
  save.image(file = rdatFile)
}

if (exists("snakemake")) {

  #save.image(file = "snakemake.Rdata")
  callDmrs(snakemake@input$meth,
           snakemake@config$samples,
           snakemake@config$group1,
           snakemake@config$group2,
           snakemake@threads,
           snakemake@config$min_cpg,
           snakemake@config$min_diff,
           snakemake@output$rdata,
           snakemake@output$csv,
           snakemake@output$pdf)
  

}
