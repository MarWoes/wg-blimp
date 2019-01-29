library(bsseq)

callDmrs <- function (methylDackelBedGraphFiles, sampleNames, group1Samples, group2Samples, threads, min_cpg, min_diff, rdatFile, csvFile, pdfFile) {

  methylationData <- read.bismark(methylDackelBedGraphFiles,
                                  sampleNames,
                                  strandCollapse = FALSE,
                                  fileType = "cov",
                                  mc.cores = threads)

  smoothedData <- BSmooth(methylationData, mc.cores = threads, verbose = TRUE)

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
  save.image(file = "snakemake.Rdata")
}
