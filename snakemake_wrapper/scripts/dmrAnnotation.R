if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")

  snakemake@source("regionAnnotation.R")
}


wgbs.annotateDMRs <- function (dmrFile, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, gzippedCoverageFiles, allowedBiotypes, promoterTSSDistances, annotatedDmrFile) {

  dmrs <- fread(dmrFile)

  dmrs <- annotation.annotateRegions(dmrs, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, allowedBiotypes, promoterTSSDistances)

  dmrLocationKeys <- paste0(dmrs$chr, ":", dmrs$start, "-", dmrs$end)

  coverages <- sapply(gzippedCoverageFiles, function (gzippedCoverageFile) {

    coverage = fread(cmd = paste("zcat", gzippedCoverageFile))
    coverageLocationKeys = paste0(coverage$V1, ":", coverage$V2, "-", coverage$V3)

    return(coverage$V4[match(dmrLocationKeys, coverageLocationKeys)])
  })

  dmrs$mean_cov <- apply(coverages, 1, mean)

  write.table(dmrs, file = annotatedDmrFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")

  return(dmrs)
}

if (exists("snakemake")) {

   # save.image(file = "dmr-annotation-snakemake.Rdata")
   wgbs.annotateDMRs(
     snakemake@input$combined_dmrs,
     snakemake@params$cgi_annotation_file,
     snakemake@params$gene_annotation_file,
     snakemake@params$repeat_masker_annotation_file,
     snakemake@input$coverages,
     snakemake@params$biotypes,
     snakemake@params$tss_distances,
     snakemake@output$annotated_dmrs
   )

}
