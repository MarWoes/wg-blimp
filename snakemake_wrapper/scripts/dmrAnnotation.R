if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")

  snakemake@source("regionAnnotation.R")
}


wgbs.annotateDMRs <- function (dmrFile, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, gzippedTranscriptionStartSiteFile, gzippedCoverageFiles, allowedBiotypes, promoterTSSDistances, annotatedDmrFile) {

  dmrs <- fread(dmrFile)

  dmrs <- annotation.annotateRegions(dmrs, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, gzippedTranscriptionStartSiteFile, allowedBiotypes, promoterTSSDistances)

  coverages <- sapply(gzippedCoverageFiles, function (gzippedCoverageFile) fread(paste("zcat", gzippedCoverageFile))$V4)

  dmrs$mean_cov <- apply(coverages, 1, mean)

  write.table(dmrs, file = annotatedDmrFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")

  return(dmrs)
}

if (exists("snakemake")) {

   # save.image(file = "scripts/dmr-annotation-snakemake.Rdata")
   wgbs.annotateDMRs(
     snakemake@input$combined_dmrs,
     snakemake@params$cgi_annotation_file,
     snakemake@params$gene_annotation_file,
     snakemake@params$repeat_masker_annotation_file,
     snakemake@params$transcript_start_site_file,
     snakemake@input$coverages,
     snakemake@params$biotypes,
     snakemake@params$tss_distances,
     snakemake@output$annotated_dmrs
   )

}
