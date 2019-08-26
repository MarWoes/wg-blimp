if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")
}

library(data.table)
library(GenomicRanges)
library(stringr)

wgbs.annotateOverlap <- function (dmrRanges, overlapRanges, geneTable) {

  intersection <- as.list(findOverlaps(dmrRanges, overlapRanges))
  isOverlapping <- sapply(intersection, length) > 0
  overlappingGeneNames <- sapply(intersection, function (indices) paste(unique(geneTable$external_gene_name[indices]), collapse = ","))

  return(data.table(
    isOverlapping,
    overlappingGeneNames
  ))
}

wgbs.annotateGenes <- function (dmrs, dmrRanges, genes) {

  geneRanges <- GRanges(seqnames = genes$chromosome_name, ranges = IRanges(genes$start_position, genes$end_position))

  dmrs[,c("gene_overlap", "gene_name")] <- wgbs.annotateOverlap(dmrRanges, geneRanges, genes)

  return(dmrs)
}

wgbs.annotateCGIslands <- function (dmrs, dmrRanges, cgis) {

  cgiRanges  <- GRanges(seqnames = cgis$chrom, ranges = IRanges(cgis$chromStart, cgis$chromEnd))

  dmrs$cgi_overlap <- countOverlaps(dmrRanges, cgiRanges) > 0

  return(dmrs)
}

wgbs.annotatePromoters <- function (dmrs, dmrRanges, transcriptStartSites, promoterTSSDistances) {

  promoterStarts <- transcriptStartSites$transcription_start_site + promoterTSSDistances[1]
  promoterEnds <- transcriptStartSites$transcription_start_site + promoterTSSDistances[2]

  promoterRanges <- GRanges(seqnames = transcriptStartSites$chromosome_name, ranges = IRanges(promoterStarts, promoterEnds))

  dmrs[,c("promoter_overlap", "promoter_name")] <- wgbs.annotateOverlap(dmrRanges, promoterRanges, transcriptStartSites)

  return(dmrs)
}

wgbs.annotateRepeats <- function (dmrs, dmrRanges, repeats, classes) {

  relevantRepeats <- repeats[repClass %in% classes]

  repeatRanges <- GRanges(relevantRepeats$genoName, IRanges(relevantRepeats$genoStart, relevantRepeats$genoEnd))

  dmrs$num_repeats <- countOverlaps(dmrRanges, repeatRanges)

  return(dmrs)
}

wgbs.annotateDMRs <- function (dmrFile, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, gzippedTranscriptionStartSiteFile, gzippedCoverageFiles, allowedBiotypes, promoterTSSDistances, annotatedDmrFile) {

  dmrs <- fread(dmrFile)

  cgis <- fread(paste("zcat", gzippedCgiFile))
  genes <- fread(paste("zcat", gzippedGeneFile))
  transcriptStartSites <- fread(paste("zcat", gzippedTranscriptionStartSiteFile))
  repeats <- fread(paste("zcat", gzippedRepeatMaskerAnnotationFile))

  cgis$chrom <- substr(cgis$chrom, 4, nchar(cgis$chrom))
  repeats$genoName <- substr(repeats$genoName, 4, nchar(repeats$genoName))

  genes <- genes[gene_biotype %in% allowedBiotypes]
  transcriptStartSites <- transcriptStartSites[gene_biotype %in% allowedBiotypes]

  dmrs$length <- dmrs$end - dmrs$start

  # remove any trailing 'chr's for range intersection
  dmrRanges  <- GRanges(seqnames = str_remove(dmrs$chr, "^chr"), ranges = IRanges(dmrs$start, dmrs$end))

  dmrs <- wgbs.annotateGenes(dmrs, dmrRanges, genes)
  dmrs <- wgbs.annotateCGIslands(dmrs, dmrRanges, cgis)
  dmrs <- wgbs.annotatePromoters(dmrs, dmrRanges, transcriptStartSites, promoterTSSDistances)
  dmrs <- wgbs.annotateRepeats(dmrs, dmrRanges, repeats, c("LINE", "SINE", "LTR", "DNA"))

  coverages <- sapply(gzippedCoverageFiles, function (gzippedCoverageFile) fread(paste("zcat", gzippedCoverageFile))$V4)

  dmrs$mean_cov <- apply(coverages, 1, mean)

  write.table(dmrs, file = annotatedDmrFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")

  return(dmrs)
}

if (exists("snakemake")) {

   # save.image(file = "scripts/dmr-annotation-snakemake.Rdata")
   wgbs.annotateDMRs(
     snakemake@input$combined_dmrs,
     snakemake@input$cgi_annotation_file,
     snakemake@input$gene_annotation_file,
     snakemake@input$repeat_masker_annotation_file,
     snakemake@input$transcript_start_site_file,
     snakemake@input$coverages,
     snakemake@params$biotypes,
     snakemake@params$tss_distances,
     snakemake@output$annotated_dmrs
   )

}
