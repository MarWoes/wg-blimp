library(data.table)
library(GenomicRanges)
library(stringr)

annotation.isExistingFileOrNone <- function (fileName) {

  if (fileName == "None") {
    return(FALSE)
  }

  if (file.exists(fileName)) {
    return(TRUE)
  }

  stop(paste0("Annotation file not found: ", fileName))

}

annotation.annotateOverlap <- function (regionRanges, overlapRanges, geneTable) {

  intersection <- as.list(findOverlaps(regionRanges, overlapRanges))
  isOverlapping <- sapply(intersection, length) > 0
  overlappingGeneNames <- sapply(intersection, function (indices) paste(unique(geneTable$external_gene_name[indices]), collapse = ","))

  return(data.table(
    isOverlapping,
    overlappingGeneNames
  ))
}

annotation.annotateGenes <- function (regions, regionRanges, genes) {

  geneRanges <- GRanges(seqnames = genes$chromosome_name, ranges = IRanges(genes$start_position, genes$end_position))

  regions[,c("gene_overlap", "gene_name")] <- annotation.annotateOverlap(regionRanges, geneRanges, genes)

  return(regions)
}

annotation.annotateCGIslands <- function (regions, regionRanges, cgis) {

  cgiRanges  <- GRanges(seqnames = cgis$chrom, ranges = IRanges(cgis$chromStart, cgis$chromEnd))

  regions$cgi_overlap <- countOverlaps(regionRanges, cgiRanges) > 0

  return(regions)
}

annotation.annotatePromoters <- function (regions, regionRanges, transcriptStartSites, promoterTSSDistances) {

  promoterStarts <- transcriptStartSites$transcription_start_site + promoterTSSDistances[1]
  promoterEnds <- transcriptStartSites$transcription_start_site + promoterTSSDistances[2]

  promoterRanges <- GRanges(seqnames = transcriptStartSites$chromosome_name, ranges = IRanges(promoterStarts, promoterEnds))

  regions[,c("promoter_overlap", "promoter_name")] <- annotation.annotateOverlap(regionRanges, promoterRanges, transcriptStartSites)

  return(regions)
}

annotation.annotateRepeats <- function (regions, regionRanges, repeats, classes) {

  relevantRepeats <- repeats[repClass %in% classes]

  repeatRanges <- GRanges(relevantRepeats$genoName, IRanges(relevantRepeats$genoStart, relevantRepeats$genoEnd))

  regions$num_repeats <- countOverlaps(regionRanges, repeatRanges)

  return(regions)
}

annotation.annotateRegions <- function (regionTable, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, gzippedTranscriptionStartSiteFile, allowedBiotypes, promoterTSSDistances) {

  regionTable$length <- regionTable$end - regionTable$start

  # remove any trailing 'chr's for range intersection
  regionRanges  <- GRanges(seqnames = str_remove(regionTable$chr, "^chr"), ranges = IRanges(regionTable$start, regionTable$end))

  if (annotation.isExistingFileOrNone(gzippedCgiFile)) {

    cgis <- fread(cmd = paste("zcat", gzippedCgiFile))
    cgis$chrom <- substr(cgis$chrom, 4, nchar(cgis$chrom))
    regionTable <- annotation.annotateCGIslands(regionTable, regionRanges, cgis)

  }

  if (annotation.isExistingFileOrNone(gzippedGeneFile)) {

    genes <- fread(cmd = paste("zcat", gzippedGeneFile))
    genes <- genes[gene_biotype %in% allowedBiotypes]
    regionTable <- annotation.annotateGenes(regionTable, regionRanges, genes)

  }

  if (annotation.isExistingFileOrNone(gzippedTranscriptionStartSiteFile)) {

    transcriptStartSites <- fread(cmd = paste("zcat", gzippedTranscriptionStartSiteFile))
    transcriptStartSites <- transcriptStartSites[gene_biotype %in% allowedBiotypes]
    regionTable <- annotation.annotatePromoters(regionTable, regionRanges, transcriptStartSites, promoterTSSDistances)

  }

  if (annotation.isExistingFileOrNone(gzippedRepeatMaskerAnnotationFile)) {

    repeats <- fread(cmd = paste("zcat", gzippedRepeatMaskerAnnotationFile))
    repeats$genoName <- substr(repeats$genoName, 4, nchar(repeats$genoName))
    regionTable <- annotation.annotateRepeats(regionTable, regionRanges, repeats, c("LINE", "SINE", "LTR", "DNA"))

  }

  return(regionTable)
}
