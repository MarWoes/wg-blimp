library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

annotation.isExistingFileOrNone <- function (fileName) {

  if (is.null(fileName)) {
    return(FALSE)
  }

  if (length(fileName) == 1 && file.exists(fileName)) {
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

annotation.annotateGenes <- function (regions, regionRanges, gtfRanges) {

  geneGtfRanges = gtfRanges[gtfRanges$type == "gene"]

  geneNames <- data.table(external_gene_name = geneGtfRanges$gene_name)

  regions[,c("gene_overlap", "gene_name")] <- annotation.annotateOverlap(regionRanges, geneGtfRanges, geneNames)

  return(regions)
}

annotation.annotateCGIslands <- function (regions, regionRanges, cgis) {

  cgiRanges  <- GRanges(seqnames = cgis$chrom, ranges = IRanges(cgis$chromStart, cgis$chromEnd))

  regions$cgi_overlap <- countOverlaps(regionRanges, cgiRanges) > 0

  return(regions)
}

annotation.annotatePromoters <- function (regions, regionRanges, gtfRanges, promoterTSSDistances) {

  transcriptRanges = gtfRanges[gtfRanges$type == "transcript"]
  geneNames = data.table(external_gene_name = transcriptRanges$gene_name)

  promoterStarts <- ifelse(
    strand(transcriptRanges) == "+",
    start(transcriptRanges) + promoterTSSDistances[1],
    end(transcriptRanges) - promoterTSSDistances[2]
  )

  promoterEnds <- ifelse(
    strand(transcriptRanges) == "+",
    start(transcriptRanges) + promoterTSSDistances[2],
    end(transcriptRanges) - promoterTSSDistances[1]
  )

  promoterRanges <- GRanges(seqnames = seqnames(transcriptRanges), ranges = IRanges(promoterStarts, promoterEnds))

  regions[,c("promoter_overlap", "promoter_name")] <- annotation.annotateOverlap(regionRanges, promoterRanges, geneNames)

  return(regions)
}

annotation.annotateRepeats <- function (regions, regionRanges, repeats, classes) {

  relevantRepeats <- repeats[repClass %in% classes]

  repeatRanges <- GRanges(relevantRepeats$genoName, IRanges(relevantRepeats$genoStart, relevantRepeats$genoEnd))

  regions$num_repeats <- countOverlaps(regionRanges, repeatRanges)

  return(regions)
}

annotation.annotateRegions <- function (regionTable, gzippedCgiFile, gzippedGTFFile, gzippedRepeatMaskerAnnotationFile, allowedBiotypes, promoterTSSDistances) {

  regionTable$length <- regionTable$end - regionTable$start

  # remove any trailing 'chr's for range intersection
  regionRanges  <- GRanges(seqnames = str_remove(regionTable$chr, "^chr"), ranges = IRanges(regionTable$start, regionTable$end))

  if (annotation.isExistingFileOrNone(gzippedCgiFile)) {

    cgis <- fread(cmd = paste("zcat", gzippedCgiFile))
    cgis$chrom <- str_remove(cgis$chrom, "^chr")
    regionTable <- annotation.annotateCGIslands(regionTable, regionRanges, cgis)

  }

  if (annotation.isExistingFileOrNone(gzippedGTFFile)) {

    gtfRanges <- import(gzippedGTFFile)
    gtfMetaData <- values(gtfRanges)

    noHgncNameIndices <- is.na(gtfRanges$gene_name)
    gtfRanges$gene_name[noHgncNameIndices] <- gtfRanges$gene_id[noHgncNameIndices]

    biotypeColumn <- str_which(colnames(gtfMetaData), "gene_(bio)?type")

    gtfRanges <- gtfRanges[gtfMetaData[,biotypeColumn] %in% allowedBiotypes]
    gtfRanges <- renameSeqlevels(gtfRanges, str_remove(seqlevels(gtfRanges), "^chr"))

    regionTable <- annotation.annotateGenes(regionTable, regionRanges, gtfRanges)
    regionTable <- annotation.annotatePromoters(regionTable, regionRanges, gtfRanges, promoterTSSDistances)
  }

  if (annotation.isExistingFileOrNone(gzippedRepeatMaskerAnnotationFile)) {

    repeats <- fread(cmd = paste("zcat", gzippedRepeatMaskerAnnotationFile))
    repeats$genoName <- str_remove(repeats$genoName, "^chr")
    regionTable <- annotation.annotateRepeats(regionTable, regionRanges, repeats, c("LINE", "SINE", "LTR", "DNA"))

  }

  return(regionTable)
}
