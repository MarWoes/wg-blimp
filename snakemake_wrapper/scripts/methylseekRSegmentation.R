if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")
}

library(MethylSeekR)
library(BSgenome)
library(stringr)
library(data.table)
library(rtracklayer)
library(parallel)

### GLOBALS

FDR_CUTOFF <- 5
METHYLATION_CUTOFF <- 0.5
SEGMENTATION_WIDTH_PX <- 1366
SEGMENTATION_HEIGHT_PX <- 768
SEGMENTATION_FONT_SIZE <- 14

### FUNCTIONS
segmentIntoUMRsAndLMRs <- function(sample, methylationRanges, cgiRanges, numThreads, targetDir, genomeSeq, seqLengths, pmdSegments, sequencesStartingWithChr) {

  png(paste0(targetDir, "/fdr.stats.png"))
  fdrStats <- calculateFDRs(
    m = methylationRanges,
    CGIs = cgiRanges,
    PMDs = pmdSegments,
    num.cores = numThreads
  )
  dev.off()

  # TODO: use other values than the default one?
  selectedN <- as.integer(names(fdrStats$FDRs[as.character(METHYLATION_CUTOFF), ][fdrStats$FDRs[as.character(METHYLATION_CUTOFF), ] < FDR_CUTOFF])[1])

  png(paste0(targetDir, "/umr-lmr-heatmap.png"))
  segments <- segmentUMRsLMRs(
    m = methylationRanges,
    meth.cutoff = METHYLATION_CUTOFF,
    nCpG.cutoff = selectedN,
    PMDs = pmdSegments,
    num.cores = numThreads,
    myGenomeSeq = genomeSeq,
    seqLengths = seqLengths
  )
  dev.off()

  png(paste0(targetDir, "/umr-lmr-segmentation.png"), width = SEGMENTATION_WIDTH_PX, height = SEGMENTATION_HEIGHT_PX, pointsize = SEGMENTATION_FONT_SIZE)
  plotFinalSegmentation(
    m = methylationRanges,
    segs = segments,
    PMDs = pmdSegments,
    meth.cutoff = METHYLATION_CUTOFF
  )
  dev.off()

  umrLmrTable <- data.table(
    sample = sample,
    chr = as.character(seqnames(segments)),
    start = start(segments),
    end = end(segments),
    type = values(segments)$type,
    num_cpgs_filtered = values(segments)$nCG.segmentation,
    num_cpgs_ref = values(segments)$nCG,
    mean_methylation = values(segments)$pmeth,
    median_methylation = values(segments)$median.meth,
    pmd_included = !(length(pmdSegments) == 1 && is.na(pmdSegments))
  )

  if (!sequencesStartingWithChr) {
    umrLmrTable$chr <- str_remove(umrLmrTable$chr, "^chr")
  }

  fwrite(umrLmrTable, file = paste0(targetDir, "/umr-lmr.csv"))
}

loadOrInstall <- function(packageName) {

  isPackageInstalled <- packageName %in% rownames(installed.packages())

  if (!isPackageInstalled) {

    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repo = "http://cran.rstudio.com/")

    BiocManager::install(packageName)
  }

  library(packageName, character.only = TRUE)
}

### INPUT

# if (exists("snakemake")) {
#  save.image("methylseekr-debug.rds")
# }
fastaRef <- snakemake@input$ref
cgiAnnotationFile <- snakemake@input$cgi_annotation_file
geneAnnotationFile <- snakemake@input$gene_annotation_file
repeatMaskerAnnotationFile <- snakemake@input$repeat_masker_annotation_file
transcriptStartSiteFile <- snakemake@input$transcript_start_site_file
samples <- snakemake@params$samples
targetDir <- snakemake@params$target_dir
methylationDir <- snakemake@params$methylation_dir
calibrationChr <- snakemake@params$calibration_chr
targetGenomeName <- snakemake@params$genome_build
allowedBiotypes <- snakemake@params$biotypes
tssDistances <- snakemake@params$tss_distances
numThreads <- snakemake@threads
targetPmdFile <- snakemake@output$pmd_all
targetUmrLmrFile <- snakemake@output$umr_lmr_all

targetGenomeBS <- paste0("BSgenome.Hsapiens.UCSC.", targetGenomeName)

loadOrInstall(targetGenomeBS)
snakemake@source("regionAnnotation.R")

### SCRIPT

referenceBS <- get(targetGenomeBS)

# since own FASTA files may not be used with methylseekr, UCSC/GRCh 'chr's need to be taken care of
calibrationChr <- str_replace(calibrationChr, "^(?!chr)", "chr")

# query CGI's
session <- browserSession()
genome(session) <- targetGenomeName
query <- ucscTableQuery(session, "cpgIslandExt")
cgiRanges <- track(query)
genome(cgiRanges) <- NA
cgiRanges <- suppressWarnings(resize(cgiRanges, 5000, fix = "center"))

for (sample in samples) {

  methylationValues <- fread(
    paste0(methylationDir, "/", sample, "_CpG.bedGraph"),
    skip = 1L,
    header = FALSE,
    col.names = c("chr", "start", "end", "methPercent", "methylated", "unmethylated")
  )

  # restore 'non-chr' notation after computation is done
  sequencesStartingWithChr <- any(str_count(methylationValues$chr, "^chr")) > 0

  # adjust possibly faulty 'chr's here too
  methylationValues$chr <- str_replace(methylationValues$chr, "^(?!chr)", "chr")

  sampleRanges <- GRanges(methylationValues$chr, ranges = methylationValues$start)
  values(sampleRanges)$T <- methylationValues$methylated + methylationValues$unmethylated
  values(sampleRanges)$M <- methylationValues$methylated

  # only keep ranges that match the reference genome
  sampleRanges <- sampleRanges[as.logical(seqnames(sampleRanges) %in% seqnames(referenceBS))]

  png(paste0(targetDir, "/", sample, "/alphaCalibration.png"))
  pmdSegments <- segmentPMDs(
    m = sampleRanges,
    chr.sel = calibrationChr,
    num.cores = numThreads,
    seqLengths = seqlengths(referenceBS)
  )
  dev.off()

  png(paste0(targetDir, "/", sample, "/alphaPMDRemoved.png"))
  plotAlphaDistributionOneChr(
    m = subsetByOverlaps(sampleRanges, pmdSegments[values(pmdSegments)$type == "notPMD"]),
    chr.sel = calibrationChr,
    num.cores = numThreads
  )
  dev.off()

  png(paste0(targetDir, "/", sample, "/pmd.png"), width = SEGMENTATION_WIDTH_PX, height = SEGMENTATION_HEIGHT_PX, pointsize = SEGMENTATION_FONT_SIZE)
  plotPMDSegmentation(
    m = sampleRanges,
    segs = pmdSegments
  )
  dev.off()

  pmdSegmentTable <- data.table(
    sample = sample,
    chr = as.character(seqnames(pmdSegments)),
    type = values(pmdSegments)$type,
    start = start(pmdSegments),
    end = end(pmdSegments),
    num_cpgs = values(pmdSegments)$nCG
  )
  pmdSegmentTable <- pmdSegmentTable[type == "PMD"]

  if (!sequencesStartingWithChr) {
    pmdSegmentTable$chr <- str_remove(pmdSegmentTable$chr, "^chr")
  }

  fwrite(pmdSegmentTable, file = paste0(targetDir, "/", sample, "/pmd-segments.csv"))

  segmentIntoUMRsAndLMRs(
    sample = sample,
    methylationRanges = sampleRanges,
    cgiRanges = cgiRanges,
    numThreads = numThreads,
    targetDir = paste0(targetDir, "/", sample, "/LMRUMRwithPMD/"),
    genomeSeq = referenceBS,
    seqLengths = seqlengths(referenceBS),
    pmdSegments = pmdSegments,
    sequencesStartingWithChr
  )

  segmentIntoUMRsAndLMRs(
    sample = sample,
    methylationRanges = sampleRanges,
    cgiRanges = cgiRanges,
    numThreads = numThreads,
    targetDir = paste0(targetDir, "/", sample, "/LMRUMRwithoutPMD/"),
    genomeSeq = referenceBS,
    seqLengths = seqlengths(referenceBS),
    pmdSegments = NA,
    sequencesStartingWithChr
  )

}

pmdFiles <- list.files(targetDir, pattern = "pmd-segments.csv", recursive = TRUE, full.names = TRUE)
umrLmrFiles <- list.files(targetDir, pattern = "umr-lmr.csv", recursive = TRUE, full.names = TRUE)

pmdContents <- lapply(pmdFiles, fread)
umrLmrContents <- lapply(umrLmrFiles, fread)

combinedPmdTable <- rbindlist(pmdContents)
combinedUmrLmrTable <- rbindlist(umrLmrContents)

annotatedPmdTable <- annotation.annotateRegions(
  combinedPmdTable,
  cgiAnnotationFile,
  geneAnnotationFile,
  repeatMaskerAnnotationFile,
  transcriptStartSiteFile,
  allowedBiotypes,
  tssDistances
)

annotatedUmrLmrTable <- annotation.annotateRegions(
  combinedUmrLmrTable,
  cgiAnnotationFile,
  geneAnnotationFile,
  repeatMaskerAnnotationFile,
  transcriptStartSiteFile,
  allowedBiotypes,
  tssDistances
)

fwrite(annotatedPmdTable, targetPmdFile)
fwrite(annotatedUmrLmrTable, targetUmrLmrFile)
