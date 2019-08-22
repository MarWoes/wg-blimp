library(MethylSeekR)
library(BSgenome)
library(stringr)
library(data.table)
library(rtracklayer)
library(parallel)


### FUNCTIONS

segmentIntoUMRsAndLMRs <- function(methylationRanges, cgiRanges, numThreads, targetDir, genomeSeq, seqLengths, pmdSegments) {

  fdrStats <- calculateFDRs(
    m = methylationRanges,
    CGIs = cgiRanges,
    PMDs = pmdSegments,
    num.cores = numThreads,
    pdfFilename = paste0(targetDir, "/fdr.stats.pdf")
  )

  # TODO: use other values than the default one?
  fdrCutoff <- 5
  methylationCutoff <- 0.5
  selectedN <- as.integer(names(fdrStats$FDRs[as.character(methylationCutoff), ][fdrStats$FDRs[as.character(methylationCutoff), ] < fdrCutoff])[1])

  segments <- segmentUMRsLMRs(
    m = methylationRanges,
    meth.cutoff = methylationCutoff,
    nCpG.cutoff = selectedN,
    PMDs = pmdSegments,
    num.cores = numThreads,
    myGenomeSeq = genomeSeq,
    seqLengths = seqLengths,
    pdfFilename = paste0(targetDir, "/umr-lmr-heatmap.pdf")
  )

  plotFinalSegmentation(
    m = methylationRanges,
    segs = segments,
    PMDs = pmdSegments,
    meth.cutoff = methylationCutoff,
    pdfFilename = paste0(targetDir, "/umr-lmr-segmentation.pdf")
  )

  saveUMRLMRSegments(
    segments,
    TableFilename = paste0(targetDir, "/umr-lmr.tab")
  )
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
#   save.image("methylseekr-debug.rds")
# }
fastaRef <- snakemake@input$ref
samples <- snakemake@params$samples
targetDir <- snakemake@params$target_dir
methylationDir <- snakemake@params$methylation_dir
calibrationChr <- snakemake@params$calibration_chr
targetGenomeName <- snakemake@params$genome_build
numThreads <- snakemake@threads

### SCRIPT

targetGenomeBS <- paste0("BSgenome.Hsapiens.UCSC.", targetGenomeName)

loadOrInstall(targetGenomeBS)

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

  # adjust possibly faulty 'chr's here too
  methylationValues$chr <- str_replace(methylationValues$chr, "^(?!chr)", "chr")

  sampleRanges <- GRanges(methylationValues$chr, ranges = methylationValues$start)
  values(sampleRanges)$T <- methylationValues$methylated + methylationValues$unmethylated
  values(sampleRanges)$M <- methylationValues$methylated

  pmdSegments <- segmentPMDs(
    m = sampleRanges,
    chr.sel = calibrationChr,
    num.cores = numThreads,
    seqLengths = seqlengths(referenceBS),
    pdfFilename = paste0(targetDir, "/", sample, "/alphaCalibration.pdf")
  )

  plotAlphaDistributionOneChr(
    m = subsetByOverlaps(sampleRanges, pmdSegments[values(pmdSegments)$type == "notPMD"]),
    chr.sel = calibrationChr,
    num.cores = numThreads,
    pdfFilename = paste0(targetDir, "/", sample, "/alphaPMDRemoved.pdf")
  )

  plotPMDSegmentation(
    m = sampleRanges,
    segs = pmdSegments,
    pdfFilename =  paste0(targetDir, "/", sample, "/pmd.pdf")
  )

  savePMDSegments(pmdSegments, TableFilename = paste0(targetDir, "/", sample, "/pmd-segments.tab"))

  segmentIntoUMRsAndLMRs(
    methylationRanges = sampleRanges,
    cgiRanges = cgiRanges,
    numThreads = numThreads,
    targetDir = paste0(targetDir, "/", sample, "/LMRUMRwithPMD/"),
    genomeSeq = referenceBS,
    seqLengths = seqlengths(referenceBS),
    pmdSegments = pmdSegments
  )

  segmentIntoUMRsAndLMRs(
    methylationRanges = sampleRanges,
    cgiRanges = cgiRanges,
    numThreads = numThreads,
    targetDir = paste0(targetDir, "/", sample, "/LMRUMRwithoutPMD/"),
    genomeSeq = referenceBS,
    seqLengths = seqlengths(referenceBS),
    pmdSegments = NA
  )

}


file.create(snakemake@output$alpha_distribution)
