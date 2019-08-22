library(MethylSeekR)
library(BSgenome)
library(stringr)
library(data.table)
library(rtracklayer)
library(parallel)
# if (exists("snakemake")) {
#   save.image("methylseekr-debug.rds")
# }
fastaRef <- snakemake@input$ref
samples <- snakemake@params$samples
targetDir <- snakemake@params$target_dir
methylationDir <- snakemake@params$methylation_dir
calibrationChr <- snakemake@params$calibration_chr
cgiTargetGenome <- snakemake@params$genome_build


numThreads <- snakemake@threads

sapply(paste0(targetDir, "/", samples), dir.create, showWarnings = FALSE)

referenceBS <- readDNAStringSet(fastaRef)

# according to FASTA defaults, only the ID before the first whitespace is considered. Also remove any 'chr's
names(referenceBS) <- str_remove(names(referenceBS), "\\s.*")
names(referenceBS) <- str_remove(names(referenceBS), "^chr")
calibrationChr <- str_remove(calibrationChr, "^chr")

# query CGI's
session <- browserSession()
genome(session) <- cgiTargetGenome
query <- ucscTableQuery(session, "cpgIslandExt")
cgiRanges <- track(query)
genome(cgiRanges) <- NA
cgiRanges <- suppressWarnings(resize(cgiRanges, 5000, fix = "center"))
cgiRanges <- renameSeqlevels(cgiRanges, str_remove(as.character(seqlevels(cgiRanges)), "^chr"))

for (sample in samples) {
  
  methylationValues <- fread(
    paste0(methylationDir, "/", sample, "_CpG.bedGraph"),
    skip = 1L,
    header = FALSE,
    col.names = c("chr", "start", "end", "methPercent", "methylated", "unmethylated")
  )
  
  sampleRanges <- GRanges(methylationValues$chr, ranges = methylationValues$start)
  values(sampleRanges)$T <- methylationValues$methylated + methylationValues$unmethylated
  values(sampleRanges)$M <- methylationValues$methylated
  
  pmdSegments <- segmentPMDs(
    m = sampleRanges,
    chr.sel = calibrationChr,
    num.cores = numThreads,
    seqLengths = seqlengths(referenceBS),
    pdfFilename = paste0(targetDir, "/", sample, "/alpha.", calibrationChr, ".pdf")
  )

  plotAlphaDistributionOneChr(
    m = subsetByOverlaps(sampleRanges, pmdSegments[values(pmdSegments)$type == "notPMD"]),
    chr.sel = calibrationChr,
    num.cores = numThreads,
    pdfFilename = paste0(targetDir, "/", sample, "/alpha.pmd.removed.", calibrationChr, ".pdf")
  )
  
  plotPMDSegmentation(
    m = sampleRanges,
    segs = pmdSegments,
    pdfFilename =  paste0(targetDir, "/", sample, "/pmd.pdf")
  )
  
  savePMDSegments(pmdSegments, TableFilename = paste0(targetDir, "/", sample, "/pmd-segments.tab"))
  
  fdrStats <- calculateFDRs(
    m = sampleRanges,
    CGIs = cgiRanges,
    PMDs = pmdSegments,
    num.cores = numThreads,
    pdfFilename = paste0(targetDir, "/", sample, "/fdr.stats.pdf")
  )
  
  # TODO: replace with non-default values?
  fdrCutoff <- 5
  selectedMethylation <- 0.5
  selectedNumber <- as.integer(names(fdrStats$FDRs[as.character(selectedMethylation), ][fdrStats$FDRs[as.character(selectedMethylation), ] < fdrCutoff])[1])
  

}


file.create(snakemake@output$alpha_distribution)
