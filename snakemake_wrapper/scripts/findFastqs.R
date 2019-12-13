if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")
}

library(data.table)

### FUNCTIONS

guessFastqFiles <- function (rawDir, forwardReadRegex, reverseReadRegex, targetSample, allSamples) {

  availableForwardReads <- list.files(rawDir, pattern = forwardReadRegex, full.names = TRUE, recursive = TRUE)
  availableReverseReads <- list.files(rawDir, pattern = reverseReadRegex, full.names = TRUE, recursive = TRUE)

  # prevent overlapping names being assigned wrong, e.g. 'sample11_1.fq.gz' being assigned to 'sample1'
  superStringSamples <- allSamples[allSamples != targetSample & grepl(targetSample, allSamples, fixed = TRUE)]

  for (superStringSample in superStringSamples) {

    availableForwardReads <- availableForwardReads[!grepl(superStringSample, basename(availableForwardReads), fixed = TRUE)]
    availableReverseReads <- availableReverseReads[!grepl(superStringSample, basename(availableReverseReads), fixed = TRUE)]

  }

  forwardReadMatches <- grepl(targetSample, basename(availableForwardReads), fixed = TRUE)
  reverseReadMatches <- grepl(targetSample, basename(availableReverseReads), fixed = TRUE)


  # sort read files to prevent any mismatching forward/reverse read situations
  sortedForwardReads <- sort(availableForwardReads[forwardReadMatches])
  sortedReverseReads <- sort(availableReverseReads[reverseReadMatches])

  if (length(sortedForwardReads) == 0 || length(sortedReverseReads) == 0) {
    stop(paste0("No forward or reverse reads found for sample ", targetSample, "."))
  }

  if (length(sortedForwardReads) != length(sortedReverseReads)) {
    stop(paste0("Number of forward and reverse read files differs for sample ", targetSample, "."))
  }

  return(list(
    forward = sortedForwardReads,
    reverse = sortedReverseReads
  ))
}

readFastqsFromCsv <- function (csvFile, targetSample) {

  fastqMap <- fread(csvFile)

  targetSampleReadFiles <- fastqMap[sample == targetSample]

  forwardFiles <- targetSampleReadFiles$forward
  reverseFiles <- targetSampleReadFiles$reverse

  for (file in c(forwardFiles, reverseFiles)) {

    if (!file.exists(file))
      stop(paste0("Could not find file ", file))

  }

  return(list(
    forward = forwardFiles,
    reverse = reverseFiles
  ))

}

### SCRIPT

rawDir <- snakemake@params$rawdir
allSamples <- snakemake@params$samples
targetSample <- snakemake@wildcards$sample

forwardOutputFile <- snakemake@output$first
reverseOutputFile <- snakemake@output$second

forwardReadRegex <- snakemake@params$first_regex
reverseReadRegex <- snakemake@params$second_regex

sampleFastqCsv <- snakemake@params$sample_fastq_csv

if (length(sampleFastqCsv) == 1 && is.character(sampleFastqCsv) && file.exists(sampleFastqCsv)) {

  sampleFastqFiles <- readFastqsFromCsv(sampleFastqCsv, targetSample)

} else {

  sampleFastqFiles <- guessFastqFiles(rawDir, forwardReadRegex, reverseReadRegex, targetSample, allSamples)

}

writeLines(sampleFastqFiles$forward, forwardOutputFile)
writeLines(sampleFastqFiles$reverse, reverseOutputFile)