library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)

wgbs.extractSequences <- function (cgis, indexFn) {
  
  simulatedSequences <- apply(cgis, 2, function (cgi) {
    
    randomSequence <- BSgenome.Hsapiens.UCSC.hg38[[cgi$chrom]][indexFn(cgi$chromStart, cgi$chromEnd)]
    
    return(as.character(randomSequence))
  })
  
  return(simulatedSequences)
}

wgbs.saveSequencesAsFasta <- function (simulatedSequences, sequenceNames, fileName) {
  
  fastaIdentifier <- paste(">", sequenceNames, sep = "")
  
  combinedFasta <- c(rbind(fastaIdentifier, simulatedSequences))
  
  writeLines(combinedFasta, fileName)
  
}

wgbs.simulateReference <- function (cgiFile, dmrFastaFile, flankingFastaFile, chosenCgiFile, minSize = 1000, healthy_flanking_regions = 50000) {
  
  cgis <- fread(cgiFile)
  
  chromosomesOfInterest <- paste("chr", c(1:21, "X", "Y"), sep = "")
  
  randomCgis <- data.frame(sapply(chromosomesOfInterest, function (chr) {
    
    chromosomeCandidates <- cgis[chrom == chr & length >= minSize]
    
    randomCgi <- chromosomeCandidates[sample(seq_len(nrow(chromosomeCandidates)), size = 1)]
    
    return(randomCgi)
  }))
  
  simulatedDmrSequences <- wgbs.extractSequences(randomCgis, function (start, end) start:end)
  simulatedFlankingSequences <- wgbs.extractSequences(randomCgis, function (start, end) c((start - healthy_flanking_regions):start, end:(end + healthy_flanking_regions)))
  
  wgbs.saveSequencesAsFasta(simulatedDmrSequences, chromosomesOfInterest, dmrFastaFile)
  wgbs.saveSequencesAsFasta(simulatedFlankingSequences, chromosomesOfInterest, flankingFastaFile)
  
  write.table(t(randomCgis), file = chosenCgiFile, col.names = TRUE, row.names = FALSE, sep = ";")
  
  return(simulatedDmrSequences)
}