if (exists("snakemake")) {
  logFile <- file(snakemake@log[[1]])

  sink(logFile, append = TRUE)
  sink(logFile, append = TRUE, type = "message")
}

library(tidyverse)
library(data.table)

benchmarkDir <- snakemake@params$benchmark_dir
targetPlot <- snakemake@output$runtime_report

benchmarkFiles <- list.files(benchmarkDir, recursive = TRUE, pattern = ".tsv")

benchmarkReports <- lapply(benchmarkFiles, function (file) {

  task <- str_remove_all(file, fixed(".tsv"))
  fullPath <- paste0(benchmarkDir, "/", file)

  benchmarkReport <- fread(fullPath, na.strings = "-")
  benchmarkReport$task <- task

  return(benchmarkReport)
})

fullBenchmarkReport <- rbindlist(benchmarkReports, fill = TRUE)

pdf(targetPlot)

ggplot(fullBenchmarkReport, aes(x = task, y = s)) +
  ggtitle("Task run times in seconds") +
  xlab("task") +
  ylab("time in seconds") +
  geom_bar(stat = "identity") +
  coord_flip()

ggplot(fullBenchmarkReport, aes(x = task, y = max_rss)) +
  ggtitle("Maximum task RAM usage") +
  xlab("task") +
  ylab("maximum RAM usage in megabytes") +
  geom_bar(stat = "identity") +
  coord_flip()


dev.off()
