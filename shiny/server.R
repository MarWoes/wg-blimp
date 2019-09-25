
shiny.wgbs.datasetsFile <- ifelse(exists("projectsFileArgument"), projectsFileArgument, "projects.txt")
shiny.wgbs.multiqcPath <- ifelse(exists("multiqcDirArgument"), multiqcDirArgument, "multiqc")

shiny.wgbs.datasets <- shiny.wgbs.loadDatasets(shiny.wgbs.datasetsFile)

shiny.wgbs.serveMultiqcReports(shiny.wgbs.multiqcPath, shiny.wgbs.datasets)
shiny.wgbs.serveQualimapReports(shiny.wgbs.datasets)
shiny.wgbs.serveBamFiles(shiny.wgbs.datasets)

shiny.wgbs.getIgvLocusName <- function (reactiveTableHandle) {

  dmrTable <- reactiveTableHandle()

  chr <- dmrTable$chr
  regionStart <- dmrTable$start
  regionEnd <- dmrTable$end

  return(paste(chr, ":", regionStart, "-", regionEnd, sep = ""))
}

shiny.wgbs.getIgvCommandLinks <- function (commandName, arguments) {

  igvLinkAddress <- paste0("http://localhost:60151/", commandName, "?", arguments)

  # link to local igv running, but do nothing when clicked (except performing igv http request)
  links <- paste("<a href='", igvLinkAddress,  "' target='_blank'>IGV</a>", sep = "")

  return(links)
}

shiny.wgbs.plotHistogram <- function(reactiveTableHandle, columnName, stat = "bin") {

  req(reactiveTableHandle())

  ggplot(reactiveTableHandle(), aes_string(x = columnName)) +
    geom_histogram(bins = 50, fill = "steelblue1", colour = "steelblue4", stat = stat)

}

shiny.wgbs.getPrecomputedSegmentationPlot <- function (input, selectedDataset, subPath) {

  req(selectedDataset())
  req(input$segmentationSampleSelect)

  imagePath <- paste(shiny.wgbs.datasets[[selectedDataset()]]$segmentationDirectory, input$segmentationSampleSelect, subPath, sep = "/")

  list(
    src = imagePath
  )
}

shinyServer(function(input, output, session) {

  selectedDataset <- reactiveVal()
  dmrTable <- reactiveVal()

  observe({

    req(selectedDataset())

    progress <- shiny::Progress$new()
    progress$set(message = "Loading dataset...", value = 0.1)

    selectedTable <- shiny.wgbs.datasets[[selectedDataset()]]$dmrs

    dmrTable(selectedTable)

    rawGeneList <- as.character(unique(c(selectedTable$gene_name, selectedTable$promoter_name)))

    genesInTable <- unlist(strsplit(rawGeneList, ","))

    updateSelectizeInput(session, "geneFilter", choices = genesInTable, selected = NULL)
    updateSelectizeInput(session, "toolFilter", choices = unique(selectedTable$tool), selected = NULL)

    minCpG <- min(selectedTable$num_cpgs)
    maxCpG <- max(selectedTable$num_cpgs)

    minLength <- min(selectedTable$length)
    maxLength <- max(selectedTable$length)

    maxNumberOfRepeats <- max(selectedTable$num_repeats)

    updateSliderInput(session, "lengthFilter", min = minLength, max = maxLength, value = c(minLength, maxLength))
    updateSliderInput(session, "cpgFilter", min = minCpG, max = maxCpG, value = c(minCpG, maxCpG))
    updateSliderInput(session, "cpgFilter", min = minCpG, max = maxCpG, value = c(minCpG, maxCpG))
    updateSliderInput(session, "diffFilter", min = 0, max = 1, value = c(0,1))
    updateSliderInput(session, "repeatFilter", min = 0, max = maxNumberOfRepeats, value = c(0, maxNumberOfRepeats))
    updateSliderInput(session, "qFilter", value = c(0,1))

    updateCheckboxInput(session, "onlyCgiFilter", value = FALSE)

    progress$close()

  })

  observe({

    req(selectedDataset())

    selectedTable <- shiny.wgbs.datasets[[selectedDataset()]]$dmrs

    req(selectedTable)

    genePattern <- paste0("\\b", input$geneFilter, "\\b", collapse = "|")
    toolPattern <- paste0("\\b", input$toolFilter, "\\b", collapse = "|")

    selectedTable <- selectedTable[
      input$cpgFilter[1] <= num_cpgs & num_cpgs <= input$cpgFilter[2] &
      input$diffFilter[1] <= abs(diff) & abs(diff) <= input$diffFilter[2] &
      input$lengthFilter[1] <= length & length <= input$lengthFilter[2] &
      input$repeatFilter[1] <= num_repeats & num_repeats <= input$repeatFilter[2] &
      input$covFilter <= mean_cov &
      (is.na(qValue) | input$qFilter[1] <= qValue & qValue <= input$qFilter[2]) &
      (is.null(input$geneFilter) | seq_len(nrow(selectedTable)) %in% grep(genePattern, paste(gene_name, promoter_name))) &
      (is.null(input$toolFilter) | seq_len(nrow(selectedTable)) %in% grep(toolPattern, tool)) &
      (!input$onlyCgiFilter | `cgi_overlap`)
    ]

    dmrTable(selectedTable)
  })

  output$datasetSelection <- renderUI({

    lapply(names(shiny.wgbs.datasets), function (datasetName) {

      selectionButtonId <- paste(datasetName, "SelectButton", sep = "")
      selectionButton <- actionButton(selectionButtonId, label = "Select", class = "btn btn-sm")

      numberOfSamples <- nrow(shiny.wgbs.datasets[[datasetName]]$summary)
      infoBoxLabel <- paste(numberOfSamples, "samples")

      return(infoBox(datasetName, list(selectionButton, infoBoxLabel), icon = icon("th-list")))

    })

  })

  observers <- sapply(names(shiny.wgbs.datasets), function (datasetName) {

    buttonId <- paste(datasetName, "SelectButton", sep = "")

    observeEvent({input[[buttonId]]}, {

      selectedDataset(datasetName)
    })

  })

  output$selectedDataset <- renderText({
    selectedDataset()
  })

  output$summaryTable <- renderTable({

    req(selectedDataset())

    shiny.wgbs.datasets[[selectedDataset()]]$summary
  }, striped = TRUE, bordered = TRUE)

  output$bamIgvLinks <- renderTable({

    req(selectedDataset())

    # shamelessly stolen from https://stackoverflow.com/questions/45943034/in-a-deployed-shinyapp-how-to-get-the-name-of-the-account
    hostname <- session$clientData$url_hostname
    port <- session$clientData$url_port

    url <- paste0(session$clientData$url_protocol,
                  "//",
                  hostname,
                  if(port != "") paste0(":", port))

    sampleNames <- shiny.wgbs.datasets[[selectedDataset()]]$summary$sample
    bamUrls <- paste0(url, "/", selectedDataset(), "/", sampleNames, ".bam")

    igvLoadingLinks <- shiny.wgbs.getIgvCommandLinks("load", paste0("file=",encodeURIComponent(bamUrls)))

    linkTable <- data.table(
      sample = sampleNames,
      URL = bamUrls,
      `IGV link` = igvLoadingLinks
    )

    return(linkTable)
  }, bordered = TRUE, striped = TRUE, sanitize.text.function = function(x) x)

  output$multiqcLink <- renderUI({

    req(selectedDataset())

    list(
      a("Full MultiQC report", href = paste("multiqc/", selectedDataset(), ".html", sep = ""), class = "btn btn-primary", target = "_blank"),
      hr(),
      div("For every project, a multiqc report is created. It contains detailed information including sequencing and alignment qualities and coverages and GC content.")
    )
  })

  output$qualimapLinks <- renderUI({

    req(selectedDataset())

    samples <- shiny.wgbs.datasets[[selectedDataset()]]$summary$sample

    lapply(samples, function (sample) div(a(sample, href = paste0("qualimap-", selectedDataset(), "/", sample, "/qualimapReport.html"), target = "_blank")))

  })

  output$dmrTable <- DT::renderDataTable({

    req(selectedDataset())
    req(dmrTable())

    dmrTableWithLinks <- dmrTable()

    dmrIgvLoci <- shiny.wgbs.getIgvLocusName(dmrTable)
    dmrTableWithLinks$igv <- shiny.wgbs.getIgvCommandLinks("goto", paste0("locus=", dmrIgvLoci))

    return(DT::datatable(dmrTableWithLinks,
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("compact", "hover", "stripe")),
                         escape = FALSE
          ) %>% formatRound(columns = "diff")
    )
  })

  output$dmrDownload <- downloadHandler(
    filename = 'dmrs.tsv',
    content = function (file) {
      write.table(dmrTable(), file = file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    }
  )

  output$pipelineConfigCode <- renderUI({

    req(selectedDataset())

    tagList(
      pre(code(shiny.wgbs.datasets[[selectedDataset()]]$config, class = "language-yaml")),
      tags$script("Prism.highlightAll()")
    )

  })

  observe({

    req(selectedDataset())

    availableSamples <- shiny.wgbs.datasets[[selectedDataset()]]$summary$sample

    updateSelectInput(session, "segmentationSampleSelect", choices = availableSamples)

  })

  output$segmentationUmrLmrTable <- DT::renderDataTable({

    req(selectedDataset())
    req(input$segmentationSampleSelect)

    umrLmrTable <- shiny.wgbs.datasets[[selectedDataset()]]$umrLmrAll

    req(umrLmrTable)

    subsetUmrLmrTable <- umrLmrTable[sample == input$segmentationSampleSelect & (pmd_included == input$segmentationWithPMD)]

    return(DT::datatable(subsetUmrLmrTable,
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("compact", "hover", "stripe")),
                         escape = FALSE
    ) %>% formatRound(columns = c("mean_methylation", "median_methylation")))
  })

  output$segmentationPmdTable <- DT::renderDataTable({

    req(selectedDataset())
    req(input$segmentationSampleSelect)
    req(input$segmentationWithPMD)

    allPmdTable <- shiny.wgbs.datasets[[selectedDataset()]]$pmd

    req(allPmdTable)

    subsetPmdTable <- allPmdTable[sample == input$segmentationSampleSelect]

    return(DT::datatable(subsetPmdTable,
                         selection = "single",
                         style = "bootstrap",
                         class = DT:::DT2BSClass(c("compact", "hover", "stripe")),
                         escape = FALSE
    ))
  })

  output$segmentationPosteriorAlphaImage <- renderImage({ shiny.wgbs.getPrecomputedSegmentationPlot(input, selectedDataset, "alphaCalibration.png") }, deleteFile = FALSE)
  output$segmentationCpgMedianMethylationWithPMD <- renderImage({ shiny.wgbs.getPrecomputedSegmentationPlot(input, selectedDataset, "LMRUMRwithPMD/umr-lmr-heatmap.png") }, deleteFile = FALSE)
  output$segmentationFdrStatsWithPMD <- renderImage({ shiny.wgbs.getPrecomputedSegmentationPlot(input, selectedDataset, "LMRUMRwithPMD/fdr.stats.png") }, deleteFile = FALSE)

  output$segmentationPosteriorPMDRemovedAlphaImage <- renderImage({ shiny.wgbs.getPrecomputedSegmentationPlot(input, selectedDataset, "alphaPMDRemoved.png") }, deleteFile = FALSE)
  output$segmentationCpgMedianMethylationWithoutPMD <- renderImage({ shiny.wgbs.getPrecomputedSegmentationPlot(input, selectedDataset, "LMRUMRwithoutPMD/umr-lmr-heatmap.png") }, deleteFile = FALSE)
  output$segmentationFdrStatsWithoutPMD <- renderImage({ shiny.wgbs.getPrecomputedSegmentationPlot(input, selectedDataset, "LMRUMRwithoutPMD/fdr.stats.png") }, deleteFile = FALSE)

  output$numCpGHist <- renderPlot({ shiny.wgbs.plotHistogram(dmrTable, "num_cpgs") })
  output$diffHist <- renderPlot({ shiny.wgbs.plotHistogram(dmrTable, "diff") })
  output$lengthHist <- renderPlot({ shiny.wgbs.plotHistogram(dmrTable, "length") })
  output$toolHist <- renderPlot({ shiny.wgbs.plotHistogram(dmrTable, "tool", stat = "count") })
  output$cgiHist <- renderPlot({ shiny.wgbs.plotHistogram(dmrTable, "cgi_overlap", stat = "count") })
  output$repeatHist <- renderPlot({ shiny.wgbs.plotHistogram(dmrTable, "num_repeats") })
  output$qHist <- renderPlot({ shiny.wgbs.plotHistogram(dmrTable, "qValue") })
  output$covHist <- renderPlot({ shiny.wgbs.plotHistogram(dmrTable, "mean_cov") })

})
