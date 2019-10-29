# load prism for yaml config highlighting, see
# https://stackoverflow.com/questions/47445260/how-to-enable-syntax-highlighting-in-r-shiny-app-with-htmloutput
prismDependencies <- tags$head(
  tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/prism/1.15.0/prism.min.js"),
  tags$link(rel = "stylesheet", type = "text/css",
            href = "https://cdnjs.cloudflare.com/ajax/libs/prism/1.15.0/themes/prism.min.css")
)

prismYamlDefinition <- tags$head(
  tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/prism/1.15.0/components/prism-yaml.min.js")
)

shinyUI(
  dashboardPage(
    dashboardHeader(title = "WGBS browser"),
    dashboardSidebar(
      sidebarMenuOutput("sidebarTabs")
    ),
    dashboardBody(
      prismDependencies,
      prismYamlDefinition,
      tabItems(
        tabItem(tabName = "datasetTab", uiOutput("datasetSelection")),
        tabItem(tabName = "summaryTab",
          fluidRow(
            box(title = "Summary", tableOutput("summaryTable"), footer = "This box shows some basic alignment statistics for each sample in the dataset."),
            box(title = "MultiQC report", uiOutput("multiqcLink"))
          ),
          fluidRow(
            box(
              title = "IGV links",
              tableOutput("bamIgvLinks"),
              footer = "You may copy/paste the sample's BAM links into IGV (File -> Load From URL) or use the direct IGV link."
            ),
            box(
              title = "Qualimap reports",
              uiOutput("qualimapLinks"),
              footer = "Detailed Qualimap reports for all samples are available here."
            )
          )
        ),
        tabItem(tabName = "parameterTab",
          fluidRow(
            box(
              title = "Pipeline configuration",
              width = 12,
              "This box shows the snakemake configuration used for the pipeline run for this dataset.",
              htmlOutput("pipelineConfigCode")
            )
          )
        ),
        tabItem(tabName = "dmrTab",
          fluidRow(
            box(
              title = "Filters (CpGs, difference, length, repeats)",
              width = 3,
              sliderInput("cpgFilter", label = "Number of CpGs:", value = c(1, 100), min = 1, max = 100),
              sliderInput("diffFilter", label = "Absolute diff.:", value = c(0, 1), min = 0, max = 1, step = 0.01),
              sliderInput("lengthFilter",  label = "DMR length:", value = c(0, 1000000), min = 0, max = 1000000)
            ),
            box(
              title = "Filters (repeats, q-values, coverage)",
              width = 3,
              sliderInput("repeatFilter",  label = "Number of repetitive elements:", value = c(0, 10), min = 0, max = 10),
              sliderInput("qFilter", label = "q value (metilene):", value = c(0,1), min = 0, max = 1, step = 0.01),
              sliderInput("covFilter", label = "Minimum mean coverage:", value = 0, min = 0, max = 30, step = 1)
            ),
            box(
              title = "Filters (genes, tools, CG islands)",
              width = 3,
              selectizeInput("geneFilter", label = "Neighbouring genes:", choices = "None", multiple = TRUE),
              selectizeInput("toolFilter", label = "Include calls by tools:", choices = "None", multiple = TRUE),
              checkboxInput("onlyCgiFilter", "Only DMRs in CG islands", value = FALSE)
            ),
            box(
              width = 3,
              title = "Download",
              downloadButton("dmrDownload", label = "Download filtered DMRs", icon = icon("download")),
              footer = "Note: Only DMRs passing the filters are provided in the output"
            )
          ),
          fluidRow(
            box(
              title = "DMR table",
              width = 12,
              DT::dataTableOutput("dmrTable")
            )
          ),
          fluidRow(
            box(
              title = "Number of CpGs",
              width = 3,
              plotOutput("numCpGHist")
            ),
            box(
              title = "Average difference",
              width = 3,
              plotOutput("diffHist")
            ),
            box(
              title = "DMR length",
              width = 3,
              plotOutput("lengthHist")
            ),
            box(
              title = "Tools",
              width = 3,
              plotOutput("toolHist")
            )
          ),
          fluidRow(
            box(
              title = "CG island overlap",
              width = 3,
              plotOutput("cgiHist")
            ),
            box(
              title = "Number of repeats",
              width = 3,
              plotOutput("repeatHist")
            ),
            box(
              title = "q-value (metilene)",
              width = 3,
              plotOutput("qHist")
            ),
            box(
              title = "mean coverage (mapq >= 10)",
              width = 3,
              plotOutput("covHist")
            )
          )
        ),
        tabItem(tabName = "segmentationTab",
          fluidRow(
            box(
              title = "Select sample",
              width = 4,
              selectInput("segmentationSampleSelect", "Show segmentation for sample:", choices = NULL, selectize = FALSE),
              checkboxInput("segmentationWithPMD", "Include PMDs in computation", value = FALSE)
            ),
            box(
              title = "Download UMRs/LMRs",
              width = 4,
              downloadButton("segmentationUmrLmrDownload", "Download UMR/LMR table", icon = icon("download"))
            ),
            conditionalPanel(
              condition = "input.segmentationWithPMD",
              box(
                title = "Download PMDs",
                width = 4,
                downloadButton("segmentationPmdDownload", "Download PMDs", icon = icon("download"))
              )
            )
          ),
          conditionalPanel(
            "input.segmentationWithPMD",
            fluidRow(
              box(
                title = "Posterior mean of alpha (PMDs removed)",
                width = 4,
                imageOutput("segmentationPosteriorPMDRemovedAlphaImage", height = "480px")
              ),
              box(
                title = "UMR/LMR heatmap with PMDs",
                width = 4,
                imageOutput("segmentationCpgMedianMethylationWithPMD", height = "480px")
              ),
              box(
                title = "FDR stats with PMDs",
                width = 4,
                imageOutput("segmentationFdrStatsWithPMD", height = "480px")
              )
            ),
            fluidRow(
              box(
                title = "PMDs",
                width = 12,
                DT::dataTableOutput("segmentationPmdTable")
              )
            )
          ),
          conditionalPanel(
            "!input.segmentationWithPMD",
            fluidRow(
              box(
                title = "Posterior mean of alpha",
                width = 4,
                imageOutput("segmentationPosteriorAlphaImage", height = "480px")
              ),
              box(
                title = "UMR/LMR heatmap without PMDs",
                width = 4,
                imageOutput("segmentationCpgMedianMethylationWithoutPMD", height = "480px")
              ),
              box(
                title = "FDR stats without PMDs",
                width = 4,
                imageOutput("segmentationFdrStatsWithoutPMD", height = "480px")
              )
            )
          ),
          fluidRow(
            box(
              title = "UMRs/LMRs",
              width = 12,
              DT::dataTableOutput("segmentationUmrLmrTable")
            )
          )
        )
      )
    )
  )
)
