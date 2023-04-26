#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Libraries ----
library(shiny)
library(shinymeta)
library(purrr)
require(ggplot2)
require(GenomicRanges)
require(UpSetR)
require(eulerr)
require(ggplotify)

# Default values ----
colorList <- 
  list(
    c("#1B9E77", "#D95F02", "#bdbdbd"),
    RColorBrewer::brewer.pal(9, "Set1"),
    RColorBrewer::brewer.pal(8, "Set2"),
    RColorBrewer::brewer.pal(12, "Set3"),
    RColorBrewer::brewer.pal(9, "Pastel1"),
    RColorBrewer::brewer.pal(8, "Pastel2"),
    RColorBrewer::brewer.pal(12, "Paired"),
    RColorBrewer::brewer.pal(8, "Dark2"),
    RColorBrewer::brewer.pal(8, "Accent")
  )

# UI ----
# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Make Euler or Upset plot from multiple bed files"),
  HTML(
    "
<a href='https://zenodo.org/badge/latestdoi/375552189'><img src='https://zenodo.org/badge/375552189.svg' alt='DOI'></a>
<a href='https://github.com/sansamcl/BedUpsetShinyApp'>
<img alt='GitHub release (latest by date)' src='https://img.shields.io/github/v/release/sansamcl/BedUpsetShinyApp?display_name=tag&label=GitHub%20Release'>
</a>
"
  ),

# Sidebar with a slider input for number of bins
sidebarLayout(sidebarPanel(wellPanel(
  fileInput(
    inputId = "bedFiles",
    label = "BedFiles",
    multiple = T
  )
)),
mainPanel(wellPanel(
  uiOutput("sampleChoices")
))),
sidebarLayout(
  sidebarPanel(
    radioButtons(
      "plotType",
      "Plot Type",
      choices = c("Euler", "Upset"),
      selected = "Upset"
    ),
    sliderInput(
      inputId = "minFractOverlap",
      label = "Minimum Fraction Overlap",
      value = 0,
      min = 0,
      max = 1,
      step = 0.05
    ),
    numericInput(
      inputId = "fontScale",
      label = "Font Scale",
      value = 1
    ),
    numericInput(
      inputId = "plotHeight",
      label = "Plot Height",
      value = 400
    ),
    wellPanel(uiOutput("sampleLabelInputPanel")),
    wellPanel(
      uiOutput("sampleColorInputPanel"),
      sliderInput(
        "transparency",
        "Adjust color transparency",
        min = 0.1,
        max = 0.9,
        step = 0.1,
        value = 0.7
      )
    )
  ),
  mainPanel(
    wellPanel(plotOutput("upsetPlot"),
              uiOutput("downloadPlotButton")),
    wellPanel(verbatimTextOutput("upsetPLotCode"))
  )
)
)

# Server ----
# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize = 100 * 1024 ^ 2)
  
  # Reactives ----
  sampleLabels <- metaReactive2({
    req(input$bedFileChoices)
    metaExpr({data.frame(
      "files" = ..(input$bedFileChoices),
      "labels" = purrr::map_chr(..(input$bedFileChoices), ~ ..(input[[.x]]) %||% "")
    )})
  })
  
  sampleColors <-
    metaReactive2({
      req(input$bedFileChoices)
      metaExpr({
        data.frame(
          "files" = ..(input$bedFileChoices),
          "colors" = purrr::map_chr(..(input$bedFileChoices), ~ input[[paste(.x, "color", sep = "_")]] %||% "")
        )
      })
    })

  datapaths <- 
    metaReactive2({
      req(input$bedFiles)
      req(input$bedFileChoices)
      metaExpr({
        ..(input$bedFiles$datapath[match(input$bedFileChoices, input$bedFiles$name)])})
  })
  
  Beds.Df <- metaReactive2({
    req(datapaths())
    metaExpr({lapply(..(datapaths()), read.table)})
  })
    
  Beds.gr <-
    metaReactive({
      lapply(
        ..(Beds.Df()),
        GenomicRanges::makeGRangesFromDataFrame,
        seqnames.field = "V1",
        start.field = "V2",
        end.field = "V3"
      )}
    ) 
  
  AllBeds.gr <- 
    metaReactive({GenomicRanges::reduce(do.call("c", ..(Beds.gr())))}) 
  
  nmes <-
    metaReactive({
      ..(sampleLabels())$labels[match(..(sampleLabels())$file, basename(..(input$bedFileChoices)))] %>% gsub("\\.[^.]*$", "", .)})
    
  List <-
    metaReactive({
        lapply(..(Beds.gr()),
               function(i) {
                 all.gr <- ..(AllBeds.gr())
                 minOverlap = ..(input$minFractOverlap)
                 hits <- GenomicRanges::findOverlaps(all.gr, i)
                 overlaps <-
                   GenomicRanges::pintersect(all.gr[queryHits(hits)],
                                             i[subjectHits(hits)])
                 percentOverlap <-
                   GenomicRanges::width(overlaps) /
                   GenomicRanges::width(all.gr[queryHits(hits)])
                 sigHits <- hits[percentOverlap > minOverlap]
                 GRovrlps <- all.gr[queryHits(sigHits)]
                 paste(
                   GenomicRanges::seqnames(GRovrlps),
                   GenomicRanges::start(GRovrlps),
                   GenomicRanges::end(GRovrlps),
                   sep = "_"
                 )
                 }
               )
        })
  
  lst <- metaReactive({
    set_names(..(List()), 
              ..(nmes()))})
  
  upsetList <- metaReactive({UpSetR::fromList(..(lst()))})
  
  Plot <-
    metaReactive2({
      req(input$bedFiles$name)
      metaExpr({
        if (..(input$plotType) == "Upset") {
          UpSetR::upset(
            ..(upsetList()),
            nsets = length(..(lst())),
            order.by = "freq",
            text.scale = ..(input$fontScale)
          )
        } else {
          plot(
            eulerr::euler(..(upsetList())),
            quantities = T,
            fills = ..(sampleColors())$colors,
            alpha = ..(input$transparency)
          )
        }})}
    )
  
  # Outputs ----
  output$sampleChoices <- renderUI({
    req(input$bedFiles$name)
    #checkboxGroupInput
    shinyWidgets::multiInput(
      "bedFileChoices",
      "Choose Bed Files",
      input$bedFiles$name,
      selected = input$bedFiles$name
    )
  })
  
  output$upsetPlot <-
    metaRender(renderPlot, {
     ..(Plot())}
    )
  
  output$upsetPLotCode <- 
    renderPrint({
      expandChain(output$upsetPlot())
    })
  
  output$sampleLabelInputPanel <- renderUI({
    purrr::map(input$bedFileChoices, ~ textInput(.x, .x, .x))
  })
  
  output$sampleColorInputPanel <- renderUI({
    purrr::map(
      input$bedFileChoices,
      ~ shinyWidgets::spectrumInput(paste(.x, "color", sep = "_"),
                                    .x,
                                    choices = colorList)
    )
  })
  
  output$downloadPlotButton <- renderUI({ 
    req(Plot())
    downloadButton("downloadPlot",
                   "Download Plot")
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      "plot.pdf"
    },
    content = function(file) {
      if (input$plotType == "Euler") {
        ggplot2::ggsave(file, plot = Plot(), device = "pdf")
      } else {
        ggplot2::ggsave(file,
                        plot = ggplotify::as.ggplot(Plot()),
                        device = "pdf")
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
