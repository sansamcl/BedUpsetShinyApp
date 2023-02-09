#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(purrr)
library(shinymeta)
require(GenomicRanges)
require(ggplot2)
require(UpSetR)
require(eulerr)
require(ggplotify)

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

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2) 
  Beds.Df <- metaReactive2({
    req(datapaths())
    lapply(datapaths(), read.table)
  })
  
  makePlotFromBeds <-
    function(BedFilenames,
             plotChoice,
             sampleLabelsDF,
             colrs,
             alfa,
             fontScale,
             minFracOverlap) {
      Beds.df <- Beds.Df()
      Beds.gr <- lapply(
        Beds.df,
        GenomicRanges::makeGRangesFromDataFrame,
        seqnames.field = "V1",
        start.field = "V2",
        end.field = "V3"
      )
      AllBeds.gr <- GenomicRanges::reduce(do.call("c", Beds.gr))
      GetOverlapsWithAll <- function(bed.gr, all.gr, minOverlap) {
        hits <- GenomicRanges::findOverlaps(all.gr, bed.gr)
        overlaps <- GenomicRanges::pintersect(all.gr[queryHits(hits)],
                                              bed.gr[subjectHits(hits)])
        percentOverlap <-
          GenomicRanges::width(overlaps) / GenomicRanges::width(all.gr[queryHits(hits)])
        sigHits <- hits[percentOverlap > minOverlap]
        GRovrlps <- all.gr[queryHits(sigHits)]
        paste(
          GenomicRanges::seqnames(GRovrlps),
          GenomicRanges::start(GRovrlps),
          GenomicRanges::end(GRovrlps),
          sep = "_"
        )
      }
      lst <-
        lapply(Beds.gr, GetOverlapsWithAll, all.gr=AllBeds.gr, minOverlap=minFracOverlap)
      #nmes <- gsub(".bed","",basename(BedFilenames))
      nmes <-
        sampleLabelsDF$labels[match(sampleLabelsDF$files, basename(BedFilenames))]
      nmes <- gsub("\\.[^.]*$", "", nmes)
      names(lst) <- nmes
      if (max(sapply(lst,length)) > 0 ){
        upsetList <- UpSetR::fromList(lst)
        eulerPlot <- plot(eulerr::euler(upsetList), quantities = T, fills = colrs$colors, alpha = alfa)
        plt <- UpSetR::upset(upsetList, nsets = length(lst), order.by = "freq", text.scale = fontScale)
        if (plotChoice == "Upset") {
          plt
        } else {
          eulerPlot
        }
      } else {
        par(mar = c(0, 0, 0, 0))
        plot(x = 0:1,
             y = 0:1,
             ann = F,
             bty = "n",
             type = "n",
             xaxt = "n",
             yaxt = "n")
        text(x = 0.5,
             y = 0.5,
             "No Overlaps Detected", 
             cex = 1.8)
        }
    }
    
  
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
  
  datapaths <- metaReactive2({
    req(input$bedFiles)
    req(input$bedFileChoices)
    input$bedFiles$datapath[match(input$bedFileChoices, input$bedFiles$name)]
  })
  
  eulerPlot <- metaReactive2({
    req(Beds.Df())
    req(datapaths())
    req(input$bedFileChoices)
    req(input$minFractOverlap)
    metaExpr(makePlotFromBeds(
      ..(input$bedFileChoices),
      "Euler",
      ..(sampleLabels()),
      ..(sampleColors()),
      ..(input$transparency),
      ..(input$fontScale),
      ..(input$minFractOverlap)
    ))
  })
  
  upsetPLot <- metaReactive2({
    req(Beds.Df())
    req(datapaths())
    req(input$bedFileChoices)
    metaExpr(makePlotFromBeds(
      ..(input$bedFileChoices),
      "Upset",
      ..(sampleLabels()),
      ..(sampleColors()),
      ..(input$transparency),
      ..(input$fontScale),
      ..(input$minFractOverlap)
    ))
  })
  
  output$upsetPLotCode <- renderPrint({
    if (input$plotType == "Euler") {
      expandChain(quote(makePlotFromBeds <-
                          function(BedFilenames,
                                   plotChoice,
                                   sampleLabelsDF,
                                   colrs,
                                   alfa,
                                   fontScale,
                                   minFracOverlap) {
                            Beds.df <- lapply(BedFilenames, read.table)
                            Beds.gr <- lapply(
                              Beds.df,
                              GenomicRanges::makeGRangesFromDataFrame,
                              seqnames.field = "V1",
                              start.field = "V2",
                              end.field = "V3"
                            )
                            AllBeds.gr <- GenomicRanges::reduce(do.call("c", Beds.gr))
                            GetOverlapsWithAll <-
                              function(bed.gr, all.gr, minOverlap) {
                                hits <- GenomicRanges::findOverlaps(all.gr, bed.gr)
                                overlaps <- GenomicRanges::pintersect(all.gr[queryHits(hits)],
                                                                      bed.gr[subjectHits(hits)])
                                percentOverlap <-
                                  GenomicRanges::width(overlaps) / GenomicRanges::width(all.gr[queryHits(hits)])
                                sigHits <- hits[percentOverlap > minOverlap]
                                GRovrlps <- all.gr[queryHits(sigHits)]
                                paste(
                                  GenomicRanges::seqnames(GRovrlps),
                                  GenomicRanges::start(GRovrlps),
                                  GenomicRanges::end(GRovrlps),
                                  sep = "_"
                                )
                              }
                            lst <-
                              lapply(Beds.gr, GetOverlapsWithAll, all.gr=AllBeds.gr, minOverlap=minFracOverlap)
                            nmes <-
                              sampleLabelsDF$labels[match(sampleLabelsDF$files, basename(BedFilenames))]
                            nmes <- gsub("\\.[^.]*$", "", nmes)
                            names(lst) <- nmes
                            upsetList <- UpSetR::fromList(lst)
                            
                            
                            eulerPlot <-
                              plot(
                                eulerr::euler(upsetList),
                                quantities = T,
                                fills = colrs$colors,
                                alpha = alfa
                              )
                            
                            
                            
                            plt <- UpSetR::upset(
                              upsetList,
                              nsets = length(lst),
                              order.by = "freq",
                              text.scale = fontScale
                            )
                            
                            if (plotChoice == "Upset") {
                              plt
                            } else{
                              eulerPlot
                            }
                            
                          }),
                  eulerPlot())
    } else{
      expandChain(quote(makePlotFromBeds <-
                          function(BedFilenames,
                                   plotChoice,
                                   sampleLabelsDF,
                                   colrs,
                                   alfa,
                                   fontScale,
                                   minFracOverlap) {
                            Beds.df <- lapply(BedFilenames, read.table)
                            Beds.gr <- lapply(
                              Beds.df,
                              GenomicRanges::makeGRangesFromDataFrame,
                              seqnames.field = "V1",
                              start.field = "V2",
                              end.field = "V3"
                            )
                            AllBeds.gr <-
                              GenomicRanges::reduce(do.call("c", Beds.gr))
                            GetOverlapsWithAll <- function(bed.gr, all.gr,minOverlap) {
                              hits <- GenomicRanges::findOverlaps(all.gr, bed.gr)
                              overlaps <- GenomicRanges::pintersect(all.gr[queryHits(hits)],
                                                                    bed.gr[subjectHits(hits)])
                              percentOverlap <-
                                GenomicRanges::width(overlaps) / GenomicRanges::width(all.gr[queryHits(hits)])
                              sigHits <- hits[percentOverlap > minOverlap]
                              GRovrlps <- all.gr[queryHits(sigHits)]
                              paste(
                                GenomicRanges::seqnames(GRovrlps),
                                GenomicRanges::start(GRovrlps),
                                GenomicRanges::end(GRovrlps),
                                sep = "_"
                              )
                            }
                            lst <-
                              lapply(Beds.gr,
                                     GetOverlapsWithAll,
                                     all.gr=AllBeds.gr, minOverlap=minFracOverlap)
                            nmes <-
                              sampleLabelsDF$labels[match(sampleLabelsDF$files, basename(BedFilenames))]
                            nmes <- gsub("\\.[^.]*$", "", nmes)
                            names(lst) <- nmes
                            upsetList <- UpSetR::fromList(lst)
                            
                            
                            eulerPlot <-
                              plot(
                                eulerr::euler(upsetList),
                                quantities = T,
                                fills = colrs$colors,
                                alpha = alfa
                              )
                            
                            
                            
                            plt <- UpSetR::upset(
                              upsetList,
                              nsets = length(lst),
                              order.by = "freq",
                              text.scale = fontScale
                            )
                            
                            if (plotChoice == "Upset") {
                              plt
                            } else{
                              eulerPlot
                            }
                            
                          }),
                  upsetPLot())
    }
  })
  
  
  output$upsetPlot <- renderPlot({
    if (input$plotType == "Euler") {
      eulerPlot()
    } else{
      upsetPLot()
    }
    
  }, height = metaReactive2({
    input$plotHeight
  }))
  
  output$sampleLabelInputPanel <- renderUI({
    purrr::map(input$bedFileChoices, ~ textInput(.x, .x, .x))
  })
  
  sampleLabels <- metaReactive2({
    req(input$bedFileChoices)
    data.frame(
      "files" = input$bedFileChoices,
      "labels" = purrr::map_chr(input$bedFileChoices, ~ input[[.x]] %||% "")
    )
  })
  
  colorList <- metaReactive2({
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
  })
  
  output$sampleColorInputPanel <- renderUI({
    purrr::map(
      input$bedFileChoices,
      ~ shinyWidgets::spectrumInput(paste(.x, "color", sep = "_"),
                                    .x,
                                    choices =
                                      colorList())
    )
  })
  
  sampleColors <- metaReactive2({
    req(input$bedFileChoices)
    data.frame(
      "files" = input$bedFileChoices,
      "colors" = purrr::map_chr(input$bedFileChoices, ~ input[[paste(.x, "color", sep =
                                                                       "_")]] %||% "")
    )
  })

  
  output$downloadPlotButton <- renderUI({
    req(upsetPLot())
    req(eulerPlot())
    downloadButton("downloadPlot",
                   "Download Plot")
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      "plot.pdf"
    },
    content = function(file) {
      if (input$plotType == "Euler") {
        ggplot2::ggsave(file, plot = eulerPlot(), device = "pdf")
      } else {
        ggplot2::ggsave(file, plot = ggplotify::as.ggplot(upsetPLot()),device="pdf")
        }
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)
