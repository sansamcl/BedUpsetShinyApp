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

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Make Euler or Upset plot from multiple bed files"),
  HTML("
<img alt='GitHub release (latest by date)' src='https://img.shields.io/github/v/release/sansamcl/BedUpsetShinyApp?display_name=tag&label=GitHub%20Release'>
  "),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      wellPanel(fileInput(inputId="bedFiles",
                          label="BedFiles",
                          multiple=T))
    ),
    mainPanel(
      wellPanel(uiOutput("sampleChoices"))
    )
  ),
  sidebarLayout(
      sidebarPanel(
        radioButtons("plotType",
                     "Plot Type",
                     choices=c("Euler","Upset"),
                     selected="Upset"),
        numericInput(inputId="fontScale",
                     label="Font Scale",
                     value=1),
        numericInput(inputId="plotHeight",
                     label="Plot Height",
                     value=400),
        wellPanel(uiOutput("sampleLabelInputPanel")),
        wellPanel(uiOutput("sampleColorInputPanel"),
                  sliderInput("transparency",
                              "Adjust color transparency",
                              min=0.1,
                              max=0.9,
                              step=0.1,
                              value=0.7))
      ),
      mainPanel(
        wellPanel(plotOutput("upsetPlot"),
                  uiOutput("downloadPlotButton"))
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  makeUpsetFromBeds <- function(BedFilenamesFull,BedFilenames,plotChoice,sampleLabelsDF,colrs){
    Beds.df <- lapply(BedFilenamesFull,read.table)
    Beds.gr <- lapply(Beds.df,
                      GenomicRanges::makeGRangesFromDataFrame,
                      seqnames.field="V1",
                      start.field="V2",
                      end.field="V3")
    AllBeds.gr <- GenomicRanges::reduce(do.call("c",Beds.gr))
    GetOverlapsWithAll <- function(bed.gr,all.gr){
      GRovrlps <- all.gr[GenomicRanges::countOverlaps(all.gr,
                                                      bed.gr)>0]
      paste(GenomicRanges::seqnames(GRovrlps),
            GenomicRanges::start(GRovrlps),
            GenomicRanges::end(GRovrlps),
            sep="_")
    }
    lst <- lapply(Beds.gr,GetOverlapsWithAll,AllBeds.gr)
    #nmes <- gsub(".bed","",basename(BedFilenames))
    nmes <- sampleLabelsDF$labels[match(sampleLabelsDF$files,basename(BedFilenames))]
    nmes <- gsub("\\.[^.]*$","",nmes)
    names(lst) <- nmes
    #names(lst) <- sampleLabelsDF$labels
    upsetList <- UpSetR::fromList(lst)
    
    
    eulerPlot <- plot(eulerr::euler(upsetList),quantities=T,fills=colrs$colors,alpha=input$transparency)
    
    
    
    plt <- UpSetR::upset(upsetList,
                         nsets=length(lst),
                         order.by = "freq",
                         text.scale=input$fontScale)
    
    if(plotChoice=="Upset"){
      plt
    }else{
      eulerPlot
    }
    
  }
  
  output$sampleChoices <- renderUI({
    req(input$bedFiles$name)
    #checkboxGroupInput
    shinyWidgets::multiInput("bedFileChoices", 
                             "Choose Bed Files", 
                             input$bedFiles$name,
                             selected=input$bedFiles$name)
  })
  
  datapaths <- reactive({
    req(input$bedFiles)
    req(input$bedFileChoices)
    input$bedFiles$datapath[match(input$bedFileChoices,input$bedFiles$name)]
  })
  
  eulerPlot <- reactive({
    req(datapaths())
    req(input$bedFileChoices)
    makeUpsetFromBeds(datapaths(),
                      input$bedFileChoices,
                      "Euler",
                      sampleLabels(),
                      sampleColors())
  })
  
  upsetPLot <- reactive({
    req(datapaths())
    req(input$bedFileChoices)
    makeUpsetFromBeds(datapaths(),
                      input$bedFileChoices,
                      "Upset",
                      sampleLabels(),
                      sampleColors())
  })
  
  output$upsetPlot <- renderPlot({
    if(input$plotType=="Euler"){
      eulerPlot()
    }else{
      upsetPLot()
    }
    
  },height=reactive(input$plotHeight))
  
  output$sampleLabelInputPanel <- renderUI({
    purrr::map(input$bedFileChoices, ~ textInput(.x, .x, .x))
  })
  
  sampleLabels <- reactive({
    req(input$bedFileChoices)
    data.frame("files"=input$bedFileChoices,
               "labels"=purrr::map_chr(
                 input$bedFileChoices, ~ input[[.x]] %||% "")
               )
  })
  ###
  # output$sampleColorInputPanel <- renderUI({
  #   purrr::map(input$bedFileChoices, ~ colourpicker::colourInput(paste(.x,"color",sep="_"), 
  #                                                                .x, 
  #                                                                .x,
  #                                                                allowTransparent = TRUE,
  #                                                                palette="limited",
  #                                                                allowedCols=RColorBrewer::brewer.pal(8, "Dark2")
  #                                                                ))
  # })
  
  colorList <- reactive({
    list(
      c("#1B9E77","#D95F02","#bdbdbd"),
      RColorBrewer::brewer.pal(9, "Set1"),
      RColorBrewer::brewer.pal(8, "Set2"),
      RColorBrewer::brewer.pal(12, "Set3"),
      RColorBrewer::brewer.pal(9, "Pastel1"),
      RColorBrewer::brewer.pal(8, "Pastel2"),
      RColorBrewer::brewer.pal(12, "Paired"),
      RColorBrewer::brewer.pal(8, "Dark2"),
      RColorBrewer::brewer.pal(8, "Accent")
    )
    # lst2 <- lapply(lst,GISTools::add.alpha,input$transparency)
    # tp <- lapply(lst2,function(vr){sapply(vr,stringr::str_extract,"..$")})[[1]][1]
    # hex <- lapply(lst,function(vr){sapply(vr,stringr::str_extract,"[^#]*$")})
    # lapply(hex,function(x){paste("#",tp,x,sep="")})
    
  })
  
  output$sampleColorInputPanel <- renderUI({
    purrr::map(input$bedFileChoices, ~ shinyWidgets::spectrumInput(paste(.x,"color",sep="_"), 
                                                                 .x,
                                                                 choices=colorList()
    ))
  })
  
  sampleColors <- reactive({
    req(input$bedFileChoices)
    data.frame("files"=input$bedFileChoices,
               "colors"=purrr::map_chr(
                 input$bedFileChoices, ~ input[[paste(.x,"color",sep="_")]] %||% "")
    )
  })
  ###
  output$downloadPlotButton <- renderUI({
    req(upsetPLot())
    req(eulerPlot())
    downloadButton("downloadPlot",
                   "Download Plot")
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { "plot.pdf" },
    content = function(file) {
      ggplot2::ggsave(file, plot = eulerPlot(), device = "pdf")
      if(input$plotType=="Euler"){
        
      }else
        ggplot2::ggsave(file, plot = upsetPLot(), device = "pdf")
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
