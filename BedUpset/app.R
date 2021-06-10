#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Make an Upset plot from multiple bed files"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId="bedFiles",
                      label="BedFiles",
                      multiple=T),
            uiOutput("sampleChoices"),
            numericInput(inputId="fontScale",
                         label="Font Scale",
                         value=1),
            numericInput(inputId="plotHeight",
                         label="Plot Height",
                         value=400)
        ),
        mainPanel(
            plotOutput("upsetPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    makeUpsetFromBeds <- function(BedFilenamesFull,BedFilenames){
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
        nmes <- gsub(".bed","",basename(BedFilenames))
        names(lst) <- nmes
        upsetList <- UpSetR::fromList(lst)
        plt <- UpSetR::upset(upsetList,
                             nsets=length(lst),
                             order.by = "freq",
                             text.scale=input$fontScale)
        plt
    }
    
    output$sampleChoices <- renderUI({
        #checkboxGroupInput
        shinyWidgets::multiInput("bedFileChoices", 
                           "Choose Bed Files", 
                           input$bedFiles$name,
                           selected=input$bedFiles$name)
    })
    
    datapaths <- reactive({
        input$bedFiles$datapath[match(input$bedFileChoices,input$bedFiles$name)]
        })
    
    output$upsetPlot <- renderPlot({
        req(input$bedFiles)
        makeUpsetFromBeds(datapaths(),
                          input$bedFileChoices)
    },height=reactive(input$plotHeight))

}

# Run the application 
shinyApp(ui = ui, server = server)
