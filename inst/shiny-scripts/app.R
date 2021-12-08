library(shiny)

# Adapted from MPLNClust by Anjali Silva
# https://github.com/anjalisilva/MPLNClust/

ui <- fluidPage(
    # Create title header
    titlePanel(tags$h1(tags$b("EpiGPlot:"),"Epigenetic factors expression visualization")),
    # Create layout for two panels, sidebar for loading and plotting purposes, main bar for outputs
    sidebarLayout(
        sidebarPanel(
            tags$p("Given a data set of a gene/epigenetic factor, expression values is plotted against quantile over all genes."),
            tags$div(
                style="background-color: #DCDCDC; border-radius:8px; padding: 5px",
                tags$body(
                    p("Import a data set of a sample gene/epigenetic factor from .csv/.rda files"),
                    p("Rows should correspond to the sample's interaction with other samples."),
                    p("Columns should correspond to the information of each samples."),
                    p(".csv files must contain column names in the first row."),
                    p("Only .rda and .csv files are accepted.")
                )
            ),
            tags$div(
                style="background-color: #DCDCDC; border-radius:8px; padding: 5px; margin-top: 8px",
                tags$body(
                    fileInput(inputId = "csvInput",
                              label = "Select a dataset to import:",
                              placeholder = ".csv/.rda files",
                              accept = c(".csv", ".rda")),

                    actionButton(inputId = "button1",
                                 label = "Start Plotting"),
                )
            ),


        ),
        mainPanel(
            # Create tabs on the main panel for settings, input summary, plot, and calculated data
            tabsetPanel(type="tabs",
                        tabPanel("Plot Settings",
                                 tags$div(
                                     style="background-color: #DCDCDC;
                                     border-radius:8px; overflow: hidden;
                                     padding: 5px; margin-top: 8px",
                                     tags$div(
                                         style="width: 48%; float:left; margin-left: 8px;",
                                         tags$h5("Plot features"), hr(),
                                         checkboxInput(inputId = "pDisplayLegend",
                                                       "Display Legend", value=TRUE),
                                         checkboxInput(inputId = "pPlotGreyscale",
                                                       "Greyscale Plot", value=FALSE),
                                         checkboxInput(inputId = "pNormalizeData",
                                                       "Normalize Data", value=FALSE),
                                         numericInput(inputId =  "pArrowSize",
                                                      "Arrow Size", 1),
                                         numericInput(inputId =  "pLineSize",
                                                      "Linear regression line size", 1)
                                     ),
                                     tags$div(
                                         style="width: 48%; float:left; margin-left: 8px;",
                                         tags$h5("Plot's Labels settings"), hr(),
                                         numericInput(inputId =  "pTopLabels",
                                                      "Display number of labels of top expressions values", 0),
                                         numericInput(inputId =  "pLabelsSize",
                                                      "Label size of top expression values", 3),
                                         textAreaInput(inputId="pLabelsWhich", "Specify sample names to be displayed(separated by new line)",value="serous cystadenocarcinoma cell line:HTOA")
                                     )
                                 ),
                                 tags$div(
                                     style="background-color: #DCDCDC;
                                     border-radius:8px; overflow: hidden;
                                     padding: 5px; margin-top: 8px",
                                     tags$div(
                                         style="width: 48%; float:left; margin-left: 8px;",
                                         tags$h5("Specify data subset for plot"), hr(),
                                         textAreaInput(inputId="pSampleClass", "Specify sample classes to displayed",
                                                       width="100%", resize = "none",
                                                       value="cell_line,fractionation,primary_cell,timecourse,tissue"),
                                         textAreaInput(inputId="pSamplePred", "Perform linear regression on sample classes",
                                                       width="100%", resize = "none",
                                                       value="cell_line")
                                     ),
                                     tags$div(style="width: 48%; float:left; margin-left: 8px;",
                                         tags$h5("Specify sample colours for plot"), hr(),
                                         tags$div(style="width: 48%; float:left; margin-left: 4px;",
                                                  textInput("pCol1", "Colour 1","#1b9e77")),
                                         tags$div(style="width: 48%; float:left; margin-left: 4px;",
                                                  textInput("pCol2", "Colour 2","#d95f02")),
                                         tags$div(style="width: 48%; float:left; margin-left: 4px;",
                                                  textInput("pCol3", "Colour 3","#7570b3")),
                                         tags$div(style="width: 48%; float:left; margin-left: 4px;",
                                                  textInput("pCol4", "Colour 4","#e7298a")),
                                         tags$div(style="width: 48%; float:left; margin-left: 4px;",
                                                  textInput("pCol5", "Colour 5","#66a61e"))
                                     ),
                                     # Create previewable colour palette selection by user
                                     tags$div(
                                         style="float:right",
                                         htmlOutput("pCol1Prev",style="width: 18%; float:left; margin-right:2px"),
                                         htmlOutput("pCol2Prev",style="width: 18%; float:left; margin-right:2px"),
                                         htmlOutput("pCol3Prev",style="width: 18%; float:left; margin-right:2px"),
                                         htmlOutput("pCol4Prev",style="width: 18%; float:left; margin-right:2px"),
                                         htmlOutput("pCol5Prev",style="width: 18%; float:left; margin-right:2px"))

                                 )),
                        tabPanel("Input Summary", verbatimTextOutput("textOut")),
                        tabPanel("Output Plot", plotOutput("outputPlot")),
                        tabPanel("Calculated results", htmlOutput("calcOut"))
                        )
        )
    )
)

server <- function(input, output, session) {
    # Initialize a reactive for data loading
    dataReactive <- reactive({
        if(!is.null(input$csvInput)){
            if(grepl("\\.csv|\\.rda", input$csvInput$datapath, ignore.case = TRUE)){
                if(grepl("\\.csv", input$csvInput$datapath, ignore.case = TRUE)){
                    EpiGPlot::parseEpigeneticData(input$csvInput$datapath)
                }
                else{
                    EpiGPlot::loadEpigeneticData(input$csvInput$datapath)
                }
            }
        }
    })

    # Initialize a reactive for plot layout creation
    startLayoutCreate <- reactive({
        layout <- dataReactive()
        if(!is.null(layout)){
            EpiGPlot::layoutEpigeneticEV(dataReactive(),
                                         normalized=as.logical(input$pNormalizeData),
                                         sample.class=strsplit(input$pSampleClass, split=",")[[1]],
                                         class.colour=c(input$pCol1, input$pCol2, input$pCol3, input$pCol4, input$pCol5))
        }
    })

    # Initialize a reactive event for plot generation
    startPlotCreate <- eventReactive(
        eventExpr = input$button1,
        {
            withProgress(
                message = "Generating Visualization Plot",
                value = 0,
                {
                    EpiGPlot::plotEpigeneticEV(
                        startLayoutCreate(),
                        normalize = input$pNormalizeData,
                        colour.legend = input$pDisplayLegend,
                        colour.grey = input$pPlotGreyscale,
                        size.arrow = input$pArrowSize,
                        size.pred = input$pLineSize,
                        labels.which = strsplit(input$pLabelsWhich, split =
                                                    "\n")[[1]],
                        labels.top = input$pTopLabels,
                        labels.size = input$pLabelsSize,
                        sample.class = strsplit(input$pSampleClass, split =
                                                    ",")[[1]],
                        sample.pred = strsplit(input$pSamplePred, split =
                                                   ",")[[1]]
                    )

        })
    })

    # Output information of the raw input data
    output$textOut <- renderPrint({
        if (!is.null(input$csvInput)){
            inputDataframe <- dataReactive()
            cat("Number of rows of input data:", nrow(inputDataframe), "\n")
            cat("Number of columns of input data:", ncol(inputDataframe), "\n")
            cat("Loaded column names:", paste(colnames(inputDataframe), collapse=", "), "\n")
            cat("Unique Sample Classes:", paste(unique(inputDataframe[,1]), collapse=","), "\n")
            cat("summary():\n")
            summary(inputDataframe)
        }
    })

    # Output information of processed input data
    output$calcOut <- renderUI({
        t <- startLayoutCreate()
        t <- cbind(t, Sample=rownames(t))
        renderTable(t)
    })

    # Create plot through initialized reactive
    output$outputPlot <- renderPlot({
        startPlotCreate()
    })


    # Preview user's colour palette by rendering a rectangle svg element
    output$pCol1Prev <- renderUI({
        tag("svg", varArgs = list(width="20px", height="20px",
                                  tags$rect(width="20px", height="20px",
                                            style=paste0("fill: ", input$pCol1, "; "))))
    })
    output$pCol2Prev <- renderUI({
        tag("svg", varArgs = list(width="20px", height="20px",
                                  tags$rect(width="20px", height="20px",
                                            style=paste0("fill: ", input$pCol2, "; "))))
    })
    output$pCol3Prev <- renderUI({
        tag("svg", varArgs = list(width="20px", height="20px",
                                  tags$rect(width="20px", height="20px",
                                            style=paste0("fill: ", input$pCol3, "; "))))
    })
    output$pCol4Prev <- renderUI({
        tag("svg", varArgs = list(width="20px", height="20px",
                                  tags$rect(width="20px", height="20px",
                                            style=paste0("fill: ", input$pCol4, "; "))))
    })
    output$pCol5Prev <- renderUI({
        tag("svg", varArgs = list(width="20px", height="20px",
                                  tags$rect(width="20px", height="20px",
                                            style=paste0("fill: ", input$pCol5, "; "))))
    })
}
shiny::shinyApp(ui = ui, server = server)
