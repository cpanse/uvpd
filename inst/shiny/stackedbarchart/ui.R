#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    # Application title
    titlePanel(paste("On Analysing Ultraviolet Photodissociation Fragment Spectra", "- version", packageVersion('uvpd'))),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            tagList(
                htmlOutput('selectRData'),
                htmlOutput('selectCompound'),
                checkboxInput("removePC", "remove precursor items", value = FALSE, width = NULL),
                checkboxInput("negativeIonType", "negative/positive ion type", value = FALSE, width = NULL),
                selectInput("ppmerror", "ppm error cut-off", c(1, 5, 10, 15, 20, 50, 100), multiple = FALSE, selected = 10),
                selectInput("epserror", "absolute error cut-off", c(0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.5), multiple = FALSE, selected = 0.002),
                # Button
                downloadButton("downloadScores", "Download Scores"),
                downloadButton("downloadFreq", "Download Frequency"),
                 htmlOutput('selectCluster')
            )),
        
        # Show a plot of the generated distribution
        mainPanel(
            tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"),
            tabsetPanel(
            
            tabPanel("stacked", list(
                column(width = 10,
                       plotOutput("scorePlot"),
                       htmlOutput("stackedBarChartText"),
                       plotOutput("stackedBarChart"),
                       plotOutput("stackedBarChartIonType"),
                       plotOutput("top3"),
                       plotOutput("bwplot")
                ))),
            tabPanel("stacked group", list(
                column(width = 10,
                       plotOutput("stackedBarChartGroup")
                       
                ))),
            tabPanel("error", list(
                column(width = 10,       
                       tableOutput('tableFreq'),
                       plotOutput("distPlot")
                ))),
            tabPanel("ms2-table", list(
                column(width = 10,  
                       DT::dataTableOutput('tableMS2')
                ))),
            tabPanel("ms2", list(
                column(width = 10,             
                       plotOutput("xyplot")
                ))),
            # plotOutput("barchart")
            tabPanel("matchen in-silico cluster table", list(
                column(width = 10,   
                       tableOutput('ThermoUVPD')
                ))),
            tabPanel("data", list(
                column(width = 10,   
                       tableOutput('tableFilteredData')
                ))),
            tabPanel("score", list(
                column(width = 10,   
                       DT::DTOutput('tableScore')
                ))),
            tabPanel("freqTable", list(
                column(width = 10,   
                       DT::DTOutput('TableFreqAll')
                ))),
            tabPanel("scorePlot", list(
                column(width = 10,   
                       plotOutput("score1Plot"),
                       plotOutput("score2Plot"),
                       plotOutput("score3Plot")
                ))),
            tabPanel("in-silico", list(
                column(width = 10,   
                       tableOutput('tableinSilicoFragmentIon')
                ))),
            tabPanel("summary", list(
                column(width = 10,   
                       tableOutput('summary')
                ))),
            tabPanel("Session Info", verbatimTextOutput("sessionInfo"))
            
        ) 
        )
    )
))
