#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(uvpd)
helpfile <- file.path(path.package(package = 'uvpd'), 'shiny', 'stackedbarchart', 'help.md')
# Define UI for application that draws a histogram
shinyUI(fluidPage(
    
    # Application title
    titlePanel(paste("On Analysing Ultraviolet Photodissociation Fragment Spectra", "- version", packageVersion('uvpd'))),
    
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            tagList(
                h3("input:"),
                htmlOutput('selectRData'),
                hr(),
                h3("filter settings:"),
                htmlOutput('selectCompound'),
                
                checkboxInput("removePC", "remove precursor items", value = FALSE, width = NULL),
                checkboxInput("negativeIonType", "negative/positive ion type", value = FALSE, width = NULL),
                hr(),
                selectInput("ppmerror", "ppm error cut-off", c(1, 5, 10, 15, 20, 50, 100), multiple = FALSE, selected = 10),
                selectInput("epserror", "absolute error cut-off", c(0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.05, 0.1, 0.5), multiple = FALSE, selected = 0.002),
                
                htmlOutput('selectCluster')
            )),
        
        # Show a plot of the generated distribution
        mainPanel(
            tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"),
            tabsetPanel(
               
                tabPanel("stacked fragments", list(
                    column(width = 10,
                           tableOutput("stackedBarChartText"),
                           plotOutput("scorePlot"),
                           helpText("The stacked bar charts draw the logarithmically transformed fragment ion intensities for the different ion fragmentations and types, respectively. If marked in the checkbox (left), the unfragmented precursor peaks are removed."),
                           
                           plotOutput("stackedBarChart"),
                           plotOutput("stackedBarChartIonType"),
                           helpText("Each panel plots the top three most abundant total ion count (TIC) versus the master intensity of one raw file. The color linking indicates for the raw file grouping."),
                           plotOutput("top3"),
                           helpText("The boxplots draw the absolute error (in Dalton) distribution."),
                           plotOutput("bwplot")
                    ))),
                tabPanel("summary", list(
                    column(width = 10,      
                           helpText("Some statistics of the overall data and the applied filter setting is shown below:"),
                           tableOutput('summary'),
                           helpText("In the frequency table, the frequency of each compound and fragmentation type is determined under consideration of the ionization mode (positive or negative) and the filter settings (ppm error and absolute error)."),
                           tableOutput('tableFreq'),
                           helpText("The histograms display the distribution of the ppm and the absolute error over the entire data set and the selected compound. The red curve depicts the maximum-likelihood fitting assuming an underlying normal distribution."),
                           plotOutput("distPlot")
                    ))),
                tabPanel("ms2", list(
                    column(width = 10,     
                           helpText("The table lists the fragment ion for each fragmentation type."),
                           hr(),
                           DT::dataTableOutput('tableMS2'),
                           helpText("Each panel displays the fragment ion spectrum for each fragmentation type. The color indicates the file origin."),
                           plotOutput("xyplot")
                    ))),
                tabPanel("data", list(
                    column(width = 10,   
                           helpText("lists the entire data set applying the filter setting."),
                           DT::dataTableOutput('tableFilteredData')
                    ))),
                tabPanel("scores", list(
                    column(width = 10,   
                           helpText("This page provides access to the computed scores. The scores are defined as follow:"),
                           HTML("<p><dfn>score1</dfn> experimental fragments matched &frasl; theoretically possible.</p>"),
                           HTML("<p><dfn>score2</dfn> experimental fragments matched &frasl;  all fragments in spectrum.</p>"),
                           HTML("<p><dfn>score3</dfn> experimental matched fragment intensities &frasl; master.intensity.</p>"),
                           downloadButton("downloadScores", "Download Scores"),
                           hr(),
                           DT::DTOutput('tableScore'),
                           hr(),
                           plotOutput("score1Plot"),
                           plotOutput("score2Plot"),
                           plotOutput("score3Plot")
                    ))),
                tabPanel("frequencies", list(
                    column(width = 10,   
                           helpText("In the frequency table, the frequency of each compound and fragmentation type is determined under consideration of the ionization mode (positive or negative) and the filter settings (ppm error and absolute error)."),
                           downloadButton("downloadFreq", "Download the frequency table"),
                           HTML("<hr>"),
                           DT::DTOutput('TableFreqAll')
                    ))),
#                tabPanel("cluster", list(
#                    column(width = 10,   
#                           helpText("The stacked bar charts draw the logarithmically transformed fragment ion intensities for the different ion fragmentations of the same cluster assignments.
# The table shows the cluster assignments of all compounds. The clustering ware performed on a preliminary investigation. The to be displayed cluster ID can be set in the left panel."),
#                           plotOutput("stackedBarChartGroup"),
#                           tableOutput('ThermoUVPD')
#                    ))),
                tabPanel("predicted ion", list(
                    column(width = 10,   
                           helpText("By using the method `metfRag::frag.generateFragments` derived fragment ions are listed below. 
Those ions represent the entire search space for the peak assignment."),
                           DT::dataTableOutput('tableinSilicoFragmentIon')
                    ))),
                #tabPanel("R session information", verbatimTextOutput("sessionInfo"))
                tabPanel("help", list(
                    column(width = 10,
                           helpText(paste("Help uvpd -", "version", packageVersion('uvpd'))),
                           includeMarkdown(helpfile),
                           hr(),
                           HTML("<h3>Howto cite?</h3>"),
                           verbatimTextOutput("citation"),
                           HTML("<h3>Session Information</h3>"),
                           helpText("below is the output of the method call sessionInfo()"),
                           verbatimTextOutput("sessionInfo")
                           
                    )))
            ) 
        )
    )
))
