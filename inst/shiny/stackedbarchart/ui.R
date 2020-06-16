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
    titlePanel("Ultra HRMS in combination with UVPD fragments"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            tagList(
                htmlOutput('selectRData'),
            htmlOutput('selectCompound'),
            selectInput("ppmerror", "ppm error cut-off",c(1, 5, 10, 15, 20, 50, 100), multiple = FALSE, selected = 10)
        )),

        # Show a plot of the generated distribution
        mainPanel(
           
                column(width = 10,
            plotOutput("stackedBarChart"),
            plotOutput("bwplot"),
            hr(),
            tableOutput('table'),
            plotOutput("distPlot"),
            plotOutput("top3"),
            plotOutput("xyplot")
           # plotOutput("barchart")
           
            ) 
        )
    )
))
