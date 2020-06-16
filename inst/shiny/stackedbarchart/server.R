#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(readr)
library(shiny)
library(ggplot2)
library(DT)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$selectRData <- renderUI({
    selectInput("rdata", "rdata", c('uvpd.20200522.RData', 'uvpd.20200612.RData'),
                 multiple = FALSE, selected = 'uvpd.20200612.RData')
  })
  

  getData <- reactive({
    e <- new.env()
    
    fn <- file.path(system.file(package = 'uvpd'), "extdata", input$rdata)
   
    load(fn, envir = e)
    X<- e$X.top3.master.intensity
    Y <- Y <- do.call('rbind', do.call('rbind',  e$X.top3.master.intensity.MS2))
    
    X$m <- paste(X$file, X$scan)
    Y$m <- paste(Y$file, Y$scan)
    XY <- base::merge(X, Y, by="m")
    
    XY$fragmode <- gsub("uvpd50.00", "uvpd050.00", XY$fragmode)
    XY$fragmode <- gsub("uvpd25.00", "uvpd025.00", XY$fragmode)
    XY
  })
  
  
  getThermoUVPD_feb2019 <- reactive({
    
    fn <- file.path(system.file(package = 'uvpd'), "extdata/ThermoUVPD_feb2019.csv")
    ThermoUVPD_feb2019 <- read_csv(fn, na = c("", "#N/A"))
    
    ThermoUVPD_feb2019[ThermoUVPD_feb2019$Compound %in% getData()$Compound,
                       c("Compound", "Bruto formula", "Cluster number", "Group")]
  })
  
  
  
  getFilteredData <- reactive({
    DF <- getData()
    
    DF[DF$compound %in% input$compound & DF$ppmerror < as.numeric(input$ppmerror), ]
  })
  
  getAggregatedData <- reactive({
    aggregate(intensity ~ mZ * file.y * fragmode * compound * formula * Group * mode,
              data=getFilteredData(), FUN=sum)
  })
  
  getLevels <- reactive({
    DF <- getData()
    sort(unique(DF$fragmode))
  })
  
  getCompound <-  reactive({
    message("getCompound")
    DF <- getData()
    rv <- sort(unique(DF$compound))
    message(rv)
    rv
  })
  
  output$selectCompound <- renderUI({
    selectInput("compound", "compound", getCompound(), multiple = FALSE, selected = "Triadimenol")
  })
  
  output$distPlot <- renderPlot({
    par(mfrow=c(1, 2))
    hist(getFilteredData()$ppmerror, main=input$ppmerror)
    hist(getFilteredData()$eps)
    return
  })
  
  
  getTable <- reactive({
    DF <- getAggregatedData()
    print(names(DF))
    rv <- table(unique(subset(DF, select=c('fragmode','formula')))[,1])
    
    rv <- (as.data.frame(rv))
  })
  
  output$ThermoUVPD <- renderTable({
    getThermoUVPD_feb2019()  
  })
  
  output$table <- renderTable({
    
    getTable()
  })
  
  output$stackedBarChart <- renderPlot({
    
    gp <- ggplot(data = getAggregatedData(),
                 aes(x = factor(fragmode, levels = getLevels()),
                     y = log(intensity, 10),
                     fill=reorder(formula, mZ))) +
      geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
      scale_x_discrete(drop=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      facet_wrap(~ compound * Group * mode, scales="free", drop=FALSE)
    
    # gp2 <- ggplot(data=unique(subset(S, select=c('fragmode','formula'))), aes(x=factor(fragmode, levels = sort(unique(DF$fragmode))), fill=(formula))) + ggplot2::geom_bar()
    gp
  },width=1000,height = 400)
  
  output$stackedBarChartGroup <- renderPlot({
    
    gp <- ggplot(data = getAggregatedData(),
                 aes(x = factor(fragmode, levels = getLevels()),
                     y = log(intensity, 10),
                     fill=reorder(formula, mZ))) +
      geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
      scale_x_discrete(drop=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      facet_wrap(~ compound * Group * mode, scales="free", drop=FALSE)
    
    # gp2 <- ggplot(data=unique(subset(S, select=c('fragmode','formula'))), aes(x=factor(fragmode, levels = sort(unique(DF$fragmode))), fill=(formula))) + ggplot2::geom_bar()
    gp
  },width=1000,height = 400)
  
  
  output$xyplot <- renderPlot({
    lattice::xyplot(intensity ~ mZ | fragmode, 
                    group=paste(file.y, scan.y),
                    scales = list(y="free"),
                    auto.key =list(space = "right"),
                    data=getFilteredData(),
                    type='h', layout = c(1,13))
  }, width=1000, height=1000)
  
  output$bwplot <- renderPlot({
    
    lattice::bwplot(eps ~  fragmode, 
                    data= getFilteredData(),
                    scales = list(x = list(rot = 45)))
  }, width=1000, height=200)
  
  output$barchart <- renderPlot({
    lattice::barchart(log(intensity) ~ fragmode | compound,
                      group=mode ,
                      data=aggregate(intensity ~ mode * fragmode * compound, data=getFilteredData(), FUN=sum),
                      scales = list(x = list(rot = 45)),
                      layout=c(6,7))
  }, width=1000, height=1000)
  
  output$top3 <- renderPlot({
    lattice::xyplot(tic ~ master.intensity | fragmode * Compound,
                    #subset = Compound == "Triadimenol",
                    group=fragmode,
                    data=getFilteredData())
  }, width=1000, height = 400)
})
