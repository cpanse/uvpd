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
    selectInput("rdata", "rdata", c('uvpd_20200522.RData', 'uvpd_20200612.RData'),
                 multiple = FALSE, selected = 'uvpd_20200612.RData')
  })
  
  output$selectCluster <- renderUI({
    selectInput("clusterid", "clusterid", getClsuterIds(),
                multiple = FALSE, selected = getClsuterIds()[1])
  })

  
  # ease for debugging
  .gd <- function(fn=file.path(system.file(package = 'uvpd'), "extdata", 'uvpd_20200612.RData')){
    e <- new.env()
    
    load(file.path(system.file(package = 'uvpd'), "extdata", "fragments.RData"), envir = e)
    pp <- data.frame(formula = e$fragments.treeDepth1$formula,
                     nPredictedPeaks = sapply(e$fragments.treeDepth1$ms2, function(x){length(unique(x$mZ))}))
    
    load(fn, envir = e)
    X <- e$X.top3.master.intensity
    #Y <- do.call('rbind', do.call('rbind',  e$X.top3.master.intensity.MS2))
    Y <- e$X.top3.master.intensity.MS2
    Y <- merge(Y, pp, by.x='formula0', by.y='formula')
    X$m <- paste(X$file, X$scan)
    Y$m <- paste(Y$file, Y$scan)
    XY <- base::merge(X, Y, by="m")
    
    XY$fragmode <- gsub("uvpd50.00", "uvpd050.00", XY$fragmode)
    XY$fragmode <- gsub("uvpd25.00", "uvpd025.00", XY$fragmode)
    XY
  }
  
  getData <- reactive({
    message("GETDATA BEGIN")
    fn <- file.path(system.file(package = 'uvpd'), "extdata", input$rdata)
    unique(.gd(fn))
  })
  
  
  getThermoUVPD_feb2019 <- reactive({
    
    fn <- file.path(system.file(package = 'uvpd'), "extdata/ThermoUVPD_feb2019.csv")
    ThermoUVPD_feb2019 <- read_csv(fn, na = c("", "#N/A"))
    
    ThermoUVPD_feb2019[ThermoUVPD_feb2019$Compound %in% getData()$Compound,
                       c("Compound", "Bruto formula", "Cluster number", "Group")]
  })
  
  getClsuterIds <- reactive({
    sort(unique(getThermoUVPD_feb2019()$`Cluster number`))
  })
  
  getFilteredData <- reactive({
    DF <- getData()
    if(input$removePC){
      META <- getThermoUVPD_feb2019()
      formula.pc <- META[META$Compound %in% input$compound, 2]
      message(formula.pc)
      DF <- DF[as.character(DF$formula) != as.character(formula.pc), ]
    }

    DF[DF$compound %in% input$compound & DF$ppmerror < as.numeric(input$ppmerror), ]
  })
  
  getFormulaPC <- reactive({
    META <- getThermoUVPD_feb2019()
    formula.pc <- META[META$Compound %in% input$compound, 2]
    formula.pc
  })
  getFilteredClusterData <- reactive({
    DF <- getData()
    if(input$removePC){
      formula.pc <- getFormulaPC()
      message(formula.pc)
      DF <- DF[as.character(DF$formula) != as.character(formula.pc), ]
    }
    S <- getThermoUVPD_feb2019()
    
    compound <- unique(S$Compound[S$`Cluster number` == input$clusterid])
    DF[DF$compound %in% compound & DF$ppmerror < as.numeric(input$ppmerror), ]
  })
  
  getAggregatedData <- reactive({
    DF <- aggregate(intensity ~ mZ * file.y * fragmode * compound * formula * Group * mode,
              data=getFilteredData(), FUN=sum)
    

    DF
  })
  
  getAggregatedClusterData <- reactive({
    aggregate(intensity ~ mZ * file.y * fragmode * compound * formula * Group * mode,
              data=getFilteredClusterData(), FUN=sum)
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
  
  
  getScoreTable <- reactive({
    
    DF <- getFilteredData()
    # DF <- .gd(); DF <- DF[DF$Compound == "Triadimenol" & DF$ppmerror < 10, ]
    
    formula <- intensity ~ file.y * scan.y * fragmode * compound *  Group * mode *  nMs2 * nPredictedPeaks * master.intensity
    
    A.length <- aggregate(formula, data=DF, FUN=length)
    
    
    A.sum <- aggregate(formula, data=DF, FUN=sum)
    A.sum$nAssignedPeaks <- A.length$intensity
    
    A.sum$score1 <- A.length$intensity / A.length$nPredictedPeaks
    A.sum$score2 <- A.sum$nAssignedPeaks / A.sum$nMs2
    A.sum$score3 <- A.sum$intensity / A.sum$master.intensity 
    
    #save(DF, A.sum, file="/tmp/ff.RData")
    A.sum[order(A.sum$fragmode), ]
    
  })
  
  output$tableScore <- renderTable({
    getScoreTable()
  })
  
  
  output$scorePlot <- renderPlot({
    A <- getScoreTable()
    
    op <- par(mfrow=c(1, 3))
    
    hist(A$score1)
    hist(A$score2)
    hist(A$score3)
  })
  
  
  output$score1Plot <- renderPlot({
    
    A <- getScoreTable()
    
    
    
    gp <- ggplot(data = A,
                 aes(x = factor(fragmode, levels = getLevels()),
                     y = score1,
                     fill=reorder(formula, mZ))) +
      geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
      scale_x_discrete(drop=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = cm) +
      facet_wrap(~ compound * Group * mode, scales="free", drop=FALSE)
    
    # gp2 <- ggplot(data=unique(subset(S, select=c('fragmode','formula'))), aes(x=factor(fragmode, levels = sort(unique(DF$fragmode))), fill=(formula))) + ggplot2::geom_bar()
    gp
  }, width=1000, height = 400)
  
  
  
  getTableFreq <- reactive({
    DF <- getAggregatedData()
    print(names(DF))
    rv <- table(unique(subset(DF, select=c('fragmode','formula')))[,1])
    
    rv <- (as.data.frame(rv))
  })
  
  output$ThermoUVPD <- renderTable({
    getThermoUVPD_feb2019()  
  })
  
  output$tableFreq <- renderTable({
    
    getTableFreq()
  })
  
  
  
  output$tableFilteredData <- renderTable({
    
    getFilteredData()
  })
  

  
  output$stackedBarChartText <- renderText({ 
    #x  "You have selected this"
    #DF <- getAggregatedData()
    
    DF <- getThermoUVPD_feb2019()
    
    paste(DF[DF$Compound %in% input$compound, ], col=',')
  })
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  output$stackedBarChart <- renderPlot({
    
    DF <- getAggregatedData()
    
    n <- length(unique(DF$formula))
    
    cm <- gg_color_hue(n)
    if (getFormulaPC() %in% as.character(DF$formula)){
      cm <- c(gg_color_hue(n-1),'grey')
    }
    
    
    gp <- ggplot(data = DF,
                 aes(x = factor(fragmode, levels = getLevels()),
                     y = log(intensity, 10),
                     fill=reorder(formula, mZ))) +
      geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
      scale_x_discrete(drop=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = cm) +
      facet_wrap(~ compound * Group * mode, scales="free", drop=FALSE)
    
    # gp2 <- ggplot(data=unique(subset(S, select=c('fragmode','formula'))), aes(x=factor(fragmode, levels = sort(unique(DF$fragmode))), fill=(formula))) + ggplot2::geom_bar()
    gp
  }, width=1000, height = 400)
  
  output$stackedBarChartGroup <- renderPlot({
    
    gp <- ggplot(data = getAggregatedClusterData(),
                 aes(x = factor(fragmode, levels = getLevels()),
                     y = log(intensity, 10),
                     fill=reorder(formula, mZ))) +
      geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
      scale_x_discrete(drop=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      facet_wrap(~ compound , scales="free", drop=FALSE)
    
    # gp2 <- ggplot(data=unique(subset(S, select=c('fragmode','formula'))), aes(x=factor(fragmode, levels = sort(unique(DF$fragmode))), fill=(formula))) + ggplot2::geom_bar()
    gp
  }, width=1000,height = 400)
  
  
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
                    group=file.x,
                    data=getFilteredData())
  }, width=1000, height = 400)
})
