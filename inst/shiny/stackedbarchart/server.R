#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(DT)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    
    getData0 <- reactive({
        e <- new.env()
        load(input$rdata, envir = e)
        return(e$X.top3.master.intensity)
    })
    
    getData <- reactive({
        load("uvpd.20200522.RData")
        
        
        X <- getData0()
        Y <- do.call('rbind', do.call('rbind', X.top3.master.intensity.MS2))
    
        X$m <- paste(X$file, X$scan)
        Y$m <- paste(Y$file, Y$scan)
        XY <- base::merge(X, Y, by="m")
       
        XY$fragmode <- gsub("uvpd50.00", "uvpd050.00", XY$fragmode)
        XY$fragmode <- gsub("uvpd25.00", "uvpd025.00", XY$fragmode)
        XY
        
    })
    
    getFilteredData <- reactive({
        DF <- getData()
        
        DF[DF$compound %in% input$compound & DF$ppmerror < as.numeric(input$ppmerror), ]
    })
    
    getAggregatedData <- reactive({
        aggregate(intensity ~ mZ * file.y * fragmode * compound * formula * Group * mode, data=getFilteredData(), FUN=sum)
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
    output$selectRData <- renderUI({
        selectInput("rdata", "rdata", c('uvpd.20200522.RData','uvpd.20200611.RData'), multiple = FALSE, selected = 'uvpd.20200522.RData')
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
    })
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
    }, height=250)
 
    output$barchart <- renderPlot({
        lattice::barchart(log(intensity) ~ fragmode | compound,
                          group=mode ,
                          data=aggregate(intensity ~ mode * fragmode * compound, data=getFilteredData(), FUN=sum),
                          scales = list(x = list(rot = 45)),
                          layout=c(6,7))
    }, width=1000, height=1000)
})
