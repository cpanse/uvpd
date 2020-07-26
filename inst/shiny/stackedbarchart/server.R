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

negativeIonTypePattern <- c("M-", "M-2H-", "M-H-")
positiveIonTypePattern <- c("M1P", "M2H1P",  "MH1P")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$selectRData <- renderUI({
    selectInput("rdata", "rdata", c('uvpd_20200522.RData',
                                    'uvpd_20200612.RData',
                                    'uvpd_20200626.RData',
                                    'uvpd_20200702.RData'),
                 multiple = FALSE, selected = 'uvpd_20200702.RData')
  })
  
  output$selectFragmentsRData <- renderUI({
    selectInput("rdata", "rdata", c('uvpd_20200522.RData',
                                    'uvpd_20200612.RData'),
                multiple = FALSE, selected = 'uvpd_20200612.RData')
  })
  
  
  output$selectCluster <- renderUI({
    selectInput("clusterid", "cluster ID", getClsuterIds(),
                multiple = FALSE, selected = getClsuterIds()[1])
  })
  
  getPredictedFragments <- reactive({
    e <- new.env()
    load(file.path(system.file(package = 'uvpd'), "extdata",
                   "fragments.20200625.RData"), envir = e)
    e$fragments.treeDepth1
  })
  
  # ease for debugging
  .gd <- function(fn=file.path(system.file(package = 'uvpd'), "extdata",
                               'uvpd_20200626.RData')){
    e <- new.env()
    
    
    load(file.path(system.file(package = 'uvpd'), "extdata",
                   "fragments.20200625.RData"), envir = e)
    
    #e$fragments.treeDepth1 <- getPredictedFragments()
    pp <- data.frame(formula = e$fragments.treeDepth1$formula,
                     nPredictedPeaks = sapply(e$fragments.treeDepth1$ms2,
                                              function(x){length(unique(x$mZ))}))
    
    load(fn, envir = e)
    X <- e$X.top3.master.intensity
    #Y <- do.call('rbind', do.call('rbind',  e$X.top3.master.intensity.MS2))
    Y <- e$X.top3.master.intensity.MS2
    Y <- merge(Y, pp, by.x='formula0', by.y='formula')
    message("DEBUG")
    message(dim(X))
    message(dim(Y))
    X$file <- X$filename
    
    X$m <- paste(X$file, X$scan)
    Y$m <- paste(Y$file, Y$scan)
    XY <- base::merge(X, Y, by="m")
    
    XY$fragmode <- gsub("uvpd50.00", "uvpd050.00", XY$fragmode)
    XY$fragmode <- gsub("uvpd25.00", "uvpd025.00", XY$fragmode)
    
    
    f.Benzocaine <- XY$compound %in% c('Benzocaine', '3-Nitrophenol', '4-Chlorobenzoic acid') & grepl("_met1", XY$file.y)
    XY <- XY[!f.Benzocaine, ]
    
    # TODO(cp): sanity check
    dd.compound <- c('4-Nitrocatechol', '2-Nitrohydroquinone', '4-Nitro-1,3-benzenediol',
            '2-Hydroxy-3-nitrobenzoic acid', '4-Hydroxy-3-nitrobenz',
            '2-Hydroxy-5-nitrobenzoic acid')
    
    dd.formula0 <- c('C7H5NO5')
    message(dim(XY))
    XY <- XY[!(XY$compound %in% dd.compound | (XY$formula0 %in% dd.formula0 & XY$Compound != "Benzoic acid, 2-hydroxy-4-nitro-")), ]
    
    # filter 2
    
    dd.compound0 <- "2-Amino-3-nitrobenzoic acid"
    XY <- XY[!(XY$compound %in% dd.compound0 & (grepl("^KWR", XY$filename) | XY$Group == "KWR"))  , ]
    
    dd.compound1 <- "4-Nitroanthranilic acid"
    XY <- XY[!(XY$compound %in% dd.compound1 & (grepl("^DBP", XY$filename) | XY$Group == "DBP"))  , ]
    
    XY
  }
  #---- getData ----
  getData <- reactive({
    message("GETDATA BEGIN")
    fn <- file.path(system.file(package = 'uvpd'), "extdata", input$rdata)
    
    #M <- M[, c("Compound", "Cas nr")]
    DF <- unique(.gd(fn))
    #DF<-unique(merge(DF, M, by="Compound"))
    #message("GETDATA END")
    DF
    
  })
  
  
  getThermoUVPD_feb2019 <- reactive({
    
    fn <- file.path(system.file(package = 'uvpd'), "extdata/ThermoUVPD_feb2019.csv")
    ThermoUVPD_feb2019 <- read_csv(fn, na = c("", "#N/A"))
    
    DF <- ThermoUVPD_feb2019[, c("Compound", "Bruto formula", "Cluster number", "Group", "Cas nr")]
    DF[ThermoUVPD_feb2019$Compound %in% getData()$Compound, ]
  })
  
  getClsuterIds <- reactive({
    sort(unique(getThermoUVPD_feb2019()$`Cluster number`))
  })
  
  #---- getFilteredData ----
  getFilteredData <- reactive({
    DF <- getData()
    if(input$removePC){
      
      formula.pc <- getFormulaPC()
      message(paste("pc formula", formula.pc))
      
      DF <- DF[!as.character(DF$formula) %in% as.character(formula.pc), ]
    }

    
    filter <- DF$ppmerror < as.numeric(input$ppmerror) | abs(DF$eps) < as.numeric(input$epserror)
    DF <- DF[DF$compound %in% input$compound & filter, ]
    
    if (input$negativeIonType){
      message("negative mode")
      DF <- DF[DF$type %in% negativeIonTypePattern, ]
    }
    else{
      message("positive mode")
      DF <- DF[DF$type %in% positiveIonTypePattern, ]
    }
    
    drops <- c("SMILES", "file.x", "m")
    unique(DF[, !(names(DF) %in% drops)])
    
  })
  
  getFormulaPC <- reactive({
    META <- getThermoUVPD_feb2019()
    formula.pc <- META$`Bruto formula`[META$Compound %in% input$compound]
    
    unique(formula.pc)
  })
  
  getFilteredClusterData <- reactive({
    DF <- getData()
    
    if(input$removePC){
      formula.pc <- getFormulaPC()
      message(paste("pc formula", formula.pc))
      DF <- DF[as.character(DF$formula) != as.character(formula.pc), ]
    }
    S <- getThermoUVPD_feb2019()
    
    compound <- unique(S$Compound[S$`Cluster number` == input$clusterid])
    DF <- DF[DF$compound %in% compound & DF$ppmerror < as.numeric(input$ppmerror) & abs(DF$eps) < as.numeric(input$epserror), ]
    
    if (input$negativeIonType){
      message("negative mode")
      DF <- DF[DF$type %in% negativeIonTypePattern, ]
    }
    else{
      message("positive mode")
      DF <- DF[DF$type %in% positiveIonTypePattern, ]
    }
    
    drops <- c("SMILES", "file.x", "m")
    unique(DF[, !(names(DF) %in% drops)])
    
  })
  
  getAggregatedData <- reactive({
    
    if (nrow(getFilteredData())>0){
      aggregate(intensity ~ mZ * file.y * fragmode * compound * formula * Group * mode * type,
                data=getFilteredData(), FUN=sum)}
    
    else{NULL}
  })
  
  getAggregatedClusterData <- reactive({
    aggregate(intensity ~ mZ * file.y * fragmode * compound * formula * Group * mode,
              data=getFilteredClusterData(), FUN=sum)
  })
  
  getLevels <- reactive({
    DF <- getData()
    if(nrow(DF)>0){
      
    
    sort(unique(DF$fragmode))
    }else{NULL}
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
  
  #------ fitdistr-----
  output$distPlot <- renderPlot({
    par(mfrow=c(3, 2))
    hist(getFilteredData()$ppmerror, sub=input$compound)
    hist(abs(getFilteredData()$eps), sub=input$compound)
    
    if(require(MASS)){
      hist(getFilteredData()$ppmerror, probability = TRUE, sub=input$compound)
      fit <- fitdistr(getFilteredData()$ppmerror, densfun="normal")
      curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)
      
      hist(abs(getFilteredData()$eps), probability = TRUE, sub=input$compound)
      fit <- fitdistr(abs(getFilteredData()$eps), densfun="normal")
      curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)
    }
    
    if(require(MASS)){
      
      DF <- getData()
      filter <- DF$ppmerror < as.numeric(input$ppmerror) | abs(DF$eps) < as.numeric(input$epserror)
      DF <- DF[filter, ]
      
      hist(DF$ppmerror, probability = TRUE, sub='all compunds')
      fit <- fitdistr(DF$ppmerror, densfun="normal")
      curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)
      
      hist(abs(DF$eps), probability = TRUE, sub='all compunds')
      fit <- fitdistr(abs(DF$eps), densfun="normal")
      curve(dnorm(x, fit$estimate[1], fit$estimate[2]), col="red", lwd=2, add=T)
    }
    
    return
  }, height=600)
  
  
  #---- score ----
  # score1: experimental fragments matched/ theoretically possible (@CP for new dataset)
  # score2: experimental fragments matched/ all fragments in spectrum
  # score3: experimental matched fragment intensities / master.intensity  
  getScoreTable <- reactive({
    
    DF <- getData()
    
    if(input$removePC){
      formula.pc <- getFormulaPC()
      message(paste("pc formula", formula.pc))
      DF <- DF[as.character(DF$formula) != as.character(formula.pc), ]
    }
    
    filter <- DF$ppmerror < as.numeric(input$ppmerror) | abs(DF$eps) < as.numeric(input$epserror)
    DF <- DF[filter, ]
    
    formula <- intensity ~ file.y * scan.y * fragmode * compound *  formula0 * mode *  nMs2 * nPredictedPeaks * master.intensity
    
    A.length <- aggregate(formula, data=DF, FUN=length)
    
    
    A.sum <- aggregate(formula, data=DF, FUN=sum)
    A.sum$nAssignedPeaks <- A.length$intensity
    
    A.sum$score1 <- A.length$intensity / A.length$nPredictedPeaks
    A.sum$score2 <- A.sum$nAssignedPeaks / A.sum$nMs2
    A.sum$score3 <- A.sum$intensity / A.sum$master.intensity 
    
    #save(DF, A.sum, file="/tmp/ff.RData")
    DF <- A.sum[order(A.sum$compound, A.sum$fragmode),
                c('compound', 'formula0', 'file.y', 'scan.y','mode', 'fragmode', 'nMs2', 'nPredictedPeaks', 'master.intensity','score1', 'score2', 'score3')]
    #save(DF, file="/tmp/score.RData")
    M <- getThermoUVPD_feb2019()[, c("Compound", "Cas nr")]
    DF <- merge(M, DF, by.y='compound', by.x="Compound")
    DF <- unique(DF)
    names(DF)[1] <- 'compound'
    DF
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadScores <- downloadHandler(
    filename = function() {
      paste("uvpd-shiny-application-scores", ".csv", sep = "")
    },
    content = function(file) {
      
      write.csv(getScoreTable(), file, row.names = FALSE)
    }
  )
  
  # Downloadable csv of selected dataset ----
  output$downloadFreq <- downloadHandler(
    filename = function() {
      paste("uvpd-shiny-application-freq", ".csv", sep = "")
    },
    content = function(file) {
      
      write.csv(getTableFreqAll(), file, row.names = FALSE)
    }
  )
  
  getFilteredScoreTable <- reactive({
    
    DF <- getFilteredData()
    
    formula <- intensity ~ file.y * scan.y * fragmode * compound *  Group * mode *  nMs2 * nPredictedPeaks * master.intensity
    
    A.length <- aggregate(formula, data=DF, FUN=length)
    
    
    A.sum <- aggregate(formula, data=DF, FUN=sum)
    A.sum$nAssignedPeaks <- A.length$intensity
    
    A.sum$score1 <- A.length$intensity / A.length$nPredictedPeaks
    A.sum$score2 <- A.sum$nAssignedPeaks / A.sum$nMs2
    A.sum$score3 <- A.sum$intensity / A.sum$master.intensity 
    
    
    DF<-A.sum[order(A.sum$fragmode), c('compound', 'mode', 'fragmode', 'score1', 'score2', 'score3')]
   
    DF
  })
  
  output$tableScore <- DT::renderDataTable({
    datatable(getScoreTable())
  })
  
  
  output$scorePlot <- renderPlot({
    A <- getFilteredScoreTable()
    
    if (input$negativeIonType){
      message("negative mode")
      A <- A[A$mode < 0, ]
    }
    else{
      message("positive mode")
      A <- A[A$mode > 0, ]
    }
    if(nrow(A) > 1){
      score1 <- A[, c('fragmode', 'score1')] ; score1$score <-"score1"; names(score1) <- c('fragmode','value','score')
      score2 <- A[, c('fragmode', 'score2')] ; score2$score <-"score2"; names(score2) <- c('fragmode','value','score')
      score3 <- A[, c('fragmode', 'score3')] ; score3$score <-"score3"; names(score3) <- c('fragmode','value','score')
      lattice::dotplot(value ~ fragmode | score, data=do.call('rbind', list(score1, score2, score3)),
                       scales = "free", layout=c(1,3), index.cond = list(rev(1:3)))
    }else{
      plot(0,0, type='n', axes=FALSE);text(0,0, "no meaningful data", cex=3)
    }
  } , width=1000)
  
  output$score1Plot <- renderPlot({
    lattice::xyplot(score1 ~ as.factor(fragmode) | mode, group=compound,
                    #auto.key = list(space = "right"),
                    data=getScoreTable(), type='b', layout=c(1,2))
  } , width=1000)
  
  output$score2Plot <- renderPlot({
     lattice::xyplot(score2 ~ as.factor(fragmode) | mode, group=compound,
                     #auto.key = list(space = "right"),
                     data=getScoreTable(), type='b', layout=c(1,2))
  } , width=1000)
  
  output$score3Plot <- renderPlot({
    lattice::xyplot(score3 ~ as.factor(fragmode) | mode, group=compound,
                    auto.key = list(space = "right"),
                    data=getScoreTable(), type='b', layout=c(1,2))
  } , width=1000, height = 640)
  
  
  output$scorePlot11 <- renderPlot({
    
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
  
  .freqstat <- function(DF, compound){
    DF <- DF[DF$compound == compound, ];
    
    if (unique(DF$mode)<0){
      message("negative mode")
      DF <- DF[DF$type %in% negativeIonTypePattern, ]
    }
    else{
      message("positive mode")
      DF <- DF[DF$type %in% positiveIonTypePattern, ]
    }
    
    rv <- aggregate(formula ~ compound * fragmode,
                    data = unique(DF[,c('compound', 'fragmode','formula')]),
                    FUN = length)
    as.data.frame(rv)
  }
  # ==== getTableFreqAll =====
  
  getTableFreqAll <- reactive({
    message("getTableFreqAll")
    DF <- getData()
    filter <- DF$ppmerror < as.numeric(input$ppmerror) | abs(DF$eps) < as.numeric(input$epserror)
    DF <- DF[filter, ]
    if(input$removePC){
      
      formula.pc <- getFormulaPC()
      message(paste("pc formula", formula.pc))
      
      DF <- DF[!as.character(DF$formula) %in% as.character(formula.pc), ]
    }
    rv <- do.call('rbind', lapply(unique(DF$compound), .freqstat, DF=DF))
    names(rv) <- c('compound', 'fragmentation mode', 'frequency')
    rv
  })
  
  output$TableFreqAll <- DT::renderDataTable({
    getTableFreqAll()
  })
  
  getTableFreq <- reactive({
    DF <- getAggregatedData()
  
    rv <- table(unique(subset(DF, select=c('fragmode','formula')))[,1])
    
    as.data.frame(rv)
  })
  
  output$ThermoUVPD <- renderTable({
    getThermoUVPD_feb2019()[, c(1,2,3,4,5)]  
  })
  
  output$tableFreq <- renderTable({
    rv <- getTableFreq()
    names(rv) <- c('fragmentation mode', 'frequency')
    rv
  })
  
  output$tableFilteredData <- DT::renderDT({
    
    getFilteredData()
  })
  
  output$stackedBarChartText <- renderTable({ 
    #x  "You have selected this"
    #DF <- getAggregatedData()
    
    DF <- getThermoUVPD_feb2019()
    tt <- DF[DF$Compound %in% input$compound, ]
    tt[1, c(1,2,3,5)]
  })
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  #---- stacked barchart ion type
  output$stackedBarChartIonType <- renderPlot({
    
    DF <- getAggregatedData()
    
    n <- length(unique(DF$formula))
    
    cm <- gg_color_hue(n)
    if (getFormulaPC() %in% as.character(DF$formula)){
      cm <- c(gg_color_hue(n-1),'grey')
    }
    
    
    gp <- ggplot(data = DF,
                 aes(x = factor(fragmode, levels = getLevels()),
                     y = log(intensity, 10),
                     fill=reorder(type, mZ))) +
      geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
      scale_x_discrete(drop=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = cm) +
      facet_wrap(~ compound * mode, scales="free", drop=FALSE)
    
    # gp2 <- ggplot(data=unique(subset(S, select=c('fragmode','formula'))), aes(x=factor(fragmode, levels = sort(unique(DF$fragmode))), fill=(formula))) + ggplot2::geom_bar()
    gp
  }, width=1000, height = 400)
  
  #---- stacked barchart fragment ion----
  output$stackedBarChart <- renderPlot({
    
    DF <- getAggregatedData()
    
    n <- length(unique(DF$formula))
    
    cm <- gg_color_hue(n)
    if (getFormulaPC() %in% as.character(DF$formula)){
      cm <- c(gg_color_hue(n-1), 'grey')
    }
    
    
    gp <- ggplot(data = DF,
                 aes(x = factor(fragmode, levels = getLevels()),
                     y = log(intensity, 10),
                     fill=reorder(formula, mZ))) +
      geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
      scale_x_discrete(drop=FALSE) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_fill_manual(values = cm) +
      facet_wrap(~ compound * mode, scales="free", drop=FALSE)
    
    # gp2 <- ggplot(data=unique(subset(S, select=c('fragmode','formula'))), aes(x=factor(fragmode, levels = sort(unique(DF$fragmode))), fill=(formula))) + ggplot2::geom_bar()
    gp
  }, width=1000, height = 400)
  
  output$stackedBarChartGroup <- renderPlot({
    if (nrow(getAggregatedClusterData()) > 0){

      gp <- ggplot(data = getAggregatedClusterData(),
                   aes(x = factor(fragmode, levels = getLevels()),
                       y = log(intensity, 10),
                       fill=reorder(formula, mZ))) +
        geom_bar(stat="identity", position = position_stack(reverse = FALSE)) +
        scale_x_discrete(drop=FALSE) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        facet_wrap(~ mode * compound  , scales="free", drop=TRUE)
      
      # gp2 <- ggplot(data=unique(subset(S, select=c('fragmode','formula'))), aes(x=factor(fragmode, levels = sort(unique(DF$fragmode))), fill=(formula))) + ggplot2::geom_bar()
      gp }
    else{''}
  }, width=1000,height = 400)
  
  #------- ms2 -------
  
  output$tableMS2 <- DT::renderDataTable({ 
    DF <- getFilteredData()
    
    DF[, c('mZ', 'intensity', 'formula', 'compound', 'type',	'eps', 	'ppmerror', 'scan.y', 'fragmode')]
  }, options=list(pageLength = 10))
  
  
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
    if (nrow(getData())>0){
      
      
      lattice::xyplot(tic ~ master.intensity | fragmode * Compound,
                      subset = Compound %in% input$compound,
                      group=file.y,
                      data=getData())
    }else{NULL}
  }, width=1000, height = 400)
  
  #---- summary ----
  
  .summary <- function(){
    DF <- getData()
    DF.filtered <- getFilteredData()
    
    rv <- data.frame(
      row.names=c("number of rawfiles", 
                  "number of scans",
                  "number of compounds",
                  "length(unique($fragmode))"),
      all = c(length(unique(DF$file.y)),
              length(unique(paste(DF$scan.y, DF$file.y))),
              length(unique(DF$Compound)),
              length(unique(DF$fragmode))),
      filtered = c(
        length(unique(DF.filtered$file.y)),
        length(unique(paste(DF.filtered$scan.y, DF$file.y))),
        length(unique(DF.filtered$Compound)),
        length(unique(DF.filtered$fragmode)))
    )
  }
  #---- in-silico predicted fragment ions ----
  output$summary <- renderTable({ 
    .summary()
  }, rownames = TRUE)
  
  output$tableinSilicoFragmentIon <- DT::renderDataTable({ 
    DF <- getThermoUVPD_feb2019()
    
    formula=DF$`Bruto formula`[DF$Compound == input$compound]
    
    fragments <- getPredictedFragments()
    
    idx <- which(fragments$formula == formula)
    
    fragments$ms2[[idx]]
    
  }, rownames = TRUE)
  
  
  #---- sessionInfo ----
  output$citation <- renderPrint({
    capture.output(citation('uvpd'))
  })
  
  output$sessionInfo <- renderPrint({
    capture.output(sessionInfo())
  })
  
})
