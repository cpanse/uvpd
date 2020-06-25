#R


#' get all fragments of a SMILES code
#'
#' @param smiles 
#'
#' @return a \code{data.frame} containing the SMILES and MH1P charged fragments.
#' @author AB,CP 2019
#' @export getFragments
#' @importFrom  rcdk get.smiles get.mol2formula
#' @importFrom metfRag frag.generateFragments
#' @seealso \code{exec/make-data.R}
#' @references \itemize{
#' \item \url{https://cran.r-project.org/package=rcdk}
#' \item \url{https://github.com/ipb-halle/MetFragR}
#' }
#' @examples
#' df <- getFragments(treeDepth = 1)
#' plot(table(df$MH1P))
getFragments <-function(smiles="CC(C)(C)C(O)C(OC1=CC=C(Cl)C=C1)N1C=NC=N1", ...){
  molecule <- parse.smiles(smiles)[[1]]
  
  
  me <- 0.00054858026 
  mH <- 1.00782504
  
  #calculate the fragments
  fragments <- frag.generateFragments(molecule, ...)
  
  for (i in fragments) {
    try(do.aromaticity(i))
    try(do.typing(i))
    try(do.isotopes(i))
  }
  
  # TODO(ab): do we have to remove e-?
  # M1P
  mZ <- sapply(fragments, rcdk::get.exact.mass)
  
  df <- data.frame(mZ=mZ,
                   type=rep('mZ', length(mZ)),
                   SMILES=as.character(sapply(fragments, rcdk::get.smiles)),
                   formula=as.character(sapply(fragments, function(x){rcdk::get.mol2formula(x)@string}))
                  )
  df <- unique(df)
  
  df.M1P <- df
  df.M1P$mZ <- df$mZ - me
  df.M1P$type <- "M1P"
  
  
  # M-
  df.M1N <- df
  df.M1N$mZ <-  df$mZ + me
  df.M1N$type <- "M-"
  
  # M-H-
  df.M1HN <- df
  df.M1HN$mZ <- df$mZ - mH + me
  df.M1HN$type <- "M-H-"
  
  # "M-2H-"
  df.M2HN <- df
  df.M2HN$mZ <- df$mZ - mH -mH + me
  df.M2HN$type <- "M-2H-"
    
  df.MH1P <- df
  df.MH1P$mZ <- df$mZ + mH - me
  df.MH1P$type <- "MH1P"
  
  df.M2H1P <- df.MH1P
  df.M2H1P$mZ <- df.MH1P$mZ + mH
  df.M2H1P$type <- "M2H1P"
  
  df <- do.call('rbind', list(df, df.MH1P, df.M2H1P))
  
  
  # consider only mols having a weight of more than 50Da.
  df <- df[df$mZ > 50, ]
  
  # filter M1Ps
  df[(!grepl("[NSOPI]|Cl|Br", df$formula) & df$type == "M1P"), 'mZ'] <- NA
  
  # write.csv(file='peaklist1.csv',
  #  data.frame(mZ=MS2Scans[[1]]$mZ,
  #     intensity=MS2Scans[[1]]$intensity),
  # row.names = FALSE)
  
  idx <- order(df$mZ)
  df[idx, ]
}

.filterMS2Scans <- function(y, ...){
  lapply(y,
         function(x){
           rv <- x; 
           idx <- x$intensity > 0; 
           x$mZ <- x$mZ[idx]; 
           x$intensity <- x$intensity[idx]; 
           x})
}


# INPUT rawfile, vector of mass
.extractXICs <- function(rawfile, smiles, mZoffset=1.007, method='rt'){
  
  molecules <- lapply(smiles, function(x){
    mol <- parse.smiles(x)[[1]]; 
    do.isotopes(mol); mol
  })
  
  mass <-  sapply(molecules, FUN=get.exact.mass) + mZoffset
  
  XIC <- readXICs(rawfile, mass, tol = 30)
  XIC <- lapply(1:length(XIC),
                function(i){x<-XIC[[i]]; x$smiles <- smiles[i]; x})
  
  # filter
  XIC <- lapply(XIC, function(x){
    if (max(x$intensities) < 1000 || 
        length(x$intensities) < 6 || 
        5 * median(x$intensities) > max(x$intensities)){
      return (NULL)
    }
    
    x
  })
  
  XIC <- XIC[which(!sapply (XIC, is.null))]
  
  if(method=='plot'){
    gp <- plot.XICs(XIC) + facet_wrap(~ mass, scales = "free") + labs(title = rawfile)
    # pdf(paste("/tmp/uvpd--", basename(rawfile),".pdf", sep=''), 19,12)
    print(gp)
  }else{
    rv <- lapply(XIC, function(x){
      data.frame(smiles0=x$smiles,
                 mass=x$mass,
                 rt.max=x$times[which(max(x$intensities) == x$intensities)[1]])
    })
    rv <- do.call('rbind', rv)
    rv$rawfile <- rawfile
    rv
  }
  # dev.off()
}

.matchFragment <- function(x, fragments, errorCutOff=0.001, plot=FALSE){
  
  if (length(x$mZ) < 1){
    warning("no mZ values.")
    return(NULL)
  }
  # compute match
  idx.NN <- findNN(x$mZ, fragments$mZ)
  
  error <- abs(x$mZ - fragments$mZ[idx.NN])
  
  hits <- (error < errorCutOff)
  if (sum(hits)==0){return(NULL)}
  mZ.hits <-fragments$mZ[idx.NN[hits]]
  intensity.hits <- x$intensity[hits]
  
  hit.idx <- idx.NN[hits]
  
  hit.df <- fragments[idx.NN[hits], ]
  hit.df$intensity <- intensity.hits
  hit.df$scan <- x$scan
  hit.df$mZerror <- (x$mZ - fragments$mZ[idx.NN])[hits]
  
  # for scoring
  hit.df$nfragments <- length(fragments$mZ)
  hit.df$nMS2 <- sum(x$intensity > 0)
  hit.df$sumMS2intensities <- sum(x$intensity)
  
  if(plot){
    plot(x$mZ, x$intensity,
         type = "h",
         main = x$title,
         xlab = "m/Z", 
         log='y',
         col='grey',
         # ylim=c(-max(x$intensity), max(x$intensity)),
         ylab = "intensity", 
         sub = x$scanType)
    
    
    abline(v = mZ.hits, col=rgb(0.8,0.3,0.3, alpha = 0.4))
    axis(3, mZ.hits, round(mZ.hits,2))
    
    
    
    text(fragments$mZ[hit.idx],
         rep(max(x$intensity)/2, length(hit.idx)),
         paste(df.frags$formula[hit.idx], df.frags$type[hit.idx], round(df.frags$mZ[hit.idx],2)),
         srt=90, pos=rep(c(2,1), length(hit.idx)/2), cex=0.4)
    
  }
  hit.df
}

#' Compute match between in-silico fragments and MS2 scans
#'
#' @param MS2Scans list of MS2 scans.
#' @param fragments \code{data.frame} of in-silico fragments.
#' @param plot \code{TRUE} or \code{FALSE}
#' @param errorCutOff mZ error cut-off.
#' @param FUN function used by aggregate.
#' @param ... 
#'
#' @return a aggregate data frame
#' @export
#'
#' @examples
matchFragment <- function(MS2Scans, fragments, plot=FALSE, errorCutOff = 0.001, FUN=sum, ...){
    
    rv <- lapply(MS2Scans, FUN = .matchFragment, fragments = fragments, errorCutOff = errorCutOff)
    S <- do.call('rbind', rv)
    try({
        
        
        if (nrow(S) > 0){
            
            
            #S <- aggregate(intensity ~ mZ + type + SMILES + formula, 
            #               data=S,
            #              FUN=FUN)
            
            #attr(S, 'formula') <- "intensity ~ mZ + type + SMILES + formula"
            S <- S[order(S$mZ),]
            rownames(S) <- 1:nrow(S)
            return(S)
        }
    }, silent = TRUE)
    NULL
}
  

.getScanType <- function(rawfile){
  RAW <- read.raw(rawfile)
  unique(unlist(lapply(RAW$ScanType[RAW$MSOrder=='Ms2'], function(x){paste("@",strsplit(strsplit(x, "@")[[1]][2], " ")[[1]][1],sep='')})))
}


.assignAPEX <- function(rawfile, fragments, mZoffset=1.007, method='rt', tolppm=30){
  
  XIC <- readXICs(rawfile, fragments$mass + mZoffset, tol = tolppm)
  n <- length(fragments$mass)
  stopifnot(length(XIC) == n)
  
  fragments$rt.max <- rep(NA, n)
  fragments$intensity.max <-rep(NA, n)
  fragments$rawfile <- rawfile
  
  for (i in 1:n){
    try({
      
      if ("intensities" %in% names(XIC[[i]]) & length(XIC[[i]]$intensities) > 0){
        intensity.max <- max(XIC[[i]]$intensities, na.rm = TRUE)
        
        if (intensity.max < 1000 || 
            length(XIC[[i]]$intensities) < 6 || 
            5 * median(XIC[[i]]$intensities) > intensity.max){
          
        }else{
          
          fragments$rt.max[i] <- XIC[[i]]$times[which(intensity.max == XIC[[i]]$intensities)[1]]
          fragments$intensity.max[i] <- intensity.max
        }
        
      }
    })
  }
  
  fragments
}

.assignMatchedFragmentIons <- function(fragments, mZoffset=1.007, eps.mZ = 0.05,
                                       eps.rt = 0.5, scanTypeFilter = "", errorCutOff = 0.001){
  n <- length(fragments$mass)
  
  RAW <- read.raw(fragments$rawfile)
  RAW <- RAW[RAW$MSOrder == "Ms2", ]
  
  res <- lapply(which(!is.na(fragments$rt.max)), function(i){
    #print(i)
    mz <- fragments$mass[i] + mZoffset
    rt.max <- fragments$rt.max[i] 
    
    filter <- grepl(scanTypeFilter, RAW$ScanType) &
      mz - eps.mZ < RAW$PrecursorMass & 
      RAW$PrecursorMass <  mz + eps.mZ &
      rt.max - eps.rt < RAW$StartTime &
      RAW$StartTime <= rt.max + eps.rt
    
    
    if (sum(filter) > 0){
      scans <- RAW$scanNumber[filter]
      
      MS2Scans <- uvpd:::.filterMS2Scans(readScans(fragments$rawfile, scans))
      
      df <- matchFragment(MS2Scans, fragments = na.omit(fragments$ms2[[i]]), errorCutOff=errorCutOff)
      df$rawfile <- fragments$rawfile
      df$scanTypeFilter <- scanTypeFilter
      df$SMILES0 <- fragments$SMILES[i]
      df$formula0 <- fragments$formula[i]
      df$mass <- fragments$mass[i]
      df$nScans <- length(scans)
      df
    }
  })
}

#' Perfom analyses run
#'
#' @param rawfile filepath to a Thermo Fisher rawfile
#' @param mZoffset 
#' @param itol fragment ion tolerance; default is 1mDa.
#' @param eps.mZ pre-cursor mass tolerance in Da; default is set to 50mDa.
#' @param eps.rt  retention time window; default is 0.5 minutes.
#' @param fragments in-silico computed fragment ions of smile codes. see \code\{\link{getFragments}}. 
#'
#' @author Christian Panse <cp@fgcz.ethz.ch>, 2019
#' @return returns a \code{data.frame} object
#' @export analyze summary.uvpd
#' @aliases summary.uvpd uvpd analyse
#' 
#' @details 
#' INPUT:
#' \itemize{
#' \item in-silico fragment ion spectra of a given set of SMILES 
#' \item ThermoFisher raw file
#' }
#' 
#' STEP1: 
#'  \itemize{
#'  \item dermines APEX(MS1) of each molecule mass return rt for each given mass; 
#' }
#' 
#' STEP2: 
#' \itemize{
#' \item extract MS2 at a given rt and precursor window
#' }
#' 
#' STEP3: 
#' \itemize{
#' \item compute match between in-silico and MS2 fragment specs 
#' }
#' 
#' OUTPUT: 
#' \itemize{
#' \item a \code{data.frame} object having the column names
#' \code{c('mZ', 'type', 'SMILES', 'formula', 'intensity', 'rawfile',
#'   'scanTypeFilter', 'SMILES0', 'formula0', 'mass', 'n'}.
#' }
#' 
#' STATISTICS:
#'   see examples
#' 
#' @examples 
#' library(uvpd)
#' 
#' # load in-silico fragments
#' load(file.path(system.file(package = 'uvpd'), "/extdata/fragments.RData"))
#' 
#' # load Thermo Fisher rawfiles
#' rawfiles <- scan(file.path(system.file(package = 'uvpd'),
#'   "/extdata/rawfiles.txt"), what = character())
#'   
#' # filter 
#' rawfiles <- rawfiles[grepl("Castell|(KWR|stds)|DBPs", rawfiles)]
#' rawfiles <- rawfiles[!grepl("AcX|blank|Blank", rawfiles)]
#' 
#' # DBPs neg
#' # Castell pos
#' # KWR std. pos and neg
#' 
#' \dontrun{
#' S1 <- analyze(rawfiles[20], fragments.treeDepth1)
#' 
#' S2 <- lapply(rawfiles[c(15:20)], analyze, fragments=fragments.treeDepth1)
#' }
#' 
#' \dontrun{
#' 
#' 
#' rawfile24 <- file.path(Sys.getenv('HOME'), "Downloads/CastellonStds_pos_HCD_20_35_60.raw")
#' rawfile29 <- file.path(Sys.getenv('HOME'), "Downloads/CastellonStds_pos_UVPD_50_300.raw")
#' 
#' rawfile24 <- rawfiles[24]
#' rawfile29 <- rawfiles[29]
#'
#' S24 <- analyze(rawfile24, fragments = fragments.treeDepth1)
#' S29 <- analyze(rawfile29, fragments = fragments.treeDepth1)
#' 
#' S24.log <- aggregate(intensity ~ ., data=S24,  FUN=function(x){log(x, 10)})
#' S24.log$intensity.norm <- (S24.log$intensity - mean(S24.log$intensity ))/ sd(S24.log$intensity )
#' 
#' S29.log <- aggregate(intensity ~ ., data=S29,  FUN=function(x){log(x, 10)})
#' S29.log$intensity.norm <- (S29.log$intensity - mean(S29.log$intensity ))/ sd(S29.log$intensity )
#' 
#' S <- rbind(S24, S29)
#' S.log <- rbind(S24.log, S29.log)
#' 
#' boxplot(intensity.norm ~ rawfile, data=S.log )
#' # STATISTICS:
#' 
#' # merge profile data
#' SS.log <- aggregate(intensity ~ mZ + type + SMILES + formula + scan + nfragments + nMS2 + sumMS2intensities + rawfile + scanTypeFilter + SMILES0 + formula0 + mass + nScans, data = S.log, FUN=sum)
#' 
#' # sum intensities
#' SSS.log<-aggregate(intensity ~ rawfile + formula0 + scanTypeFilter + scan + nfragments + nMS2 + sumMS2intensities, data=SS.log, FUN=sum)
#' 
#' library(lattice)
#' 
#' histogram(~intensity/sumMS2intensities | scanTypeFilter, data=SSS.log, type='count')
#' histogram(~intensity/nfragments | scanTypeFilter, data=SSS, type='count')
#' }
analyze <- function(rawfile, fragments, mZoffset = 1.007, itol = 0.001,
                    eps.mZ = 0.05, eps.rt = 0.5){
  
  # 1. extract APEX for a given set of precomputed SMILES fragments
  X <- uvpd:::.assignAPEX(rawfile, fragments,  mZoffset=mZoffset, tolppm=30)
  
  # 2. help function to extract all possible MS2 scanTypes of a rawfile
  scanTypes <- uvpd:::.getScanType(rawfile)
  
  # 3. compute match between in-silico fragments and extracted MS2 scans
  XX <- lapply(scanTypes, function(st){
      uvpd:::.assignMatchedFragmentIons(X, scanTypeFilter = st,
                                 mZoffset=mZoffset,
                                 eps.mZ = eps.mZ,
                                 eps.rt = eps.rt,
                                 errorCutOff=itol)})
  
  # does some list cosmetics
  XXX <- do.call('rbind', lapply(XX, function(x){ do.call('rbind',x[sapply(x, length) == 16])}))
  
  class(XXX) <- c(class(XXX), "uvpd")
  
  XXX
}

summary.uvpd <- function(object, ...){
  cat("\n\nrawfile:\n\t")
  print(table(object$rawfile))
   
  cat("\n\nformula for aggregation:\n\t")
  cat(attr(object, 'formula'))
  
  cat("\n\n")
  print(aggregate(intensity ~ scanTypeFilter, data=object, FUN=sum))
  print(aggregate(intensity ~ scanTypeFilter, data=object, FUN=quantile))
  cat("\n")
  print(aggregate(intensity ~ scanTypeFilter, data=object, FUN=length))
}
