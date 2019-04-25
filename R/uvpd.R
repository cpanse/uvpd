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
#' @examples
#' df <- getFragments(treeDepth = 1)
#' plot(table(df$MH1P))
getFragments <-function(smiles="CC(C)(C)C(O)C(OC1=CC=C(Cl)C=C1)N1C=NC=N1", ...){
  molecule <- parse.smiles(smiles)[[1]]
  
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
                   type=rep('M1P', length(mZ)),
                   SMILES=as.character(sapply(fragments, rcdk::get.smiles)),
                   formula=as.character(sapply(fragments, function(x){rcdk::get.mol2formula(x)@string}))
                  )
 
  
  df <- unique(df)
  
  df.MH1P <- df
  df.MH1P$mZ <- df$mZ  + 1.00727646677
  df.MH1P$type <- "MH1P"
  
  
  
  df.M2H1P <- df.MH1P
  df.M2H1P$mZ <- df.MH1P$mZ + 1.00782504
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


#' Title
#'
#' @param x 
#' @param eps.mZ 
#' @param eps.rt 
#' @param treeDepth 
#' @param ... 
#'
#' @return
#'

#' @importFrom rawDiag readScans read.raw
.computeMatch <- function(x, eps.mZ = 0.05, eps.rt=0.5, treeDepth=1,scanTypeFilter="", ...){
  
  rawfile <- x$rawfile[1]
  RAW <- read.raw(rawfile)
  RAW <- RAW[RAW$MSOrder == "Ms2", ]
  
  rv <- lapply(1:length(x$mass), function(i){
    
    mz <- x$mass[i]
    smiles <- as.character(x$smiles0[i])
    rt.max <- x$rt.max[i] 
    
    filter <- grepl(scanTypeFilter, RAW$ScanType) &
      mz - eps.mZ < RAW$PrecursorMass & 
      RAW$PrecursorMass <  mz + eps.mZ &
      rt.max - eps.rt < RAW$StartTime &
      RAW$StartTime <= rt.max + eps.rt
    
    (sn <- RAW$scanNumber[filter])
    
    if(length(sn) > 0){
      MS2Scans <- .filterMS2Scans(readScans(rawfile, sn))
      
      df.frags <- getFragments(smiles, treeDepth=treeDepth)
      df.frags <- na.omit(df.frags)
      
      rv.match <- matchFragment(MS2Scans, fragments = df.frags, FUN=sum)
      msg <- paste(i, "id\n",
                   length(MS2Scans), "#MS2 scans\n", 
                   nrow(df.frags), "#in-silico fragments\n",
                   nrow(rv.match), "#found matches\n", 
                   basename(rawfile), "rawfile\n", 
                   smiles, "SMILES\n",
                   mz, "mass\n\n", sep="\t")
    
      
      if(!is.null(rv.match)){
        message(msg)
        rv.match$rawfile <- rawfile
	rv.match$scanTypeFilter <- scanTypeFilter
        rv.match$mass <- mz
        rv.match$rt.max <- rt.max
        rv.match$SMILES0 <- smiles
        rv.match$nMS2Scans <- length(sn)
        
        return(rv.match)}
    }
    return (NULL)
  })
  do.call('rbind', rv)
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
  
  
  mZ.hits <-fragments$mZ[idx.NN[hits]]
  intensity.hits <- x$intensity[hits]
  
  hit.idx <- idx.NN[hits]
  
  hit.df <- fragments[idx.NN[hits], ]
  hit.df$intensity <- intensity.hits
  
  
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
  if (nrow(S)>0){
    
  
  S <- aggregate(intensity ~ mZ + type + SMILES + formula, 
                 data=S,
                 FUN=FUN)
  
  S <- S[order(S$mZ),]
  rownames(S) <- 1:nrow(S)
  return(S)
  }
  NULL
}
  

.getScanType <- function(rawfile){
  RAW <- read.raw(rawfile)
  unique(unlist(lapply(RAW$ScanType[RAW$MSOrder=='Ms2'], function(x){paste("@",strsplit(strsplit(x, "@")[[1]][2], " ")[[1]][1],sep='')})))
}


#' Analyse ms2 specs of a given vector of smiles
#'
#' @param rawfile 
#' @param smiles 
#' @param mZoffset 
#'
#' @return 
#' @export analyze 
analyze <- function(rawfile, smiles,  mZoffset=1.007){
  X <- uvpd:::.extractXICs(rawfile, smiles, mZoffset=mZoffset)
  
  M <- lapply(.getScanType(rawfile), function(stf){uvpd:::.computeMatch(X, scanTypeFilter=stf)})
  
  M <- do.call('rbind', M)
  #class(M) <- list('data.frame', 'uvpdres')
  M
}
