#!/usr/bin/env Rscript

# Christian Panse <cp@fgcz.ethz.ch>, 2019-03-26
# Christian Panse <cp@fgcz.ethz.ch>, 2019-05-10



makeFragments <- function(){
    library(readr)
    library(uvpd)
    
    ThermoUVPD_feb2019 <- read_csv(file.path(system.file(package = 'uvpd'),
                                             "/extdata/ThermoUVPD_feb2019.csv"))

    idx <- order(ThermoUVPD_feb2019$MONOISOTOPIC_MASS_1, ThermoUVPD_feb2019$SMILES)

    fragments.treeDepth1 <- list(
        SMILES = ThermoUVPD_feb2019$SMILES[idx],
        formula = ThermoUVPD_feb2019$MOLECULAR_FORMULA[idx],
        mass = ThermoUVPD_feb2019$MONOISOTOPIC_MASS_1[idx],
        ms2 = lapply(ThermoUVPD_feb2019$SMILES[idx], uvpd:::getFragments, treeDepth=1)
    )
    
    fragments.treeDepth2 <- list(
        SMILES = ThermoUVPD_feb2019$SMILES[idx],
        formula = ThermoUVPD_feb2019$MOLECULAR_FORMULA[idx],
        mass = ThermoUVPD_feb2019$MONOISOTOPIC_MASS_1[idx],
        ms2 = lapply(ThermoUVPD_feb2019$SMILES, uvpd:::getFragments, treeDepth=2)
    )
   
    save(fragments.treeDepth1,  fragments.treeDepth2,  file="fragments.20200625.RData", compression_level = 9)  
}

summary.analyze.uvpd <- function(S, group="Castell"){
    S <- S[!sapply(S, is.null)]
    
    S <- lapply(S, function(x){aggregate(intensity ~ ., data=x,  FUN=function(x){log(x, 10)})})
    
    S <- lapply(S, function(x){x$intensity.norm <- (x$intensity - mean(x$intensity))/sd(x$intensity); x})
    
    S <- lapply(S, function(x){
        xx <- aggregate(intensity.norm ~ mZ + type + SMILES + formula + scan + nfragments + nMS2 + sumMS2intensities + rawfile + scanTypeFilter + SMILES0 + formula0 + mass + nScans,
                        data = x, FUN=sum)
        xxx <- aggregate(intensity.norm ~ rawfile + formula0 + scanTypeFilter + scan + nfragments + nMS2 + sumMS2intensities, data=xx, FUN=sum)
        xxx
    })
    SS <- do.call('rbind', S)
    SS$group <- group
    SS
    
}

make_X20200612 <- function(){
    library(readr)
    ThermoUVPD_feb2019 <- read_csv(file.path(system.file(package = 'uvpd'),
                                             "/extdata/ThermoUVPD_feb2019.csv"))
    
    X20200612_uvpd <- read_csv("../inst/extdata/20200612-uvpd.csv")
    
    file_group <- read_csv("../inst/extdata/file-group.csv")
    
    FG <- file_group[paste(file_group$file, ".raw",sep='') %in% X20200612_uvpd$filename, ]
    FG$file <- paste(FG$file, '.raw', sep='')
    
    
    GS <- data.frame(group=ThermoUVPD_feb2019$Group, SMILES=ThermoUVPD_feb2019$SMILES)
    GSF <- merge(GS, FG, by='group')
    
    X20200612_uvpd$m <- paste(X20200612_uvpd$filename, X20200612_uvpd$SMILES)
    
    GSF$m <- paste(GSF$file, GSF$SMILES)
    
    
    X20200612 <- merge(GSF, X20200612_uvpd, by='m')
    save(X20200612, file='/tmp/X20200612.RData')
}


make_shiny <- function(){
    library(MsBackendRawFileReader)
    library(protViz)
    
    load(file.path(system.file(package = 'uvpd'), "/extdata/X20200612.RData"))
    X20200612 <- X20200612_uvpd
    
    #X <- X20200612[!(X20200612$Compound=="4-Chlorobenzoic acid" & !grepl("KWR", X20200612$file) & X20200612$Group != "noFrag" ), ]
    
    X20200612$file <- X20200612$filename
    X <-  X20200612
    cc <- expand.grid(as.character(unique(X$file)), as.character(unique(X$SMILES)), as.character(unique(X$fragmode)))
    names(cc) <- c('file', 'SMILES', 'fragmode')
    cc <- paste(cc$file, cc$SMILES, cc$fragmode)
    
    # this is the real world
    X$m <- paste(X$file, X$SMILES, X$fragmode)
    
    library(parallel)
    X.top3.master.intensity <- do.call('rbind', lapply(cc, function(m){
        XX <- X[X$m == m, ]
        XX <- XX[order(XX$master.intensity, decreasing = TRUE), ]
        n <- nrow(XX)
        
        if (nrow(XX) > 2){
            return(XX[1:3, ])}
        else if (n > 0){
            return(XX[1:n, ])}
        else{return(NULL)}
    }))
    


    table(t<-table(paste(X.top3.master.intensity$scan, X.top3.master.intensity$file, X.top3.master.intensity$Compound))>1)
    X.top3.master.intensity <- X.top3.master.intensity[!(X.top3.master.intensity$Compound == "4-Chlorobenzoic acid" & grepl("stds_pos", X.top3.master.intensity$file)), ]
    X.top3.master.intensity <- X.top3.master.intensity[!(X.top3.master.intensity$Compound == "4-Chlorobenzoic acid" & X.top3.master.intensity$Group == "noFrag"), ]
    table(t<-table(paste(X.top3.master.intensity$scan, X.top3.master.intensity$file, X.top3.master.intensity$Compound))>1)


    # sanity check
    # t<-table(paste(X.top3.master.intensity$scan, X.top3.master.intensity$file, X.top3.master.intensity$Compound))
   
    load(file.path(system.file(package = 'uvpd'), "/extdata/fragments.20200625.RData"))
    
    fragmentsPredicted <- as.data.frame(do.call('rbind', lapply(fragments.treeDepth1$ms2, function(x){table(x$type)})))
    
    fragmentsPredicted$formula <- fragments.treeDepth1$formula
    
   
    
    RAWFILEDIR <- '/export/03c3e7cc-5661-4600-94fa-116ccb424918/p2722'
    
    absoluteErrorCutOffInDalton <- 0.01; centroid <- TRUE
    .fragmentMatch <- function(absoluteErrorCutOffInDalton=0.01, centroid=TRUE){
        lapply(sort(unique(X.top3.master.intensity$file)), function(f){
            sn <- X.top3.master.intensity$scan[X.top3.master.intensity$file==f]
            sm <- X.top3.master.intensity$SMILES[X.top3.master.intensity$file==f]
            cm <- X.top3.master.intensity$Compound[X.top3.master.intensity$file==f]
            rawfile <- file.path(RAWFILEDIR, f)
            message(paste("processing rawfile", rawfile, "..."))
            x <- .cnew ("Rawfile", rawfile)
            
            rv <- mapply(function(scan, smile, compound){
                mZ <- x$GetSpectrumMasses(scan)
                intensity <- x$GetSpectrumIntensities(scan)
		message(paste("scan", scan))
                
                if (centroid){
                    cc <- protViz::centroid(mZ, intensity)
                }else{
                    cc <- data.frame(mZ=mZ, intensity=intensity)
                }
                
                DF <- data.frame(mZ=cc$mZ, intensity=cc$intensity)
                
                # computed by MetFrag and stored in the uvpd package
                idx <- which(fragments.treeDepth1$SMILES == as.character(smile))[1]
                insilico <- fragments.treeDepth1$ms2[[idx]]
                insilico <- insilico[!is.na(insilico$mZ),]
                
                # determine best match - assigne in-silico fragment ion to peak
                NN <- findNN(cc$mZ, insilico$mZ)
                absoluteError <- insilico$mZ[NN] - cc$mZ 
                ppmerror <- (abs(absoluteError) /  insilico$mZ[NN]) * 1000000
                
                # filtering but we do not remove yet
                absoluteErrorFilter <- abs(absoluteError) < absoluteErrorCutOffInDalton
                
                DF$nPeaks <- length(cc$mZ)
                DF$nAssignedPeaks <- sum(absoluteErrorFilter)
                DF$intensityCoverage <- round(sum(cc$intensity[absoluteErrorFilter])/sum(cc$intensity),2)
                DF$file <- f
                # TODO(cp)
                # rv$formula0 <- NA
                DF$formula0 <- fragments.treeDepth1$formula[idx]
                DF$formula <- insilico$formula[NN]
                DF$compound <- compound
                DF$type <- insilico$type[NN]
                DF$eps <- absoluteError
                DF$ppmerror <- ppmerror
                DF$scan <- scan
                #rv$SMILE <- smile
                DF[absoluteErrorFilter, ]
            }, scan=sn, smile=sm, compound=cm, SIMPLIFY = FALSE)
            do.call('rbind', rv)
        })
    }
    
    X.top3.master.intensity.MS2 <-do.call('rbind', .fragmentMatch(absoluteErrorCutOffInDalton=1.0))
    save(X.top3.master.intensity, X.top3.master.intensity.MS2, file="uvpd.20200702.RData", compression_level = 9)
}

