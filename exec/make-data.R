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
        ms2 = lapply(ThermoUVPD_feb2019$SMILES[idx], getFragments, treeDepth=1)
    )
    
    fragments.treeDepth2 <- list(
        SMILES = ThermoUVPD_feb2019$SMILES[idx],
        formula = ThermoUVPD_feb2019$MOLECULAR_FORMULA[idx],
        mass = ThermoUVPD_feb2019$MONOISOTOPIC_MASS_1[idx],
        ms2=lapply(ThermoUVPD_feb2019$SMILES, getFragments, treeDepth=2)
    )
   
    save(fragments.treeDepth1,  fragments.treeDepth2,  file="fragments.RData", compression_level = 9)  
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

pp <- c('stds_pos_neg_MS_highconc_HCD_mz100-800',
        'stds_pos_neg_MS_highconc_UVPD_100_150_mz50',
        'stds_pos_neg_MS_highconc_UVPD_200_250_mz50',
        'stds_pos_neg_MS_highconc_UVPD_50_300_mz50',
        'stds_pos_neg_MS_highconc_UVPD_400_500_mz50',
        'stds_pos_neg_MS_highconc_UVPD_25_800_mz50',
        'DBP_neg_UVPD_25_800',
        'DBP_neg_UVPD_50_300',
        'DBP_neg_UVPD_100_150',
        'DBP_neg_UVPD_200_250',
        'DBP_neg_UVPD_400_500',
        'DBP_neg_HCD_20_35_60',
        'CastellonStds_pos_UVPD_25_800',
        'CastellonStds_pos_UVPD_50_300',
        'CastellonStds_pos_UVPD_100_150',
        'CastellonStds_pos_UVPD_200_250',
        'CastellonStds_pos_UVPD_400_500',
        'CastellonStds_pos_HCD_20_35_60',
        'KWRstds_UVPD_25_800_met1',
        'KWRstds_UVPD_25_800_met2',
        'KWRstds_UVPD_50_300_met1',
        'KWRstds_UVPD_50_300_met2',
        'KWRstds_UVPD_100_150_met1',
        'KWRstds_UVPD_100_150_met2',
        'KWRstds_UVPD_200_250_met1',
        'KWRstds_UVPD_200_250_met2',
        'KWRstds_UVPD_400_500_met1',
        'KWRstds_UVPD_400_500_met2',
        'KWRstds_HCD_20_35_60_met1')


extractMs2Feature <- function(){
    library(uvpd)
    load(file.path(system.file(package = 'uvpd'), "/extdata/fragments.RData"))
    
    rawfiles <- scan(file.path(system.file(package = 'uvpd'),
                               "/extdata/rawfiles_all.txt"), what = character())
    unlist(lapply(pp, function(p){idx <- grep(p, rawfiles); rawfiles[idx]}))
    rawfiles <- rawfiles[grepl("Castell|(KWR|stds)|DBPs", rawfiles)]
    rawfiles <- rawfiles[!grepl("AcX|blank|Blank", rawfiles)]
    
    #' # DBPs neg
    #' # Castell pos
    #' # KWR std. pos and neg
    #' 

    rawfiles.Castell <- rawfiles[grepl("Castell", rawfiles)]
    S.Castell <- lapply(rawfiles.Castell[2:7], FUN=analyze, fragments = fragments.treeDepth1, mZoffset = +1.007)
    
    rawfiles.DBPs <- rawfiles[grepl("DBPs", rawfiles)]
    S.DBP <- lapply(rawfiles.DBPs[c(1, 3:6)], FUN=analyze, fragments = fragments.treeDepth1, mZoffset = +1.007, eps.rt = 1)
    
    rawfiles.KWRneg <- rawfiles[grepl("(KWR|stds)", rawfiles) &  grepl("(neg)", rawfiles)]
    rawfiles.KWRpos <- rawfiles[grepl("(KWR|stds)", rawfiles) &  grepl("(pos)", rawfiles)]
    S.KWRpos <- lapply(rawfiles.KWRpos[c(1:5, 8:16)], FUN=analyze, fragments = fragments.treeDepth1, mZoffset = +1.007)
    
    S.KWRneg <- lapply(rawfiles.KWRneg[c(1:16)], FUN=function(x){
        try({rv <- analyze(x, fragments = fragments.treeDepth1, mZoffset = -1.007); return(rv)}); NULL})
    
    #ms2Features <- do.call('rbind', list()
    uvpd.Castell <- S.Castell
    uvpd.DBP <- S.DBP
    uvpd.KWRpos <- S.KWRpos
    uvpd.KWRneg<-S.KWRneg
    save(uvpd.Castell, uvpd.DBP, uvpd.KWRpos,uvpd.KWRneg, 
         file="../inst/extdata/ms2assignments.RData")  
}
