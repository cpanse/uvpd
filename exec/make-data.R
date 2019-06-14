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

 