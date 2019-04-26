#!/usr/bin/env Rscript

# Christian Panse <cp@fgcz.ethz.ch>, 2019-03-26



makeFragments <- function(){
    library(readr)
    library(uvpd)
    
    ThermoUVPD_feb2019 <- read_csv(file.path(system.file(package = 'uvpd'),
                                             "/extdata/ThermoUVPD_feb2019.csv"))

    fragments.treeDepth1 <- lapply(ThermoUVPD_feb2019$SMILES, getFragments, treeDepth=1)
    # fragments.treeDepth2 <- lapply(ThermoUVPD_feb2019$SMILES, getFragments, treeDepth=2)
    save(fragments.treeDepth1, file="../inst/extdata/fragments.RData")  
}
