#R


#' get all fragments of a SMILES code
#'
#' @param smiles 
#'
#' @return a \code{data.frame} containing the SMILES and MH1P charged fragments.
#' @author AB,CP 2019
#' @export getFragments
#' @import rcdk
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

matchFragment <- function(x, y, ...){
  
}
  
