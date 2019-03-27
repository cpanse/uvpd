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
  
  #length(fragments)
  
  for (i in fragments) {
    try(do.aromaticity(i))
    try(do.typing(i))
    try(do.isotopes(i))
  }
  
 
  # TODO(ab): do we have to remove e-?
  df <- data.frame(M1P=sapply(fragments, rcdk::get.exact.mass) ,
                   SMILES=as.character(sapply(fragments, rcdk::get.smiles)),
                   formula=as.character(sapply(fragments, function(x){rcdk::get.mol2formula(x)@string}))
                  )
 
  # consider only mols having a weight of more than 50Da.
  df <- df[df$M1P > 50, ]
 
  df$MH1P  <- df$M1P  + 1.00727646677
  df$M2H1P <- df$MH1P + 1.00782504
  
  
  df$M1P[!grepl("[NSOPI]|Cl|Br", df$formula)] <- NA

  # write.csv(file='peaklist1.csv',
  #  data.frame(mZ=MS2Scans[[1]]$mZ,
  #     intensity=MS2Scans[[1]]$intensity),
  # row.names = FALSE)
  
  idx <- order(df$MH1P)
  df[idx, ]
}

matchFragment <- function(x, y, ...){
  
}
  
