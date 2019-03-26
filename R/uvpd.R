#R


#' get all fragments of a SMILES code
#'
#' @param smiles 
#'
#' @return \code{data.frame} containing the SMILES and fragments
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
  
  
  df <- data.frame(MH1P=sapply(fragments, rcdk::get.exact.mass) + 1.0072,
                   SMILES=as.character(sapply(fragments, rcdk::get.smiles)))
  df
}

matchFragment <- function(x, y, ...){
  
}
  