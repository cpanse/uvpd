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
  
  
  df <- data.frame(MH1P=sapply(fragments, rcdk::get.exact.mass) + 1.0072,
                   SMILES=as.character(sapply(fragments, rcdk::get.smiles)))
  
  # TODO(cp)
  # formal=as.character(sapply(fragments,rcdk::get.formal.charge)))
  
  # write.csv(file='peaklist1.csv',
  #  data.frame(mZ=MS2Scans[[1]]$mZ,
  #     intensity=MS2Scans[[1]]$intensity),
  # row.names = FALSE)
  
  idx <- order(df$MH1P)
  df[idx, ]
}

matchFragment <- function(x, y, ...){
  
}
  