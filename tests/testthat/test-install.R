#R
  
#context("package install")
library(testthat)

test_that("check if Java is performing", {
library(metfRag)

  smiles <- "CC(C)(C)C(O)C(OC1=CC=C(Cl)C=C1)N1C=NC=N1"
  molecule <- parse.smiles(smiles)[[1]]

#calculate the fragments
fragments <- frag.generateFragments(molecule, 1)

 expect_true(length(fragments) == 319)


 do.isotopes(molecule)
 expect_equal(get.exact.mass(molecule), 295.1088, tolerance = 0.001)

})

