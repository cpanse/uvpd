#R

f <- "stds_pos_neg_MS_highconc_UVPD_25_800_mz50.raw"
scan <- 1822

smiles <- "CC(C)(C)C(O)C(OC1=CC=C(Cl)C=C1)N1C=NC=N1"

formula <- "C14H18ClN3O2"

monomass <- 295.1088

mZ <- c(70.0390191982914, 74.0257514976886, 74.030434906838, 99.079684387397, 112.049218374013, 112.058273268925, 112.063100974066, 138.082174410196, 141.01000916825, 147.481643294789, 147.511671809448, 147.535027974609, 168.111623595859, 178.225542641855, 179.474373320018, 227.082535652107, 229.59835616243, 296.114679153925, 296.192315831908, 297.118660247796)
intensity <- c(17320.648803732, 500.627417750204, 745.263681649581, 1431.251431685, 102894.04907875, 948.457372731486, 629.2177057038, 1265.68594958243, 1752.41709132132, 2482.33613059091, 3585.58422323315, 1122.60063388647, 11393.8744496568, 17166.5375230878, 1653.21615357621, 6282.32124009174, 2543.97133922519, 1482661.44825281, 4068.55430343169, 9666.12433057973)

plot(mZ, intensity, type='h', main=formula, sub=smiles, lwd=3)

################################
library(metfRag)
molecule <- parse.smiles(smiles)[[1]]


matching.fragments<-frag.generateMatchingFragments(molecule, mZ, monomass, mzabs = 0.01, mzppm = 10.0, posCharge = TRUE, ionMode = 1, treeDepth = 2)
for (i in matching.fragments) {
    try(do.aromaticity(i)) 
    try(do.typing(i))
    try(do.isotopes(i))
}

matching.ions <- sapply(matching.fragments, rcdk::get.exact.mass)
matching.smiles <- as.character(sapply(matching.fragments, rcdk::get.smiles))
matching.formula <- as.character(sapply(matching.fragments, function(x){rcdk::get.mol2formula(x)@string}))

text(matching.ions, rep(0.5*max(intensity), length(matching.ions)), matching.formula, cex=1, srt=90)
     

     