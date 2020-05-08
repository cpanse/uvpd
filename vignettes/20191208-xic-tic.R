#R

library(uvpd)

library(metfRag)
library(rawDiag)

#load(file.path(system.file(package = 'uvpd'), "/extdata/fragments.RData"))

#rawfiles <- scan(file.path(system.file(package = 'uvpd'),
#                           "/extdata/rawfiles.txt"), what = character())
#rawfiles <- rawfiles[grepl("Castell|KWR|DB", rawfiles)]
(rawfile <- file.path(Sys.getenv('HOME'), "Downloads/p2722",
  "stds_pos_neg_MS_highconc_UVPD_50_300.raw"))


smiles <- "CC(C)(C)C(O)C(OC1=CC=C(Cl)C=C1)N1C=NC=N1"
molecule<-parse.smiles(smiles)[[1]]
do.isotopes(molecule)

# determine mass of protonated molecule
(mZ <- get.exact.mass(molecule) + c(1.007))

# extract ion chromatogram
XIC <- readXICs(rawfile, mZ, tol=10)
plot(XIC)

S <- read.raw(rawfile)
idxMs2 <- S$MSOrder=="Ms2"
idx <- abs(S[idxMs2, ]$PrecursorMass - mZ) < 0.1


sn <- which(idxMs2)[idx][which(max(S$TIC[which(idxMs2)[idx]]) == S$TIC[which(idxMs2)[idx]])]


ms2<-readScans(rawfile, which(idxMs2))
plot(tic1<-S$TIC[idxMs2], tic2<-sapply(ms2, function(x)sum(x$intensity)))

if(sum(idx)>0){
    plot(S[idxMs2[idx], c('StartTime', 'TIC')])

    abline(v = (XIC[[1]]$times[max(XIC[[1]]$intensities) == XIC[[1]]$intensities]), col='red')
}else{
    plot(0,0); text(0,0,"no data found", cex=4)
}
# determine mass of protonated molecule
(mZ <- get.exact.mass(molecule) + c(-1.007))


# extract ion chromatogram
XIC <- readXICs(rawfile, mZ, tol=10)
plot(XIC)


idxMs2 <- S$MSOrder=="Ms2"
idx <- abs(S[idxMs2, ]$PrecursorMass - mZ) < 0.1

if(sum(idx)>0){
    plot(S[idxMs2[idx], c('StartTime', 'TIC')])

    abline(v = (XIC[[1]]$times[max(XIC[[1]]$intensities) == XIC[[1]]$intensities]), col='red')
}else{
    plot(0,0); text(0,0,"no data found", cex=4)
}


library(readr)
library(uvpd)
f <- file.path(system.file(package = 'uvpd'), "/extdata/ThermoUVPD_feb2019.csv")
ThermoUVPD_feb2019 <- read_csv(f)
dim(ThermoUVPD_feb2019)


.xicPlot <- function(xic, RAWFILEMETADATA, idx){
        if(sum(idx) > 0){
            plot(RAWFILEMETADATA[idxMs2[idx], c('StartTime', 'TIC')])
            
            abline(v = (xic$times[max(xic$intensities) == xic$intensities]), col='red')
        }else{
            plot(0,0); text(0,0,"no data found", cex=4)
        }
}

.scanTypeParser <- function(ScanType){
  as.character(sapply(ScanType,
    function(x){strsplit(strsplit(x, "@")[[1]][2], " ")[[1]][1]}))
}   

.aggTICbasePeakIntensity <- function(S){
	if (nrow(S)==0) return(NULL)
	
	S$n <- 1
	df <- merge(aggregate( BasePeakIntensity ~ ScanType2, data=S, FUN=sum), aggregate(TIC ~ ScanType2, data=S, FUN=sum))
	df <- merge(aggregate( n ~ ScanType2, data=S, FUN=sum),df)
	df
}


.determineNegPos <- function(SMILES, mZ){
    molecule <- parse.smiles(as.character(SMILES))[[1]]
    do.isotopes(molecule)
    (mZ.cand <- get.exact.mass(molecule) + c(-1.007, 1.007))


    if  (abs(mZ - mZ.cand[1]) < 0.1){
	    return("-1")
    }else{
	    return("1")
    }


}

.xicCheck <- function(SMILES, rawfile, RAWFILEMETA, plot = FALSE, tol=0.1){
    
    molecule <- parse.smiles(as.character(SMILES))[[1]]
    do.isotopes(molecule)
    
    # determine mass of protonated molecule
    (mZ <- get.exact.mass(molecule) + c(-1.007, 1.007))
    
    #idxMs2 <-  RAWFILEMETA$MSOrder == "Ms2"

    S <- RAWFILEMETA[RAWFILEMETA$MSOrder == "Ms2",]
    S$ScanType2 <- .scanTypeParser(S$ScanType)
    
    idxList <- lapply(mZ, function(mass){
        which(abs(S[, 'PrecursorMass'] - mass) < tol)
    })
    
    rv <- lapply(idxList, function(x){.aggTICbasePeakIntensity(S[x,])})
    df <- lapply(1:length(idxList), function(i){

	    if (is.null(nrow(rv[[i]]))){}else{
	    rv[[i]]$mZ <- mZ[i]; 
    	    XIC <- readXICs(rawfile, mZ[i], tol=10)
	    rv[[i]]$XIC.max <- max(XIC[[1]]$intensities)
	    }
    	    rv[[i]]
	})

    
    df <- do.call('rbind', df)
    if (is.null(df))return(NULL)
    df$SMILES <- rep(as.character(SMILES), nrow(df))
    
    df$rawfile <- basename(rawfile)
    df
} 

(rawfile <- file.path(Sys.getenv('HOME'), "Downloads/p2722",
  "stds_pos_neg_MS_highconc_UVPD_50_300.raw"))

S <- read.raw(rawfile)

.xicCheck(SMILES =  ThermoUVPD_feb2019[48, 'SMILES'],
          rawfile = rawfile,
          RAWFILEMETA = S, plot=TRUE)


.xicCheck("CC(C)(C)C(O)C(OC1=CC=C(Cl)C=C1)N1C=NC=N1", 
          rawfile = rawfile, RAWFILEMETA = S, plot=TRUE)




.forAllRawfiles <- function(rawfile){
	S<- read.raw(rawfile)
    rv <- lapply(1:nrow(ThermoUVPD_feb2019), function(i){
    .xicCheck(SMILES =  ThermoUVPD_feb2019[i, 'SMILES'],
          rawfile = rawfile,
          RAWFILEMETA = S, plot=FALSE)
    })

        rv <- do.call('rbind', rv)
	    rv
}


f <- list.files( p<- file.path(Sys.getenv('HOME'), "Downloads/p2722"))

rv <- lapply(file.path(p, f), .forAllRawfiles)
allXic <- do.call('rbind', rv)


allXic$posNegMode <- sapply(1:nrow(allXic), function(i){.determineNegPos(allXic$SMILES[i], allXic$mZ[i])})

write_csv(allXic, path ="~/Desktop/p2722-allXic.csv")
save(allXic, file ="~/Desktop/p2722-allXic.RData")


pdf("~/Desktop/allXIC.pdf", 10,50);
dotplot(log(TIC) ~ ScanType2| SMILES, group=posNegMode, data=allXic,  scales = list(x = list(rot = 45)), layout=c(1,50));
dotplot(log(XIC.max) ~ ScanType2| SMILES, group=posNegMode, data=allXic,  scales = list(x = list(rot = 45)), layout=c(1,50));
dotplot(log(BasePeakIntensity) ~ ScanType2| SMILES, group=posNegMode, data=allXic,  scales = list(x = list(rot = 45)), layout=c(1,50));
dotplot((log(TIC)-log(XIC.max)) ~ ScanType2| SMILES, group=posNegMode, data=allXic,  scales = list(x = list(rot = 45)));
dev.off()
