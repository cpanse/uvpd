#R

# Christian Panse <cp@fgcz.ethz.ch>

# fgcz-148
# R 4.0.1

#library(uvpd)
#library(metfRag)

library(readr)

library(MsBackendRawFileReader)
library(protViz)

f <- file.path(system.file(package = 'uvpd'), "/extdata/ThermoUVPD_feb2019.csv")

ThermoUVPD_feb2019 <- read_csv(f)
dim(ThermoUVPD_feb2019)


getMass <- function(smiles = "CC(C)(C)C(O)C(OC1=CC=C(Cl)C=C1)N1C=NC=N1"){
	molecule <- parse.smiles(smiles)[[1]]
	do.isotopes(molecule)
	# determine mass of protonated molecule
	mass <- get.exact.mass(molecule)
	mass
}


m <- sapply(ThermoUVPD_feb2019$SMILES, getMass)

rawfiles <- 
c('/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190212/stds_pos_neg_MS_highconc_HCD_mz100-800.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/7stds/stds_pos_neg_MS_highconc_UVPD_100_150_mz50.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/7stds/stds_pos_neg_MS_highconc_UVPD_200_250_mz50.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/7stds/stds_pos_neg_MS_highconc_UVPD_25_800_mz50.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/7stds/stds_pos_neg_MS_highconc_UVPD_400_500_mz50.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/7stds/stds_pos_neg_MS_highconc_UVPD_50_300_mz50.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/DBPs/DBP_neg_HCD_20_35_60.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/DBPs/DBP_neg_UVPD_100_150.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/DBPs/DBP_neg_UVPD_200_250.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/DBPs/DBP_neg_UVPD_25_800.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/DBPs/DBP_neg_UVPD_400_500.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/DBPs/DBP_neg_UVPD_50_300.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/CastellonStds/CastellonStds_pos_HCD_20_35_60.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/CastellonStds/CastellonStds_pos_UVPD_100_150.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/CastellonStds/CastellonStds_pos_UVPD_200_250.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/CastellonStds/CastellonStds_pos_UVPD_25_800.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/CastellonStds/CastellonStds_pos_UVPD_400_500.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/CastellonStds/CastellonStds_pos_UVPD_50_300.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_HCD_20_35_60_met1.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_100_150_met1.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_100_150_met2.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_200_250_met1.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_200_250_met2.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_25_800_met1.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_25_800_met2.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_400_500_met1.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_400_500_met2.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_50_300_met1.raw',
'/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/KWRstds/KWRstds_UVPD_50_300_met2.raw')

testthat::expect_equal(sum(sapply(rawfiles, file.exists)), length(rawfiles))

df <- data.frame(SMILES=ThermoUVPD_feb2019$SMILES,
    Group=ThermoUVPD_feb2019$Group,
    Compound=ThermoUVPD_feb2019$Compound,
    exact.mass=m,
    mm = m - 1.007,
    mp = m + 1.007)


table(df$Group)


f <- rawfiles[grepl("KWR", rawfiles)][1]

.determineXIC <- function(x, i, ppm = 10, dt = 0.5){
  	pc <- x$GetPrecursorMz(i)
	rtime <-  x$GetRtime()[i] / 60

	rv <- x$GetXIC(pc, ppm, "ms")

	n <- length(rv)

	df <- data.frame(times = rv[seq(1, n, by = 2)],
	  intensities = rv[seq(2, n, by = 2)])

	f <-  rtime - dt  < df$times & df$times < rtime + dt
	x <- df$times[f]
	y <- df$intensities[f]


	return(protViz:::.trapez(x,y))
}


.extractXIC <- function(xic, rtime, dt=0.5){
	stopifnot(length(xic) == length(rtime))
	n <- length(xic)

	rv <- sapply(1:n, function(idx){
		xx <- xic[[idx]]
		rt <- rtime[idx] / 60
		f <-  rt - dt  < xx$times & xx$times < rt + dt
		x <- xx$times[f]
		y <- xx$intensities[f]
		protViz:::.trapez(x,y)
	})
	rv
}



.extractMasterIntensity <- function(x, i, eps=0.01){
  pc.intensity <- NA
  masterScan <- x$GetMasterScans()[i]
  # message(masterScan)

  pc <- x$GetPrecursorMz(i)
  intensity <- x$GetSpectrumIntensities(masterScan)
  mZ <- x$GetSpectrumMasses(masterScan)

  df <- centroid(mZ, intensity)

  idx <- findNN(pc, df$mZ)
  if (abs( df$mZ[idx] - pc) < eps){
	  pc.intensity <- df$intensity[idx]
  }

  pc.intensity
}

.getQuant <- function(i, x){
  stopifnot(!x$IsCentroidScan(i))

  header=x$GetScanFilter(i)
  intensity <- x$GetSpectrumIntensities(i)

  mZ <- x$GetSpectrumMasses(i)
  pc <- x$GetPrecursorMz(i)

  pl.centroid <- centroid(mZ, intensity)
  nMs2 <- nrow(pl.centroid)

  # extract PrecursorMz Intensity

  mZ.lower <- pc - 1.5
  mZ.upper <- mZ.lower + 5.0

  idx.pc <- which(mZ.lower < pl.centroid$mZ & pl.centroid$mZ < mZ.upper)

  pc.mZ <- NA
  pc.intensity <- 0.0
  pc.sum.window.intensity <- 0.0
  tic <- sum(pl.centroid$intensity)

  if (length(idx.pc) > 0){
      #  print(paste("pc", pc, "peaks"))
      #  print(pl.centroid[idx.pc, ])
      #  print(diff(pl.centroid[idx.pc, 'mZ']))
      #  print('#')

	pc.sum.window.intensity <- sum(pl.centroid$intensity[idx.pc])
  	pc.intensity <- max(pl.centroid$intensity[idx.pc])
	pc.mZ <- pl.centroid$mZ[idx.pc][which(pl.centroid$intensity[idx.pc] == pc.intensity)[1]]
  }


 data.frame(tic = tic,
   tic.wopc = tic - pc.sum.window.intensity,
   pc.sum.window.intensity = pc.sum.window.intensity,
   pc.intensity = pc.intensity,
   master.intensity = .extractMasterIntensity(x, i),
  # xic = .determineXIC(x, i),
   pc.mZ = pc.mZ,
   nMs2 = nMs2,
   header=header,
   scan=i)
}


.determineFragMode <- function(x){
	modes <- c('hcd20.00', 'hcd35.00', 'hcd60.00', 'uvpd100.00', 'uvpd150.00', 'uvpd200.00', 'uvpd25.00', 'uvpd250.00', 'uvpd300.00', 'uvpd400.00', 'uvpd50.00', 'uvpd500.00', 'uvpd800.00')

	rv <- rep(NA, length(x))

	for (m in modes){
		rv[grepl(m, x)] <- m
	}

	rv
}

pp <- function(rawfile, df, eps=0.1){
	be <- backendInitialize(MsBackendRawFileReader(), files = rawfile)
	S <- Spectra(be)
	x <- S@backend@rawfileReaderObj[[1]]
	pc <- precursorMz(S) 
	rtime <-  x$GetRtime() 

	rv.mp <- lapply(1:nrow(df), function(df.idx){
	       scan.idx <-  which (abs(pc - df$mp[df.idx]) < eps)
	       if (length(scan.idx) > 0){
		       message(paste("found+", length(scan.idx), df$SMILES[df.idx], '#', df$Group[df.idx], "in file", basename(rawfile)))
		       rv <- lapply(scan.idx, .getQuant, x=x)
		       #print(scan.idx)
		       rv <- do.call('rbind', rv)
		       rv$SMILES <- df$SMILES[df.idx]
		       rv$Compound <- df$Compound[df.idx]
		       rv$Group <- df$Group[df.idx]
		       rv$mode <- 1.007
		       rv$exact.mass <- df$exact.mass[df.idx]
		       rv
	       }
	})
	rv.mm <- lapply(1:nrow(df), function(df.idx){
	       scan.idx <-  which (abs(pc - df$mm[df.idx]) < eps)
	       if (length(scan.idx) > 0){
		       message(paste("found-", length(scan.idx), df$SMILES[df.idx], '#', df$Group[df.idx], "in file", basename(rawfile)))
		       rv <- lapply(scan.idx, .getQuant, x=x)
		       rv <- do.call('rbind', rv)
		       rv$SMILES <- df$SMILES[df.idx]
		       rv$Compound <- df$Compound[df.idx]
		       rv$Group <- df$Group[df.idx]
		       rv$mode <- -1.007
		       rv$exact.mass <- df$exact.mass[df.idx]
		       rv
	       }
	})
	rv.mm <- do.call('rbind', rv.mm)
	rv.mp <- do.call('rbind', rv.mp)
	rv <- rbind(rv.mm, rv.mp)
	rv$filename <- basename(rawfile)
	rv$rtime <- rtime[rv$scan]

	rv$fragmode <- .determineFragMode(rv$header)

	xic <- readXICs(rawfile, masses = rv$exact.mass + rv$mode, tol = 10)
	rv$xic <- .extractXIC(xic, rv$rtime)

	rv
}

message(paste("processing fileidx", fileidx, rawfiles[fileidx], "..."))
uvpdSummary <- do.call('rbind', lapply(rawfiles[fileidx], pp, df=df, eps=0.01))
write.csv(uvpdSummary, file=paste("202006l1-uvpd-summary", fileidx, "csv", sep='.'), row.names = FALSE)

#sessionInfo()

