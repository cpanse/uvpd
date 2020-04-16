#R

# Christian Panse <cp@fgcz.ethz.ch>

# fgcz-148
# R 3.6.3


library(uvpd)
library(metfRag)

library(readr)

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


rawfiles <- c('/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190213/7stds/stds_pos_neg_MS_highconc_UVPD_100_150_mz50.raw',
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
    '/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/stds_pos_neg_MS_highconc_HCD_mz100-800.raw',
    '/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190212/stds_pos_neg_MS_highconc_HCD_mz100-800.raw',
    '/srv/www/htdocs/p2722/Proteomics/LUMOS_0/Thermo_feb2019/20190214/acQ/stds_pos_neg_MS_highconc_UVPD_50_300_mz50_acQ21.raw',
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


df <- data.frame(SMILES=ThermoUVPD_feb2019$SMILES,
    Group=ThermoUVPD_feb2019$Group,
    exact.mass=m,
    mm = m - 1.007,
    mp = m + 1.007)
