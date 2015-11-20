## Script to compute annual carbon flux
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library(RNetCDF)
library(scales)
library(lubridate)
library(REddyProc)
library (dplyr)
library (plyr)
library(gtools)
library(tidyr)
library (foreach)
library (data.table)

#1. Create a df for all fluxnet sites dataframe

#Open all files
dir <- file.path(path, 'Fluxnet_Data/Daily')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = fLoadFluxNCIntoDataframe(VarList.V.s=c('NEE_f','GPP_f','Reco',
                                                             'NEE_fqcOK', 'precip_f', 'Tair_f','Rg_f'),
                                          FileName.s=list[i], NcPackage.s = 'RNetCDF')
}

#2. Compute dataframe for time since disturbance or stand age analysis
#Open csv file with flux site metadata
Site_Date<-read.csv("Input/Potential_Sites.csv", header = TRUE)

# 2.1 Convert Datetime into a Date object
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]$DateTime=as.Date(Fluxnet_Site[[i]]$DateTime)
}

# 2.2. Remove data before plantation from the list of dataframe

# Convert date of plantation into a date object
Site_Date$Plantation_Date<- as.Date(Site_Date$Plantation_Date)

# Remove data before plantation
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]<-with(Fluxnet_Site[[i]], Fluxnet_Site[[i]][(Fluxnet_Site[[i]]$DateTime > Site_Date$Plantation_Date[i]),])
}

# Remove not complete years of precipitation, temperature and Rg
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]$Tair_f <- with(Fluxnet_Site[[i]], ave(Tair_f, year, FUN = function(x) if(any(is.na(x))) NA else x))
  Fluxnet_Site[[i]]$precip_f <- with(Fluxnet_Site[[i]], ave(precip_f, year, FUN = function(x) if(any(is.na(x))) NA else x))
  Fluxnet_Site[[i]]$Rg_f <- with(Fluxnet_Site[[i]], ave(Rg_f, year, FUN = function(x) if(any(is.na(x))) NA else x))
}

#3. Analyse annual carbon flux

# Compute sum/mean annual +/- sd
Sum_Sd_Flux<-lapply(Fluxnet_Site, function (x) ddply(x, "year",
                    summarize,
                    NEE= sum(NEE_f, na.rm=T),
                    GPP= sum(GPP_f, na.rm=T),
                    Respiration= sum(Reco,na.rm=T),
                    NEP=sum(-NEE_f, na.rm=T),
                    mean_Uncert=mean(NEE_fqcOK, na.rm=T),
                    Annual_Preci= sum(precip_f, na.rm=T),
                    Tair=mean(Tair_f, na.rm=T),
                    Rg=mean(Rg_f, na.rm=T)))

#Add type of ecosystem, type of climate, type of disturbance and year after disturbance for each sites
for (i in seq_along(Sum_Sd_Flux)){
  Sum_Sd_Flux[[i]]$Ecosystem<- Site_Date$Ecosystem[i]
  Sum_Sd_Flux[[i]]$Climate<- Site_Date$Climate[i]
  Sum_Sd_Flux[[i]]$Disturbance<- Site_Date$Type_Disturbance[i]
  Sum_Sd_Flux[[i]]$Stand_Age<- Sum_Sd_Flux[[i]]$year- year(as.Date(Site_Date$Plantation_Date[i]))
  Sum_Sd_Flux[[i]]$GPP_ER<- Sum_Sd_Flux[[i]]$GPP / Sum_Sd_Flux[[i]]$Respiration
  Sum_Sd_Flux[[i]]$NEP_GPP<- Sum_Sd_Flux[[i]]$NEP / Sum_Sd_Flux[[i]]$GPP
  Sum_Sd_Flux[[i]]$Site_ID<- Site_Date$ID[i]
  Sum_Sd_Flux[[i]]$Study<- Site_Date$Study[i]
  Sum_Sd_Flux[[i]]$Stand_Replacement<- Site_Date$Stand_Replacement[i]
  Sum_Sd_Flux[[i]]$Int_Replacement<- Site_Date$Intensity_Replacement[i]
  Sum_Sd_Flux[[i]]$Lat<- Site_Date$x[i]
  Sum_Sd_Flux[[i]]$Long<- Site_Date$y[i]
}

#Combine the flux sites in one dataframe
dfAll_Sites<- do.call("rbind", Sum_Sd_Flux)

# Add climate data to the dataframe
#Au-Tum
dfAll_Sites[6, "Annual_Preci"]<- 1159
dfAll_Sites[6, "Tair"]<- 10.72

#Au-Wac
dfAll_Sites[c(7,9), "Annual_Preci"]<- 1106
dfAll_Sites[c(7,9), "Tair"]<- 12.76

#BE-Bra
dfAll_Sites[18, "Annual_Preci"]<- 743
dfAll_Sites[18, "Tair"]<- 10.01

#BR-cax
dfAll_Sites[c(21,25), "Annual_Preci"]<- 2523
dfAll_Sites[c(21,25), "Tair"]<- 25.99

#CA-CA1
dfAll_Sites[27, "Annual_Preci"]<- 1369
dfAll_Sites[27, "Tair"]<- 9.93

#CA-CA2
dfAll_Sites[36, "Annual_Preci"]<- 1474
dfAll_Sites[36, "Tair"]<- 9.86

#CA-CA3
dfAll_Sites[42, "Annual_Preci"]<- 1676
dfAll_Sites[42, "Tair"]<- 9.94

#CA-Gro
dfAll_Sites[47, "Annual_Preci"]<- 831
dfAll_Sites[47, "Tair"]<- 1.30

#CA-NS1
dfAll_Sites[c(60, 62, 63), "Annual_Preci"]<- 500
dfAll_Sites[c(60, 63), "Tair"]<- -2.89

#CA-NS2
dfAll_Sites[c(64, 68), "Annual_Preci"]<- 500
dfAll_Sites[c(64, 68), "Tair"]<- -2.88

#CA-NS3
dfAll_Sites[c(69, 71, 73), "Annual_Preci"]<- 502
dfAll_Sites[c(69, 73), "Tair"]<- -2.87

#CA-NS4
dfAll_Sites[74, "Annual_Preci"]<- 502
dfAll_Sites[74, "Tair"]<- -2.87

# CA-NS5
dfAll_Sites[c(77, 81), "Annual_Preci"]<- 500
dfAll_Sites[c(77, 81), "Tair"]<- -2.86

# CA-NS6
dfAll_Sites[c(82, 86), "Annual_Preci"]<- 495
dfAll_Sites[c(82, 86), "Tair"]<- -3.08

# CA-NS7
dfAll_Sites[c(87, 90), "Annual_Preci"]<- 483
dfAll_Sites[c(87, 90), "Tair"]<- -3.52

#CA-Qcu
dfAll_Sites[114, "Annual_Preci"]<- 950
dfAll_Sites[114, "Tair"]<- 0.13

#CA-Qfo
dfAll_Sites[120, "Annual_Preci"]<- 962
dfAll_Sites[120, "Tair"]<- -0.36

#CA-SF1
dfAll_Sites[124, "Annual_Preci"]<- 470
dfAll_Sites[124, "Tair"]<- 0.40

#CA-SF2
dfAll_Sites[129, "Annual_Preci"]<- 470
dfAll_Sites[129, "Tair"]<- 0.40

#CA-SJ1
dfAll_Sites[c(133:137), "Annual_Preci"]<- 430
dfAll_Sites[133, "Tair"]<- 0.13

#CA-SJ2
dfAll_Sites[c(138, 139), "Annual_Preci"]<- 430

#CA-SJ3
dfAll_Sites[c(141, 142), "Annual_Preci"]<- 433

#CA-TP1
dfAll_Sites[c(143), "Annual_Preci"]<- 1036

#CA-TP2
dfAll_Sites[c(145:147), "Annual_Preci"]<- 1036
dfAll_Sites[145, "Tair"]<- 8.00

#CA-TP4
dfAll_Sites[148, "Annual_Preci"]<- 1036
dfAll_Sites[148, "Tair"]<-8.00

#CN-Ku1
dfAll_Sites[c(152:153), "Annual_Preci"]<- 249
dfAll_Sites[c(152:153), "Tair"]<- 7.27

#CZ-Bk1
dfAll_Sites[c(154:160), "Annual_Preci"]<- 1026
dfAll_Sites[c(154:160), "Tair"]<- 4.72

#De-Bay
dfAll_Sites[161, "Annual_Preci"]<- 1159
dfAll_Sites[161, "Tair"]<- 5.15

#DE-Tha
dfAll_Sites[c(167, 168), "Annual_Preci"]<- 643

#DK-Sor
dfAll_Sites[183, "Annual_Preci"]<- 573
dfAll_Sites[183, "Tair"]<- 8.03

#FI-Hyy
dfAll_Sites[c(194, 198), "Annual_Preci"]<- 620
dfAll_Sites[194, "Tair"]<- 2.18

#FI-Sod
dfAll_Sites[c(205: 207), "Annual_Preci"]<- 525
dfAll_Sites[207, "Tair"]<- -1.13

#FR-Hes
dfAll_Sites[215, "Annual_Preci"]<- 793
dfAll_Sites[215, "Tair"]<- 9.24

#Il-Yat
dfAll_Sites[c(236:238), "Annual_Preci"]<- 277
dfAll_Sites[236, "Tair"]<- 17.56

#It-Col
dfAll_Sites[c(240:241, 244:248, 250), "Annual_Preci"]<- 971
dfAll_Sites[c(240,248, 250), "Tair"]<- 7.32

#It-Lav
dfAll_Sites[c(251:254, 256), "Annual_Preci"]<- 757

#It-Non
dfAll_Sites[260, "Annual_Preci"]<- 742

#It-Ro1
dfAll_Sites[c(264, 265), "Annual_Preci"]<- 764
dfAll_Sites[264, "Tair"]<- 15.35

#It-Sro
dfAll_Sites[c(276, 278, 279), "Annual_Preci"]<- 898

#JP-Tak
dfAll_Sites[c(284:289), "Annual_Preci"]<- 1024

#NL-Loo
dfAll_Sites[296, "Annual_Preci"]<- 786
dfAll_Sites[296, "Tair"]<- 9.36

#Ru-Fyo
dfAll_Sites[c(312:315, 319), "Annual_Preci"]<- 671
dfAll_Sites[312, "Tair"]<- 4.38

#SE-Fla
dfAll_Sites[c(322:323, 325:328), "Annual_Preci"]<- 616
dfAll_Sites[c(322, 325), "Tair"]<- 0.27

#SE-Nor
dfAll_Sites[c(329:338), "Annual_Preci"]<- 561
dfAll_Sites[c(331, 333:335, 337), "Tair"]<- 5.45

#SE-Sk1
dfAll_Sites[339, "Annual_Preci"]<- 567

#SE-Sk2
dfAll_Sites[340:341, "Annual_Preci"]<- 573

#Uk-Gri
dfAll_Sites[343, "Annual_Preci"]<- 1623

#UK-Ham
dfAll_Sites[354, "Annual_Preci"]<- 829
dfAll_Sites[354, "Tair"]<- 9.38

#US-Blo
dfAll_Sites[c(355:357), "Annual_Preci"]<- 1226
dfAll_Sites[c(355:357), "Tair"]<- 11.09

#US-Fmf
dfAll_Sites[c(373:374), "Annual_Preci"]<- 546
dfAll_Sites[c(373:374), "Tair"]<- 9.50

#US-Fwf
dfAll_Sites[c(375:376), "Annual_Preci"]<- 557
dfAll_Sites[c(375:376), "Tair"]<- 8.40

#US-Ha1
dfAll_Sites[377, "Annual_Preci"]<- 1071
dfAll_Sites[c(377, 379, 380, 388), "Tair"]<- 6.62

#US-Ha2
dfAll_Sites[393, "Tair"]<- 6.56

#US-Lph
dfAll_Sites[c(404, 407), "Annual_Preci"]<- 1071
dfAll_Sites[c(404, 407), "Tair"]<- 6.73

#US-Me1
dfAll_Sites[c(408, 409), "Annual_Preci"]<- 705
dfAll_Sites[c(408, 409), "Tair"]<- 7.88

#US-Oho
dfAll_Sites[415, "Annual_Preci"]<- 849
dfAll_Sites[415, "Tair"]<- 10.10

#US-Sp1
dfAll_Sites[c(417:421), "Annual_Preci"]<- 1310
dfAll_Sites[c(417:421), "Tair"]<- 20.06

#US-Sp2
dfAll_Sites[423, "Annual_Preci"]<- 1314
dfAll_Sites[423, "Tair"]<- 20.07

#US-Syv
dfAll_Sites[439, "Annual_Preci"]<- 826
dfAll_Sites[439, "Tair"]<- 3.81

#US-Umb
dfAll_Sites[440, "Tair"]<- 5.83

#Vu-Coc
dfAll_Sites[c(461, 464), "Annual_Preci"]<- 2763
dfAll_Sites[c(461, 464), "Tair"]<- 25.1

# Remove disturbed years
# sites: AU-Tum (2003= insect outbreak), BE-Bra (2000: thinning), CA-Ca1 (1997: disturbance occured), CA-Ca2 (2000/2001: forest not planted yet), CA-NS4 (2002: disturbance occured), CA-NS6 (2001), 
#CA-OJP (1999: disturbance occured), CA-SJ2 (2003\2004: outlier for GPP values), NL-Loo (2006: storm), SE-Nor (2005: thinning occured this year),UK-Gri (2006:Thinning)
# US-fmf (2006: thinnning), US-Nc2 (2006: insect attack), US-SP1 (2005: storm), US-SP3 (2000: storm), US-Wcr (2002: insect attack)
dfAll_Sites<-dfAll_Sites[-c(3, 13, 27, 36, 37, 74, 82, 107, 119, 138, 139, 306, 338, 352, 374, 413, 422, 430, 448),]

# Remove high gap filled fraction
dfAll_Sites<- dfAll_Sites[dfAll_Sites$mean_Uncert>0.80,]

# Remove years missing in each flux site
dfAll_Sites = dfAll_Sites[!is.na(dfAll_Sites$mean_Uncert),]

# Restructure dataframe
dfAll_Sites<-gather(dfAll_Sites, variable, value, -Annual_Preci, -year, 
                    -Ecosystem, -Climate, -Disturbance,
                    -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
                    -Tair, -Rg, -Study, -Lat, -Long )

#Reoder column
dfAll_Sites<- dfAll_Sites[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "value", "Annual_Preci", 
                            "Tair","Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")]

# Rename head column
colnames(dfAll_Sites)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", 
                            "Tair", "Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")

# Append forest global database
Glob_Forest_df<-read.csv("Input/Forest_Global_Database.csv", header = TRUE)
Glob_Forest_df$GPP_ER<-Glob_Forest_df$GPP / Glob_Forest_df$Reco
Glob_Forest_df$NEP= -Glob_Forest_df$NEE
Glob_Forest_df$NEP_GPP<-Glob_Forest_df$NEP / Glob_Forest_df$GPP
Glob_Forest_df<-gather(Glob_Forest_df, variable, value, -Annual_Preci, -year, 
                    -Ecosystem, -Climate, -Disturbance,
                    -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
                    -Tair, -Rg, -Study, -Lat, -Long )
Glob_Forest_df<- Glob_Forest_df[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "value", "Annual_Preci", 
                            "Tair", "Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")]
colnames(Glob_Forest_df)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", 
                         "Tair", "Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")
dfAll_Sites<-rbind(dfAll_Sites, Glob_Forest_df)

# Save dataframe in a rds file
saveRDS(dfAll_Sites, file="Output/df_Annual_Flux.rds")
