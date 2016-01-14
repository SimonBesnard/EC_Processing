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

# Flag years with uncomplete data collected
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]$NEE_f <- with(Fluxnet_Site[[i]], ave(NEE_f, year, FUN = function(x) if(any(is.na(x))) NA else x))
  Fluxnet_Site[[i]]$GPP_f <- with(Fluxnet_Site[[i]], ave(GPP_f, year, FUN = function(x) if(any(is.na(x))) NA else x))
  Fluxnet_Site[[i]]$Reco <- with(Fluxnet_Site[[i]], ave(Reco, year, FUN = function(x) if(any(is.na(x))) NA else x))
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

#CA-NS1
dfAll_Sites[c(62,63), "Annual_Preci"]<- 500

#CA-NS3
dfAll_Sites[71, "Annual_Preci"]<- 502

# CA-NS5
dfAll_Sites[77, "Annual_Preci"]<- 500

#CA-SJ1
dfAll_Sites[133, "Annual_Preci"]<- 430

#CA-SJ2
dfAll_Sites[138, "Annual_Preci"]<- 430

#CA-TP1
dfAll_Sites[143, "Annual_Preci"]<- 1036

#CA-TP2
dfAll_Sites[c(145:147), "Annual_Preci"]<- 1036

#CZ-Bk1
dfAll_Sites[c(154:157, 159:160), "Annual_Preci"]<- 1026

#FI-Sod
dfAll_Sites[c(212: 213), "Annual_Preci"]<- 525

#Il-Yat
dfAll_Sites[c(244:245), "Annual_Preci"]<- 277

#It-Col
dfAll_Sites[c(251:253), "Annual_Preci"]<- 971

#It-Lav
dfAll_Sites[c(259:260), "Annual_Preci"]<- 757

#It-Sro
dfAll_Sites[c(283, 285), "Annual_Preci"]<- 898

#JP-Tak
dfAll_Sites[c(291:296), "Annual_Preci"]<- 1024

#Ru-Fyo
dfAll_Sites[c(321:322), "Annual_Preci"]<- 671

#SE-Nor
dfAll_Sites[c(342,346), "Annual_Preci"]<- 561

#SE-Sk1
dfAll_Sites[349, "Annual_Preci"]<- 567

#SE-Sk2
dfAll_Sites[c(350:351), "Annual_Preci"]<- 573

#US-Lph
dfAll_Sites[417, "Annual_Preci"]<- 1071

#US-Syv
dfAll_Sites[473, "Annual_Preci"]<- 826

# Remove disturbed years
# sites: AU-Tum (2003: insect outbreak), BE-Bra (2000: thinning), CA-Ca2 (2001: forest not planted yet), 
#CA-SJ2 (2004: outlier for GPP values), NL-Loo (2006: storm), SE-Nor (2005: thinning occured this year), UK-Gri (2006:Thinning)
# US-Nc2 (2006: insect attack), US-SO2 (2004: wildfire), US-SP1 (2005: storm), US-Wcr (2002: insect attack)
dfAll_Sites<-dfAll_Sites[-c(3, 13, 37, 139, 313, 348, 362, 423, 438, 456, 482),]

# Remove years missing in each flux site
dfAll_Sites = dfAll_Sites[!is.na(dfAll_Sites$mean_Uncert),]

# Remove high gap filled fraction
dfAll_Sites<- dfAll_Sites[dfAll_Sites$mean_Uncert>0.80,]

# Remove uncomplete data collection within year. 
dfAll_Sites<- dfAll_Sites[!dfAll_Sites$GPP==0,]

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
Glob_Forest_df$GPP_ER<-Glob_Forest_df$GPP / Glob_Forest_df$Respiration
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
