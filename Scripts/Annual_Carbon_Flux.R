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

#1. Create a df for all fluxnet sites dataframe

#Open all files
dir <- file.path(path, 'Fluxnet_Data/Daily')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = fLoadFluxNCIntoDataframe(VarList.V.s=c('NEE_f','GPP_f','Reco',
                                                             'NEE_fqcOK', 'precip_f', 'Tair_f', 'Tsoil_f','Rg_f', 'Rn_f', 'LE_f'),
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
                    Tsoil=mean(Tsoil_f, na.rm=T),
                    LE=sum(LE_f, na.rm=T),
                    Rn=sum(Rn_f, na.rm=T),
                    Rg=mean(Rg_f, na.rm=T)))

#Add type of ecosystem, type of climate, type of disturbance and year after disturbance for each sites
for (i in seq_along(Sum_Sd_Flux)){
  Sum_Sd_Flux[[i]]$Ecosystem<- Site_Date$Ecosystem[i]
  Sum_Sd_Flux[[i]]$Climate<- Site_Date$Climate[i]
  Sum_Sd_Flux[[i]]$Disturbance<- Site_Date$Type_Disturbance[i]
  Sum_Sd_Flux[[i]]$Stand_Age<- Sum_Sd_Flux[[i]]$year- year(as.Date(Site_Date$Plantation_Date[i]))
  Sum_Sd_Flux[[i]]$GPP_ER<- Sum_Sd_Flux[[i]]$GPP / Sum_Sd_Flux[[i]]$Respiration
  Sum_Sd_Flux[[i]]$NEP_GPP<- Sum_Sd_Flux[[i]]$NEP / Sum_Sd_Flux[[i]]$GPP
  Sum_Sd_Flux[[i]]$GPP_NEP<- Sum_Sd_Flux[[i]]$GPP / Sum_Sd_Flux[[i]]$NEP
  Sum_Sd_Flux[[i]]$Site_ID<- Site_Date$ID[i]
  Sum_Sd_Flux[[i]]$Study<- Site_Date$Study[i]
  Sum_Sd_Flux[[i]]$Stand_Replacement<- Site_Date$Stand_Replacement[i]
  Sum_Sd_Flux[[i]]$Int_Replacement<- Site_Date$Intensity_Replacement[i]
  Sum_Sd_Flux[[i]]$ET<-Sum_Sd_Flux[[i]]$LE/Sum_Sd_Flux[[i]]$Rn
  Sum_Sd_Flux[[i]]$Lat<- Site_Date$x[i]
  Sum_Sd_Flux[[i]]$Long<- Site_Date$y[i]
  
}

#Combine the flux sites in one dataframe
dfAll_Sites<- do.call("rbind", Sum_Sd_Flux)

# Remove outliers values
# sites: BE-Bra (2000: thinning), CA-Ca1 (1997: disturbance occured), CA-Ca2 (2000/2001: forest not planted yet) CA-NS4 (2002: disturbance occured), CA-NS6 (2001), 
#CA-OJP (1999: disturbance occured), CA-SJ2 (2003\2004: outlier for GPP values), NL-Loo (2006: storm), SE-Nor (2005: thinning occured this year), US-Nc2 (2006: insect attack), US-SP1 (2005: storm), 
#US-SP3 (2000: storm), US-Wcr (2002: insect attack)
dfAll_Sites<-dfAll_Sites[-c(10, 18, 27, 28, 65, 73, 98, 129, 130, 236, 254, 318, 325, 333, 346),]

# Remove years missing in each flux site
dfAll_Sites = dfAll_Sites[!is.na(dfAll_Sites$mean_Uncert),]

# Remove high gap filled fraction
dfAll_Sites<- dfAll_Sites[dfAll_Sites$mean_Uncert>0.80,]

# Restructure dataframe
dfAll_Sites<-gather(dfAll_Sites, variable, value, -Annual_Preci, -year, 
                    -Ecosystem, -mean_Uncert, -Climate, -Disturbance,
                    -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
                    -Tair, -Tsoil, -Rg, -Rn, -LE, -ET, -Study, -Lat, -Long )

#Reoder column
dfAll_Sites<- dfAll_Sites[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "value", "Annual_Preci", 
                            "mean_Uncert", "Tair", "Tsoil","Rg","Rn", "LE", "ET", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")]

# Rename head column
colnames(dfAll_Sites)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", 
                            "mean_Uncert", "Tair", "Tsoil","Rg", "Rn", "LE", "ET", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")

# Save dataframe in a rds file
saveRDS(dfAll_Sites, file="Output/df_Annual_Flux.rds")
