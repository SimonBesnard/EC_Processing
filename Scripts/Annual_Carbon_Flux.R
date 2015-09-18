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
dir <- file.path(path, 'Fluxnet_Data')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = fLoadFluxNCIntoDataframe(VarList.V.s=c('NEE_f','GPP_f','Reco','NEE_fqcOK', 'precip_f', 'Tair_f', 'Tsoil_f','Rg_f', 'SWC1_f', 'SWC2_f'),
                                          FileName.s=list[i], NcPackage.s = 'RNetCDF')
}

#2. Compute dataframe for time since disturbance or stand age analysis
#Open csv file with flux site metadata
Site_Date<-read.csv("Input/Potential_Sites.csv", header = TRUE)

# Convert Datetime into a Date object
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]$DateTime=as.Date(Fluxnet_Site[[i]]$DateTime)
}

# 2.1. Option 1: Remove data before disturbance from the list of dataframe
# Convert date of disturbance into a date object
#Site_Date$Date_Disturbance<- as.Date(Site_Date$Date_Disturbance)

#Remove data before disturbance
# for (i in seq_along(Fluxnet_Site)){
# Fluxnet_Site[[i]]<-with(Fluxnet_Site[[i]], Fluxnet_Site[[i]][(Fluxnet_Site[[i]]$DateTime > Site_Date$Date_Disturbance[i]),])
# }

# 2.2. Option 2: Remove data before plantation from the list of dataframe
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
                    mean_NEE= mean(NEE_f, na.rm=T),
                    mean_GPP= mean(GPP_f, na.rm=T),
                    mean_TER= mean(Reco,na.rm=T),
                    mean_NEP= mean(-NEE_f, na.rm=T),
                    sd_NEE= sd(NEE_f, na.rm=T),
                    sd_GPP= sd(GPP_f, na.rm=T),
                    sd_Reco= sd(Reco, na.rm=T),
                    sd_NEP= sd(-NEE_f, na.rm=T),
                    mean_Uncert=mean(NEE_fqcOK, na.rm=T), 
                    Annual_Preci= sum(precip_f, na.rm=T),
                    Tair=mean(Tair_f, na.rm=T),
                    Tsoil=mean(Tsoil_f, na.rm=T),
                    SWC1=mean(SWC1_f, na.rm=T),
                    SWC2=mean(SWC2_f, na.rm=T),
                    Rg=mean(Rg_f, na.rm=T)))

#Add type of ecosystem, type of climate, type of disturbance and year after disturbance for each sites
for (i in seq_along(Sum_Sd_Flux)){
  Sum_Sd_Flux[[i]]$Ecosystem<- Site_Date$Ecosystem[i]
  Sum_Sd_Flux[[i]]$Climate<- Site_Date$Climate[i]
  Sum_Sd_Flux[[i]]$Disturbance<- Site_Date$Type_Disturbance[i]
  Sum_Sd_Flux[[i]]$Stand_Age<- Sum_Sd_Flux[[i]]$year- year(as.Date(Site_Date$Plantation_Date[i]))
  Sum_Sd_Flux[[i]]$GPP_ER<- Sum_Sd_Flux[[i]]$GPP / Sum_Sd_Flux[[i]]$Respiration
  Sum_Sd_Flux[[i]]$Site_ID<- Site_Date$ID[i]
  Sum_Sd_Flux[[i]]$Stand_Replacement<- Site_Date$Stand_Replacement[i]
  Sum_Sd_Flux[[i]]$Int_Replacement<- Site_Date$Intensity_Replacement[i]
  Sum_Sd_Flux[[i]]$SWC<-(Sum_Sd_Flux[[i]]$SWC1 + Sum_Sd_Flux[[i]]$SWC2)/2 
}

#Remove outlier sites (US-Blo, CZ-Bk1, US-DK3, IT_Ro1, IT-Ro2, US_Bn2, US-NC2)
Sum_Sd_Flux<- Sum_Sd_Flux[-c(28,36,37,49,50,51,53,63)]

#Combine the flux sites in one dataframe
dfAll_Sites<- do.call("rbind", Sum_Sd_Flux)
dfAll_Sites<-dfAll_Sites[-c(127),]# the measurements seems to be an outlier

# Remove years missing in each flux site
dfAll_Sites = dfAll_Sites[!is.na(dfAll_Sites$mean_NEE),]

# Remove high gap filled fraction
dfAll_Sites<- dfAll_Sites[dfAll_Sites$mean_Uncert>0.85,]

# Restructure dataframe
dfAll_Sites<-gather(dfAll_Sites, variable, value, -Annual_Preci, -year, 
                    -Ecosystem, -mean_Uncert, -Climate, -Disturbance,
                    -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
                    -Tair, -Tsoil, -SWC, -Rg)

#Reoder column
dfAll_Sites<- dfAll_Sites[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "value", "Annual_Preci", 
                            "mean_Uncert", "Tair", "Tsoil","SWC","Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem")]

# Rename head column
colnames(dfAll_Sites)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", 
                            "mean_Uncert", "Tair", "Tsoil","SWC","Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem")

# Reclassify climate classification
dfAll_Sites$Climate<-ifelse((dfAll_Sites$Climate =="Af" | dfAll_Sites$Climate =="Am"), "Tropical",
                    ifelse((dfAll_Sites$Climate =="Cfa" | dfAll_Sites$Climate =="Cfb" | dfAll_Sites$Climate =="Cfc" | dfAll_Sites$Climate =="Csa" | dfAll_Sites$Climate =="Csb"), "Temperate",
                           ifelse((dfAll_Sites$Climate =="Dfb" | dfAll_Sites$Climate =="Dfc" | dfAll_Sites$Climate =="Dsc"), "Continental", NA)))

# Save dataframe in a rds file
saveRDS(dfAll_Sites, file="Output/df_Annual_Flux.rds")