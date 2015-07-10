## Script to analyse IAV of annual carbon flux
## Author: Simon Besnard
## 9.07.2015
###################################
## Load the necessary packages
# install.packages('ggplot2')
# install.packages("scales") 
# install.packages("lubridate")
# install.packages("gridExtra")
# install.packages("REddyProc", repos="http://R-Forge.R-project.org", type="source")
# install.packages('RNetCDF')
# install.packages('minpack.lm')
# install.packages("dplyr")
# install.packages("testthat")
# install.packages('gtools')
# install.packages("FactoMineR")
# install.packages("leaps")
# install.packages("gamm4")
# install.packages("tidyr")
# install.packages("manipulate")
# install.packages("fitdistrplus")
# install.packages("evd")
# install.packages("flexmix")
# install.packages("mosaic")
# install.packages("bootstrap")
# install.packages("xts")

library(RNetCDF)
library (ggplot2)
library(scales)
library(lubridate)
library(gridExtra)
library(REddyProc)
library (dplyr)
library (plyr)
library(tidyr)
library(manipulate)
library(MASS)
library(mosaic)
library(xts)

#1. Create a df for all fluxnet sites dataframe

#Open all files
dir <- file.path(path, 'Fluxnet_Data/Data_4disturbed')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = fLoadFluxNCIntoDataframe(VarList.V.s=c('NEE_f','GPP_f','Reco','NEE_fqcOK'),
                                               FileName.s=list[i], NcPackage.s = 'RNetCDF')
}

# Remove data before disturbance from the list of dataframe
Site_Date<-read.csv("Input/Potential_Sites.csv", header = TRUE)
Site_Date$Measure_Date<- as.Date(Site_Date$Measure_Date)
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]$DateTime=as.Date(Fluxnet_Site[[i]]$DateTime)
}

for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]<-with(Fluxnet_Site[[i]], Fluxnet_Site[[i]][(Fluxnet_Site[[i]]$DateTime > Site_Date$Measure_Date[i]),])
}

#2.Analyse IAV of carbon fluxes since disturbances 

# 2.1. Aggregate flux data per week

# Add week number to dataframe
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]$Week=floor_date(Fluxnet_Site[[i]]$DateTime, "week")
}

# compute mean values per week
Flux_Week<-list()
for (i in seq_along(Fluxnet_Site)){
  Flux_Week[[i]]<-ddply(Fluxnet_Site[[i]], "Week", summarise, 
         GPP_Week = mean(GPP_f),
         Reco_Week= mean(Reco),
         NEE_Week=mean(NEE_f))
}

# Compute sum/mean annual +/- sd
Mean_Annual_Flux<-lapply(Fluxnet_Site, function (x) ddply(x, "year",
                                                     summarize,
                                                     Annual_NEE= mean(NEE_f, na.rm=T),
                                                     Annual_GPP= mean(GPP_f, na.rm=T),
                                                     Annual_Reco= mean(Reco,na.rm=T)))

#Add ancillary data to the list of dataframe
Site_Date<-read.csv("Input/Potential_Sites.csv", header = TRUE)
Site_Date$Measure_Date<-as.Date(Site_Date$Measure_Date)

for (i in seq_along(Flux_Week)){
  Flux_Week[[i]]$Site_ID<- Site_Date$ID[i]
  Flux_Week[[i]]$year<- year(Flux_Week[[i]]$Week)
  Flux_Week[[i]]$Ecosytem<- Site_Date$Ecosystem[i]
  Flux_Week[[i]]$Climate<- Site_Date$Climate[i]
  Flux_Week[[i]]$Disturbance<- Site_Date$Type_Disturbance[i]
  Flux_Week[[i]]$Species<- Site_Date$Species[i]
  Flux_Week[[i]]$Day_Disturbance<- as.numeric(Flux_Week[[i]]$Week- Site_Date$Measure_Date[i])
  Flux_Week[[i]]$Year_Disturbance<- Flux_Week[[i]]$Day_Disturbance/365
  Flux_Week[[i]]<-merge(Flux_Week[[i]], Mean_Annual_Flux[[i]], by="year")
  Flux_Week[[i]]$Norm_GPP<-Flux_Week[[i]]$GPP_Week - Flux_Week[[i]]$Annual_GPP
  Flux_Week[[i]]$Norm_NEE<-Flux_Week[[i]]$NEE_Week - Flux_Week[[i]]$Annual_NEE
  Flux_Week[[i]]$Norm_Reco<-Flux_Week[[i]]$Reco_Week - Flux_Week[[i]]$Annual_Reco
}

# combine fluxnet sites in one dataframe
dfAll_Sites<- do.call("rbind", Flux_Week)

# Restructure dataframe
dfAll_Sites<-gather(dfAll_Sites, Type_Flux, values, -year, -Species, -Ecosytem, -Climate, -Disturbance, 
                    -Year_Disturbance, -Site_ID, -Day_Disturbance,-Week)

#.2.2.Decompose seasonality
ts<-xts(dfAll_Sites$values, order.by=as.POSIXct(dfAll_Sites$Year_Disturbance))

# 2.3. Plot IAV since disturbance

# Clean Dataframe
dfAll_Sites<- dfAll_Sites[!dfAll_Sites$values>1e+34,] #Remove outliers
dfAll_Sites<- dfAll_Sites[!dfAll_Sites$values<(-3e10),] #Remove outliers
dfAll_Sites<- dfAll_Sites[dfAll_Sites$Type_Flux %in% c("Norm_GPP", "Norm_NEE", "Norm_Reco"),]
IAV_Harvest<- dfAll_Sites[dfAll_Sites$Disturbance %in% c("Harvest"),]
IAV_Fire<- dfAll_Sites[dfAll_Sites$Disturbance %in% c("Fire"),]

#Plot
gg1<-ggplot(transform(IAV_Harvest,
                              Type_Flux=factor(Type_Flux,levels=c("Norm_NEE","Norm_GPP","Norm_Reco"))))+
    geom_point(aes(Year_Disturbance, values), size=1, shape=3) +
  facet_grid(Type_Flux~Disturbance, scales = "free")+
  stat_function(fun=function(x) sin(2*pi*x) + cos(2*pi*x), colour="red")+ 
  xlab("Day since Disturbance") + ylab(" Carbon flux (g.m-2.d-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
