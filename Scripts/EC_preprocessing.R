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
library(raster)
library(rgdal)
library(bfastSpatial)
library(SPEI)
library (soiltexture)
library(GSIF)

#1. Create a df for all fluxnet sites dataframe

#Open all files
dir <- file.path(path, 'Fluxnet_Data/Daily')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = fLoadFluxNCIntoDataframe(VarList.V.s=c('NEE_f','GPP_f','Reco',
                                                             'NEE_fqcOK'),
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
                                                     mean_Uncert=mean(NEE_fqcOK, na.rm=T)))

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

# Remove disturbed years
# sites: AU-Tum (2003: insect outbreak), BE-Bra (2000: thinning), CA-Ca2 (2001: forest not planted yet), 
#CA-SJ2 (2004: outlier for GPP values), JP-Tef (2002: year before plantation) NL-Loo (2006: storm), SE-Nor (2005: thinning occured this year), UK-Gri (2006:Thinning)
# US-Nc2 (2006: insect attack), US-SO2 (2004: wildfire), US-SP1 (2005: storm), US-Wcr (2002: insect attack)
dfAll_Sites<-dfAll_Sites[-c(3, 13, 64, 167, 371, 388, 423, 437, 508, 523, 541, 573),]

# Remove years missing in each flux site
dfAll_Sites = dfAll_Sites[!is.na(dfAll_Sites$mean_Uncert),]

# Remove high gap filled fraction
dfAll_Sites<- dfAll_Sites[dfAll_Sites$mean_Uncert>0.80,]

# Remove uncomplete data collection within year. 
dfAll_Sites<- dfAll_Sites[!dfAll_Sites$GPP==0,]

# Restructure dataframe
dfAll_Sites<-gather(dfAll_Sites, variable, value, -year, 
                    -Ecosystem, -Climate, -Disturbance,
                    -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement, -Study, -Lat, -Long )

#Reoder column
dfAll_Sites<- dfAll_Sites[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "value", 
                            "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")]

# Rename head column
colnames(dfAll_Sites)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values",
                         "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")


# Append forest global database
Glob_Forest_df<-read.csv("Input/Forest_Global_Database.csv", header = TRUE)
Glob_Forest_df$GPP_ER<-Glob_Forest_df$GPP / Glob_Forest_df$Respiration
Glob_Forest_df$NEP= -Glob_Forest_df$NEE
Glob_Forest_df$NEP_GPP<-Glob_Forest_df$NEP / Glob_Forest_df$GPP
Glob_Forest_df<- gather(Glob_Forest_df, variable, value, -year, 
                        -Ecosystem, -Climate, -Disturbance,
                        -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
                        -Study, -Lat, -Long, -x1, -y1, -x2, -y2, -x3, -y3, -x4, -y4 )
Glob_Forest_df<- Glob_Forest_df[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "value", 
                                  "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")]
colnames(Glob_Forest_df)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", 
                            "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")
dfAll_Sites<-rbind(dfAll_Sites, Glob_Forest_df)

# 4. Compute SPEI index from CRU dataset

# 4.1 Import precipitation and PET data from CRU

#PET
dir <- file.path(path, 'CRU/PET')
list <- list.files(dir, pattern=glob2rx('*.txt'), full.names=TRUE)
PET_Site <- list()
for(i in seq_along(list)) {
  PET_Site[[i]] = read.table(list[i], header=T,  check.names=FALSE, sep=",")
  PET_Site[[i]]<- PET_Site[[i]][-c(2),]
  paste(names(PET_Site[[i]]), PET_Site[[i]][1,])
  names(PET_Site[[i]]) <- paste(names(PET_Site[[i]]), PET_Site[[i]][1,])
  PET_Site[[i]] <- tail(PET_Site[[i]], -1)
  PET_Site[[i]]<- PET_Site[[i]] %>% gather(date, value, -`year 82`)
  PET_Site[[i]]<- PET_Site[[i]] %>% separate(date, c("year", "month"))
  colnames(PET_Site[[i]])<- c("ID", "year", "month", "PET")
  PET_Site[[i]]$ID <- as.character(PET_Site[[i]]$ID)
  PET_Site[[i]]$ID <- sapply(PET_Site[[i]]$ID, function(id) as.numeric(strsplit(id,"[.]")[[1]][2]))
  PET_Site[[i]] <- PET_Site[[i]][with(PET_Site[[i]], order(ID)),]
}

# Precipitation
dir <- file.path(path, 'CRU/Precipitation')
list <- list.files(dir, pattern=glob2rx('*.txt'), full.names=TRUE)
MAP_Site <- list()
for(i in seq_along(list)) {
  MAP_Site[[i]] = read.table(list[i], header=T,  check.names=FALSE, sep=",")
  MAP_Site[[i]]<- MAP_Site[[i]][-c(2),]
  paste(names(MAP_Site[[i]]), MAP_Site[[i]][1,])
  names(MAP_Site[[i]]) <- paste(names(MAP_Site[[i]]), MAP_Site[[i]][1,])
  MAP_Site[[i]] <- tail(MAP_Site[[i]], -1)
  MAP_Site[[i]]<- MAP_Site[[i]] %>% gather(date, value, -`year 82`)
  MAP_Site[[i]]<- MAP_Site[[i]] %>% separate(date, c("year", "month"))
  colnames(MAP_Site[[i]])<- c("ID", "year", "month", "MAP")
  MAP_Site[[i]]$ID <- as.character(MAP_Site[[i]]$ID)
  MAP_Site[[i]]$ID <- sapply(MAP_Site[[i]]$ID, function(id) as.numeric(strsplit(id,"[.]")[[1]][2]))
  MAP_Site[[i]] <- MAP_Site[[i]][with(MAP_Site[[i]], order(ID)),]
}

# 4.2 Compute SPEI from CRU

# Compute SPI
MAP_Site_df<- do.call("rbind", MAP_Site)
colnames(MAP_Site_df)<- c("ID", "year", "month", "MAP")
PET_Site_df<- do.call("rbind", PET_Site)
SPEI_Site_df<- cbind(MAP_Site_df, PET_Site_df$PET)
SPEI_Site_df[SPEI_Site_df == -9999.0 ] <- NA # replace missing values with NA
for(id in unique(SPEI_Site_df$ID)){
  SPEI_Site<- spei(SPEI_Site_df$MAP - SPEI_Site_df$PET, 9, na.rm=T)
  SPEI_Site_df$SPEI_CRU<- SPEI_Site$fitted
}

# Average SPEI per year
SPEI_Site_df<- ddply(SPEI_Site_df, .(ID, year),
                     summarize,
                     SPEI_CRU= mean(SPEI_CRU, na.rm=T))

# 4.3 Add Site ID to dataset
dir <- file.path(path, 'GPP_Extraction')
list <- list.files(dir, pattern=glob2rx('*.csv'), full.names=TRUE)
Site_ID<- read.csv(list, header=T,  check.names=FALSE, sep=",")
SPEI_Site_df<- merge(SPEI_Site_df, Site_ID, by.x="ID", by.y="ID", all.x=TRUE)
SPEI_Site_df<- SPEI_Site_df[c("Site_ID", "year", "SPEI_CRU")]

# 4.4 Include SPI from CRU into flux dataset
dfAll_Sites<- merge(dfAll_Sites, SPEI_Site_df, by=c("Site_ID", "year"), all.x=TRUE)

# 5. Compute SPI index from CRU dataset

# 5.1 Import precipitation data from CRU
dir <- file.path(path, 'CRU/Precipitation')
list <- list.files(dir, pattern=glob2rx('*.txt'), full.names=TRUE)
MAP_Site <- list()
for(i in seq_along(list)) {
  MAP_Site[[i]] = read.table(list[i], header=T,  check.names=FALSE, sep=",")
  MAP_Site[[i]]<- MAP_Site[[i]][-c(2),]
  paste(names(MAP_Site[[i]]), MAP_Site[[i]][1,])
  names(MAP_Site[[i]]) <- paste(names(MAP_Site[[i]]), MAP_Site[[i]][1,])
  MAP_Site[[i]] <- tail(MAP_Site[[i]], -1)
  MAP_Site[[i]]<- MAP_Site[[i]] %>% gather(date, value, -`year 82`)
  MAP_Site[[i]]<- MAP_Site[[i]] %>% separate(date, c("year", "month"))
  colnames(MAP_Site[[i]])<- c("ID", "year", "month", "MAP")
  MAP_Site[[i]]$ID <- as.character(MAP_Site[[i]]$ID)
  MAP_Site[[i]]$ID <- sapply(MAP_Site[[i]]$ID, function(id) as.numeric(strsplit(id,"[.]")[[1]][2]))
  MAP_Site[[i]] <- MAP_Site[[i]][with(MAP_Site[[i]], order(ID)),]
}

# 5.2 Compute SPI from CRU

# Compute SPI
SPI_Site_df<- do.call("rbind", MAP_Site)
colnames(SPI_Site_df)<- c("ID", "year", "month", "MAP")
for(id in unique(SPI_Site_df$ID)){
  SPI_Site<- spi(SPI_Site_df$MAP, 9, na.rm=T)
  SPI_Site_df$SPI_CRU<- SPI_Site$fitted
}

# Average SPI per year
SPI_Site_df<- ddply(SPI_Site_df, .(ID, year),
                    summarize,
                    SPI_CRU= mean(SPI_CRU, na.rm=T))

# 5.3 Add Site ID to dataset
SPI_Site_df<- merge(SPI_Site_df, Site_ID, by.x="ID", by.y="ID", all.x=TRUE)
SPI_Site_df<- SPI_Site_df[c("Site_ID", "year", "SPI_CRU")]

# 5.4 Include SPI from CRU into flux dataset
dfAll_Sites<- merge(dfAll_Sites, SPI_Site_df, by=c("Site_ID", "year"), all.x=TRUE)

# 6. Compute temperature anomalies from CRU data

# 6.1 Import precipitation data from CRU
dir <- file.path(path, 'CRU/Mean_Temp')
list <- list.files(dir, pattern=glob2rx('*.txt'), full.names=TRUE)
MAT_Site <- list()
for(i in seq_along(list)) {
  MAT_Site[[i]] = read.table(list[i], header=T,  check.names=FALSE, sep=",")
  MAT_Site[[i]]<- MAT_Site[[i]][-c(2),]
  paste(names(MAT_Site[[i]]), MAT_Site[[i]][1,])
  names(MAT_Site[[i]]) <- paste(names(MAT_Site[[i]]), MAT_Site[[i]][1,])
  MAT_Site[[i]] <- tail(MAT_Site[[i]], -1)
  MAT_Site[[i]]<- MAT_Site[[i]] %>% gather(date, value, -`year 82`)
  MAT_Site[[i]]<- MAT_Site[[i]] %>% separate(date, c("year", "month"))
  colnames(MAT_Site[[i]])<- c("ID", "year", "month", "MAT")
  MAT_Site[[i]]$ID <- as.character(MAT_Site[[i]]$ID)
  MAT_Site[[i]]$ID <- sapply(MAT_Site[[i]]$ID, function(id) as.numeric(strsplit(id,"[.]")[[1]][2]))
  MAT_Site[[i]] <- MAT_Site[[i]][with(MAT_Site[[i]], order(ID)),]
}

# Compute average temperature per site since 1900
Mean_MAT_Site<- ddply(do.call("rbind", MAT_Site), .(ID, year),
                      summarize,
                      MAT_Mean= mean(MAT, na.rm=T))

# Add Site ID
Mean_MAT_Site$Site_ID<- Site_ID$Site_ID

# Add average temperature per site since 1900 in flux dataframe
dfAll_Sites<- merge(dfAll_Sites, Mean_MAT_Site, by=c("Site_ID"), all.x=TRUE)

# Compute average temperature per site per year
MAT_CRU_Site<- ddply(do.call("rbind", MAT_Site), .(ID, year),
                     summarize,
                     MAT_CRU= mean(MAT, na.rm=T))

MAT_CRU_Site<- merge(MAT_CRU_Site, Site_ID, by.x="ID", by.y="ID", all.x=TRUE)
MAT_CRU_Site<- MAT_CRU_Site[c("Site_ID", "year", "MAT_CRU")]

# Add average temperature per site per yeat in flux dataframe
dfAll_Sites<- dfAll_Sites[c("Site_ID", "year.x", "Type_Flux", "values", "Stand_Age", "Climate", "Ecosystem", 
                            "Study" , "Lat", "Long", "SPEI_CRU","SPI_CRU","MAT_Mean")]
colnames(dfAll_Sites)<- c("Site_ID", "year", "Type_Flux", "values", "Stand_Age", "Climate", "Ecosystem", 
                          "Study" , "Lat", "Long", "SPEI_CRU","SPI_CRU","MAT_Mean") 
dfAll_Sites<- merge(dfAll_Sites, MAT_CRU_Site, by=c("Site_ID", "year"), all.x=TRUE)

# Compute temperature anomalies
dfAll_Sites$MAT_An<- dfAll_Sites$MAT_CRU - dfAll_Sites$MAT_Mean

# 6. Append CRU precipitation  to dataset

# Compute sum precipitation per site per year
MAP_CRU_Site<- ddply(do.call("rbind", MAP_Site), .(ID, year),
                     summarize,
                     MAP_CRU= sum(MAP, na.rm=T))

MAP_CRU_Site<- merge(MAP_CRU_Site, Site_ID, by.x="ID", by.y="ID", all.x=TRUE)
MAP_CRU_Site<- MAP_CRU_Site[c("Site_ID", "year", "MAP_CRU")]

# Add average precipitation per site per yeat in flux dataframe
dfAll_Sites<- merge(dfAll_Sites, MAP_CRU_Site, by=c("Site_ID", "year"), all.x=TRUE)

# 7. Compute CEC

# 7.1 Import Site coordinate
NEP<- readRDS("Output/NEP.rds")
Site_list<-ddply(NEP, .(Site_ID),
                 summarise,
                 y= mean(Lat, na.rm=T),
                 x=mean(Long, na.rm=T))
Site.xy<- Site_list[c("x", "y")]

#7.2 SoilGrid 1km raster - CEC
dir <- file.path(path, 'Soil_Data/CEC/1km')
file <- list.files(dir, pattern=glob2rx('*.tif'), full.names=TRUE)
Soil_1km_CEC<- stack(file)

# Extract CEC value per site
Soil_1km_Data_Site <- data.frame(coordinates(Site.xy),
                                 Site_list$Site_ID, 
                                 extract(Soil_1km_CEC, Site.xy))

# Compute CEC total
wt<-c(0.05, 0.10, 0.15, 0.3, 0.4, 1)
Soil_1km_Data_Site$CEC_Total_1km<-apply(Soil_1km_Data_Site[,-c(1:3)],1,weighted.mean,w=wt)
Soil_1km_Data_Site<- Soil_1km_Data_Site[c("Site_list.Site_ID", "CEC_Total_1km")]
colnames(Soil_1km_Data_Site)<- c("Site_ID", "CEC_Total_1km")

# Add CEC per site in flux dataframe
dfAll_Sites<- merge(dfAll_Sites, Soil_1km_Data_Site, by=c("Site_ID"), all.x=TRUE)

# 8. Compute clay content

#8.1 SoilGrid 1km raster- Clay content
dir <- file.path(path, 'Soil_Data/Clay_Content/1km')
file <- list.files(dir, pattern=glob2rx('*.tif'), full.names=TRUE)
Soil_1km_Clay<- stack(file)

# Extract clay value per site

# Extract data and create dataframe
Soil_1km_Data_Site <- data.frame(coordinates(Site.xy),
                                 Site_list$Site_ID, 
                                 extract(Soil_1km_Clay, Site.xy))

# Compute clay total
wt<-c(0.05, 0.10, 0.15, 0.3, 0.4, 1)
Soil_1km_Data_Site$Clay_1km<-apply(Soil_1km_Data_Site[,-c(1:3)],1,weighted.mean,w=wt)
Soil_1km_Data_Site<- Soil_1km_Data_Site[c("Site_list.Site_ID", "Clay_1km")]
colnames(Soil_1km_Data_Site)<- c("Site_ID", "Clay_1km")

# Add clay per site in flux dataframe
dfAll_Sites<- merge(dfAll_Sites, Soil_1km_Data_Site, by=c("Site_ID"), all.x=TRUE)

# 9. Compute Silt content

#9.1 SoilGrid 1km raster- Silt content
dir <- file.path(path, 'Soil_Data/Silt_Content/1km')
file <- list.files(dir, pattern=glob2rx('*.tif'), full.names=TRUE)
Soil_1km_Silt<- stack(file)

# Extract Silt value per site

# Extract data and create dataframe
Soil_1km_Data_Site <- data.frame(coordinates(Site.xy),
                                 Site_list$Site_ID, 
                                 extract(Soil_1km_Silt, Site.xy))

# Compute Silt total
wt<-c(0.05, 0.10, 0.15, 0.3, 0.4, 1)
Soil_1km_Data_Site$Silt_1km<-apply(Soil_1km_Data_Site[,-c(1:3)],1,weighted.mean,w=wt)
Soil_1km_Data_Site<- Soil_1km_Data_Site[c("Site_list.Site_ID", "Silt_1km")]
colnames(Soil_1km_Data_Site)<- c("Site_ID", "Silt_1km")

# Add Silt per site in flux dataframe
dfAll_Sites<- merge(dfAll_Sites, Soil_1km_Data_Site, by=c("Site_ID"), all.x=TRUE)

# 10. Compute Sand content

#10.1 SoilGrid 1km raster- Sand content
dir <- file.path(path, 'Soil_Data/Sand_Content/1km')
file <- list.files(dir, pattern=glob2rx('*.tif'), full.names=TRUE)
Soil_1km_Sand<- stack(file)

# Extract Sand value per site

# Extract data and create dataframe
Soil_1km_Data_Site <- data.frame(coordinates(Site.xy),
                                 Site_list$Site_ID, 
                                 extract(Soil_1km_Sand, Site.xy))

# Compute Sand total
wt<-c(0.05, 0.10, 0.15, 0.3, 0.4, 1)
Soil_1km_Data_Site$Sand_1km<-apply(Soil_1km_Data_Site[,-c(1:3)],1,weighted.mean,w=wt)
Soil_1km_Data_Site<- Soil_1km_Data_Site[c("Site_list.Site_ID", "Sand_1km")]
colnames(Soil_1km_Data_Site)<- c("Site_ID", "Sand_1km")

# Add Sand per site in flux dataframe
dfAll_Sites<- merge(dfAll_Sites, Soil_1km_Data_Site, by=c("Site_ID"), all.x=TRUE)

#11. Perform soil texture
tdf <- NEP_Mean_Site[,c("Clay_1km", "Silt_1km", "Sand_1km")]
tdf<- tdf[77,]
tdf$Sum = rowSums(tdf)
for(i in c("Clay_1km", "Silt_1km", "Sand_1km")) { tdf[,i] <- tdf[,i]/tdf$Sum * 100 }
names(tdf)[1:3] <- c("CLAY", "SILT", "SAND")
TT.plot(class.sys = "USDA.TT", tri.data = tdf, grid.show = FALSE, pch="+", cex=1, col="red")

# 11. Create dataframe for post-processing

# 11.1 Append soil texture to dataframe
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),] #keep the site used in the study
dir <- file.path(path, 'Soil_Data/Soil_Texture')
file <- list.files(dir, pattern=glob2rx('*.csv'), full.names=TRUE)
Soil_text<- read.csv(file)
Flux_High<- dfAll_Sites<- merge(Flux_High, Soil_text, by=c("Site_ID"), all.x=TRUE)
Flux_High$Clay_Silt<- Flux_High$Clay_1km + Flux_High$Silt_1km 

# Subset dataframe according to type of flux
NEP<- Flux_High[Flux_High$Type_Flux %in% c("NEP"),]
GPP<- Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco<- Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
Ratio_NEP_GPP<- Flux_High[Flux_High$Type_Flux %in% c("NEP_GPP"),]
Ratio_GPP_Reco<- Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]

#11.2 Compute GPPclimax based on the lieth model
params <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$MAT_CRU))))
GPP$GPPp<-with(as.list(params), GPP1000*((1-exp(-k*GPP$MAP_CRU))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))

# 8.4 Create average site dataframe
NEP_Mean_Site<-ddply(NEP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     MAP_CRU=mean(MAP_CRU, na.rm=T),
                     MAT_CRU=mean(MAT_CRU, na.rm=T), 
                     SPI_CRU=mean(SPI_CRU, na.rm=T),
                     SPEI_CRU=mean(SPEI_CRU, na.rm=T),
                     MAT_An=mean(MAT_An, na.rm=T),
                     CEC_Total_1km= mean(CEC_Total_1km, na.rm=T),
                     Clay_Silt= mean(Clay_Silt, na.rm=T),
                     Soil_Quality= mean(Soil_Quality, na.rm=T))

GPP_Mean_Site<-ddply(GPP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     MAP_CRU=mean(MAP_CRU, na.rm=T),
                     MAT_CRU=mean(MAT_CRU, na.rm=T), 
                     SPI_CRU=mean(SPI_CRU, na.rm=T),
                     SPEI_CRU=mean(SPEI_CRU, na.rm=T),
                     MAT_An=mean(MAT_An, na.rm=T),
                     CEC_Total_1km= mean(CEC_Total_1km, na.rm=T),
                     Clay_Silt= mean(Clay_Silt, na.rm=T),
                     Soil_Quality= mean(Soil_Quality, na.rm=T))

Ratio_NEP_GPP_Mean_Site<-ddply(Ratio_NEP_GPP, .(Site_ID, Type_Flux),
                               summarise,
                               Stand_Age= mean(Stand_Age, na.rm=T),
                               values=mean(values, na.rm=T),
                               MAP_CRU=mean(MAP_CRU, na.rm=T),
                               MAT_CRU=mean(MAT_CRU, na.rm=T), 
                               SPI_CRU=mean(SPI_CRU, na.rm=T),
                               SPEI_CRU=mean(SPEI_CRU, na.rm=T),
                               MAT_An=mean(MAT_An, na.rm=T),
                               CEC_Total_1km= mean(CEC_Total_1km, na.rm=T),
                               Clay_Silt= mean(Clay_Silt, na.rm=T),
                               Soil_Quality= mean(Soil_Quality, na.rm=T))

# 8.5 Save all dataframe in a rds files
saveRDS(dfAll_Sites, file="Output/df_Annual_Flux.rds")
saveRDS(NEP, file="Output/NEP.rds")
saveRDS(NEP_Mean_Site, file="Output/NEP_Mean_Site.rds")
saveRDS(GPP, file="Output/GPP.rds")
saveRDS(GPP_Mean_Site, file="Output/GPP_Mean_Site.rds")
saveRDS(Ratio_NEP_GPP, file="Output/Ratio_NEP_GPP.rds")
saveRDS(Ratio_NEP_GPP_Mean_Site, file="Output/Ratio_NEP_GPP_Mean_Site.rds")
