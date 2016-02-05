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

# Remove disturbed years
# sites: AU-Tum (2003: insect outbreak), BE-Bra (2000: thinning), CA-Ca2 (2001: forest not planted yet), 
#CA-SJ2 (2004: outlier for GPP values), NL-Loo (2006: storm), SE-Nor (2005: thinning occured this year), UK-Gri (2006:Thinning)
# US-Nc2 (2006: insect attack), US-SO2 (2004: wildfire), US-SP1 (2005: storm), US-Wcr (2002: insect attack)
dfAll_Sites<-dfAll_Sites[-c(3, 13, 48, 150, 358, 393, 407, 470, 485, 503, 529),]

# Remove years missing in each flux site
dfAll_Sites = dfAll_Sites[!is.na(dfAll_Sites$mean_Uncert),]

# Remove high gap filled fraction
dfAll_Sites<- dfAll_Sites[dfAll_Sites$mean_Uncert>0.80,]

# Remove uncomplete data collection within year. 
dfAll_Sites<- dfAll_Sites[!dfAll_Sites$GPP==0,]

# Add climate data to the dataframe
rownames(dfAll_Sites) <- seq(length=nrow(dfAll_Sites)) 

#CA-NS1
dfAll_Sites[c(40), "Annual_Preci"]<- 500

#CA-NS3
dfAll_Sites[44, "Annual_Preci"]<- 502

#CZ-Bk1
dfAll_Sites[c(89), "Annual_Preci"]<- 1026

#FI-Sod
dfAll_Sites[c(143), "Annual_Preci"]<- 525

#Il-Yat
dfAll_Sites[c(173,174), "Annual_Preci"]<- 277

#It-Col
dfAll_Sites[c(177), "Annual_Preci"]<- 971

#It-Lav
dfAll_Sites[c(184:185), "Annual_Preci"]<- 757

#It-Sro
dfAll_Sites[c(201, 203), "Annual_Preci"]<- 898

#JP-Tak
dfAll_Sites[c(208:213), "Annual_Preci"]<- 1024

#Ru-Fyo
dfAll_Sites[c(227:228), "Annual_Preci"]<- 671

#SE-Nor
dfAll_Sites[c(240:241), "Annual_Preci"]<- 561

#SE-Sk1
dfAll_Sites[242, "Annual_Preci"]<- 567

#SE-Sk2
dfAll_Sites[c(243:244), "Annual_Preci"]<- 573

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
Glob_Forest_df<- gather(Glob_Forest_df, variable, value, -Annual_Preci, -year, 
                    -Ecosystem, -Climate, -Disturbance,
                    -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
                    -Tair, -Rg, -Study, -Lat, -Long, -x1, -y1, -x2, -y2, -x3, -y3, -x4, -y4 )
Glob_Forest_df<- Glob_Forest_df[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "value", "Annual_Preci", 
                            "Tair", "Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")]
colnames(Glob_Forest_df)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", 
                         "Tair", "Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")
dfAll_Sites<-rbind(dfAll_Sites, Glob_Forest_df)

#4.  Get SPEI index from RS product and inlcude in the dataset

# 4.1 Import SPEI dataset
dir <- file.path(path, 'SPEI')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)
SPEI_nc<- read.nc(open.nc(list[1]))
SPEI_b <- brick(list, varname="spei")

# 4.2. Import Site coordinates
NEP<- readRDS("Output/NEP.rds")
Site_list<-ddply(NEP, .(Site_ID),
                 summarise,
                 Lat= mean(Lat, na.rm=T),
                 Long=mean(Long, na.rm=T))

# 4.3 Extract SPEI per site
point <- data.frame(cbind(Site_list$Long, Site_list$Lat))
SPEI_Site_ts<- zooExtract(SPEI_b, point)
SPEI_Site_ts <- aggregate(SPEI_Site_ts, year, mean, na.rm=T)
SPEI_Site_df<- data.frame(SPEI_Site_ts)
SPEI_Site_df <- cbind(rownames(SPEI_Site_df), SPEI_Site_df)
rownames(SPEI_Site_df) <- NULL
names(SPEI_Site_df)[1]<- "Year"
SPEI_Site_df<-data.frame(t(SPEI_Site_df))
names(SPEI_Site_df)[1:113] <- paste(1901:2013)
SPEI_Site_df <- SPEI_Site_df[-c(1), ]
SPEI_Site_df$Site_ID<- Site_list$Site_ID
SPEI_Site_df<- gather(SPEI_Site_df, variable, value, -Site_ID)
colnames(SPEI_Site_df)<- c("Site_ID", "year", "SPEI_RS") 
SPEI_Site_df$year <- as.numeric(levels(SPEI_Site_df$year))
SPEI_Site_df$SPEI_RS <- as.numeric(SPEI_Site_df$SPEI_RS)

# 4.4 Include SPEI data into flux dataframe
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
  MAP_Site[[i]]<- MAP_Site[[i]] %>% gather(date, value, -`year 72`)
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
SPI_Site<- spi(SPI_Site_df$MAP, 6, na.rm=T)
SPI_Site_df$SPI_CRU<- SPI_Site$fitted
}

# Average SPI per year
SPI_Site_df<- ddply(SPI_Site_df, .(ID, year),
                      summarize,
                      SPI_CRU= mean(SPI_CRU, na.rm=T))

# Add Site ID to dataset
dir <- file.path(path, 'GPP_Extraction')
list <- list.files(dir, pattern=glob2rx('*.csv'), full.names=TRUE)
Site_ID<- read.csv(list, header=T,  check.names=FALSE, sep=",")
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
  MAT_Site[[i]]<- MAT_Site[[i]] %>% gather(date, value, -`year 72`)
  MAT_Site[[i]]<- MAT_Site[[i]] %>% separate(date, c("year", "month"))
  colnames(MAT_Site[[i]])<- c("ID", "year", "month", "MAT")
  MAT_Site[[i]]$ID <- as.character(MAT_Site[[i]]$ID)
  MAT_Site[[i]]$ID <- sapply(MAT_Site[[i]]$ID, function(id) as.numeric(strsplit(id,"[.]")[[1]][2]))
  MAT_Site[[i]] <- MAT_Site[[i]][with(MAT_Site[[i]], order(ID)),]
}

# Compute average temperature per site since 1900
Mean_MAT_Site<- ddply(do.call("rbind", MAT_Site), .(ID),
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
dfAll_Sites<- merge(dfAll_Sites, MAT_CRU_Site, by=c("Site_ID", "year"), all.x=TRUE)

# Compute temperature anomalies
dfAll_Sites$MAT_An<- dfAll_Sites$MAT_CRU - dfAll_Sites$MAT_Mean

# 7. Compute CEC

# 7.1 Import Site coordinate
NEP<- readRDS("Output/NEP.rds")
Site_list<-ddply(NEP, .(Site_ID),
                 summarise,
                 y= mean(Lat, na.rm=T),
                 x=mean(Long, na.rm=T))
Site.xy<- Site_list[c("x", "y")]

#7.2 SoilGrid 1km raster
dir <- file.path(path, 'Soil_Data/1km')
file <- list.files(dir, pattern=glob2rx('*.tif'), full.names=TRUE)
Soil_1km_Raster<- stack(file)

# Extract CEC value per site

# Extract data and create dataframe
Soil_1km_Data_Site <- data.frame(coordinates(Site.xy),
                   Site_list$Site_ID, 
                   extract(Soil_1km_Raster, Site.xy))

# Compute CEC total
wt<-c(0.05, 0.10, 0.15, 0.3, 0.4, 1)
Soil_1km_Data_Site$CEC_Total_1km<-apply(Soil_1km_Data_Site[,-c(1:3)],1,weighted.mean,w=wt)
Soil_1km_Data_Site<- Soil_1km_Data_Site[c("Site_list.Site_ID", "CEC_Total_1km")]
colnames(Soil_1km_Data_Site)<- c("Site_ID", "CEC_Total_1km")

# Add CEC per site in flux dataframe
dfAll_Sites<- merge(dfAll_Sites, Soil_1km_Data_Site, by=c("Site_ID"), all.x=TRUE)

# 8. Create dataframe for post-processing

# 8.1 Subset original data set into fluxes
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]
NEP<- Flux_High[Flux_High$Type_Flux %in% c("NEP"),]
GPP<- Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco<- Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
Ratio_NEP_GPP<- Flux_High[Flux_High$Type_Flux %in% c("NEP_GPP"),]
Ratio_GPP_Reco<- Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]

#8.2 Compute GPPclimax based on the lieth model
params <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$Tair))))
GPP$GPPp<-with(as.list(params), GPP1000*((1-exp(-k*GPP$Annual_Preci))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))

# 8.3 Compute ratio NEP-GPPclimax ratio
GPP$NEP_GPPmax<- NEP$values/GPP$GPPmax
Ratio_NEP_GPPmax<- GPP
Ratio_NEP_GPPmax<-Ratio_NEP_GPPmax[, !(colnames(Ratio_NEP_GPPmax) %in% c("values", "Type_Flux", "GPPmat", "GPPmax", "GPPp"))]
Ratio_NEP_GPPmax<-gather(Ratio_NEP_GPPmax, variable, values, -Annual_Preci, -year, 
                         -Ecosystem, -Climate, -Disturbance,
                         -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
                         -Tair, -Rg, -Study, -Lat, -Long, -SPI_CRU, -SPEI_RS, -MAT_An, -CEC_Total_1km)
Ratio_NEP_GPPmax<- Ratio_NEP_GPPmax[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "values", "Annual_Preci", 
                                      "Tair","Rg", "Stand_Age", "SPI_CRU", "SPEI_RS", "MAT_An", "CEC_Total_1km", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")]
colnames(Ratio_NEP_GPPmax)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", 
                              "Tair", "Rg", "Stand_Age", "SPI_CRU", "SPEI_RS", "MAT_An", "CEC_Total_1km", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")

# 8.4 Create average site dataframe
NEP_Mean_Site<-ddply(NEP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Annual_Preci=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T), 
                     SPI_CRU=mean(SPI_CRU, na.rm=T),
                     SPEI_RS=mean(SPEI_RS, na.rm=T),
                     MAT_An=mean(MAT_An, na.rm=T),
                     CEC_Total_1km= mean(CEC_Total_1km, na.rm=T))

GPP_Mean_Site<-ddply(GPP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Annual_Preci=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T), 
                     SPI_CRU=mean(SPI_CRU, na.rm=T),
                     SPEI_RS=mean(SPEI_RS, na.rm=T),
                     MAT_An=mean(MAT_An, na.rm=T),
                     CEC_Total_1km= mean(CEC_Total_1km, na.rm=T))

Reco_Mean_Site<-ddply(Reco, .(Site_ID, Type_Flux),
                      summarise,
                      Stand_Age= mean(Stand_Age, na.rm=T),
                      values=mean(values, na.rm=T),
                      Annual_Preci=mean(Annual_Preci, na.rm=T),
                      Tair=mean(Tair, na.rm=T), 
                      SPI_CRU=mean(SPI_CRU, na.rm=T),
                      SPEI_RS=mean(SPEI_RS, na.rm=T),
                      MAT_An=mean(MAT_An, na.rm=T),
                      CEC_Total_1km= mean(CEC_Total_1km, na.rm=T))

Ratio_NEP_GPP_Mean_Site<-ddply(Ratio_NEP_GPP, .(Site_ID, Type_Flux),
                               summarise,
                               Stand_Age= mean(Stand_Age, na.rm=T),
                               values=mean(values, na.rm=T),
                               Annual_Preci=mean(Annual_Preci, na.rm=T),
                               Tair=mean(Tair, na.rm=T), 
                               SPI_CRU=mean(SPI_CRU, na.rm=T),
                               SPEI_RS=mean(SPEI_RS, na.rm=T),
                               MAT_An=mean(MAT_An, na.rm=T),
                               CEC_Total_1km= mean(CEC_Total_1km, na.rm=T))

Ratio_GPP_Reco_Mean_Site<-ddply(Ratio_GPP_Reco, .(Site_ID, Type_Flux),
                                summarise,
                                Stand_Age= mean(Stand_Age, na.rm=T),
                                values=mean(values, na.rm=T),
                                Annual_Preci=mean(Annual_Preci, na.rm=T),
                                Tair=mean(Tair, na.rm=T), 
                                SPI_CRU=mean(SPI_CRU, na.rm=T),
                                SPEI_RS=mean(SPEI_RS, na.rm=T),
                                MAT_An=mean(MAT_An, na.rm=T),
                                CEC_Total_1km= mean(CEC_Total_1km, na.rm=T))

Ratio_NEP_GPPmax_Mean_Site<-ddply(Ratio_NEP_GPPmax, .(Site_ID, Type_Flux),
                                  summarise,
                                  Stand_Age= median(Stand_Age, na.rm=T),
                                  values=mean(values, na.rm=T),
                                  Annual_Preci=mean(Annual_Preci, na.rm=T),
                                  Tair=mean(Tair, na.rm=T), 
                                  SPI_CRU=mean(SPI_CRU, na.rm=T),
                                  SPEI_RS=mean(SPEI_RS, na.rm=T),
                                  MAT_An=mean(MAT_An, na.rm=T),
                                  CEC_Total_1km= mean(CEC_Total_1km, na.rm=T))

# 8.5 Save all dataframe in a rds files
saveRDS(dfAll_Sites, file="Output/df_Annual_Flux.rds")
saveRDS(NEP, file="Output/NEP.rds")
saveRDS(NEP_Mean_Site, file="Output/NEP_Mean_Site.rds")
saveRDS(GPP, file="Output/GPP.rds")
saveRDS(GPP_Mean_Site, file="Output/GPP_Mean_Site.rds")
saveRDS(Reco, file="Output/Reco.rds")
saveRDS(Reco_Mean_Site, file="Output/Reco_Mean_Site.rds")
saveRDS(Ratio_GPP_Reco, file="Output/Ratio_GPP_Reco.rds")
saveRDS(Ratio_GPP_Reco_Mean_Site, file="Output/Ratio_GPP_Reco_Mean_Site.rds")
saveRDS(Ratio_NEP_GPP, file="Output/Ratio_NEP_GPP.rds")
saveRDS(Ratio_NEP_GPP_Mean_Site, file="Output/Ratio_NEP_GPP_Mean_Site.rds")
saveRDS(Ratio_NEP_GPPmax, file="Output/Ratio_NEP_GPPmax.rds")
saveRDS(Ratio_NEP_GPPmax_Mean_Site, file="Output/Ratio_NEP_GPPmax_Mean_Site.rds")