## Script to analyse annual carbon flux
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library(RNetCDF)
library (ggplot2)
library(scales)
library(lubridate)
library(gridExtra)
library(REddyProc)
library (dplyr)
library (plyr)
library(gtools)
library (RColorBrewer)
library(gamm4)
library(tidyr)
library(manipulate)
library(MASS)
library(fitdistrplus)
library(mosaic)
library (maps)
library(ggsubplot)
library (ggmap)
library(OpenStreetMap)
library(devtools)
library (ggvis)
library(gridExtra)

#1. Create a df for all fluxnet sites dataframe

#Open all files
dir <- file.path(path, 'Fluxnet_Data')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = fLoadFluxNCIntoDataframe(VarList.V.s=c('NEE_f','GPP_f','Reco','NEE_fqcOK', 'precip_f_Night_m1', 'precip_f_Day_m1'),
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
                    mean_NEE= mean(NEE_f, na.rm=T),
                    mean_GPP= mean(GPP_f, na.rm=T),
                    mean_TER= mean(Reco,na.rm=T),
                    sd_NEE= sd(NEE_f, na.rm=T),
                    sd_GPP= sd(GPP_f, na.rm=T),
                    sd_Reco= sd(Reco, na.rm=T),
                    mean_Uncert=mean(NEE_fqcOK, na.rm=T), 
                    Annual_Pre_Day= sum(precip_f_Day_m1, na.rm=T),
                    Annual_Pre_Night= sum(precip_f_Night_m1, na.rm=T)))

#Add type of ecosystem, type of climate, type of disturbance and year after disturbance for each sites
for (i in seq_along(Sum_Sd_Flux)){
  Sum_Sd_Flux[[i]]$Ecosystem<- Site_Date$Ecosystem[i]
  Sum_Sd_Flux[[i]]$Climate<- Site_Date$Climate[i]
  Sum_Sd_Flux[[i]]$Disturbance<- Site_Date$Type_Disturbance[i]
  Sum_Sd_Flux[[i]]$Year_Disturbance<- Sum_Sd_Flux[[i]]$year- year(as.Date(Site_Date$Date_Disturbance[i]))
  Sum_Sd_Flux[[i]]$Stand_Age<- Sum_Sd_Flux[[i]]$year- year(as.Date(Site_Date$Plantation_Date[i]))
  Sum_Sd_Flux[[i]]$GPP_ER<- Sum_Sd_Flux[[i]]$GPP / Sum_Sd_Flux[[i]]$Respiration
  Sum_Sd_Flux[[i]]$Site_ID<- Site_Date$ID[i]
  Sum_Sd_Flux[[i]]$Stand_Replacement<- Site_Date$Stand_Replacement[i]
  Sum_Sd_Flux[[i]]$Int_Replacement<- Site_Date$Intensity_Replacement[i]
  Sum_Sd_Flux[[i]]$Age_Min<- Site_Date$Age_Forest_Min[i]
  Sum_Sd_Flux[[i]]$Age_Max<- Site_Date$Age_Forest_Max[i]
  Sum_Sd_Flux[[i]]$Annual_Preci<-Sum_Sd_Flux[[i]]$Annual_Pre_Day + Sum_Sd_Flux[[i]]$Annual_Pre_Night
}

#Remove outlier sites (US-Bn, US-Blo, CZ-Bk1, IT-Sro, US-Me1, US-DK3)
Sum_Sd_Flux<- Sum_Sd_Flux[-c(28,49,50,51,52,53,61)]

#Combine the flux sites in one dataframe
dfAll_Sites<- do.call("rbind", Sum_Sd_Flux)
dfAll_Sites<-dfAll_Sites[-c(127),]# the measurements seems to be an outlier

# Remove years missing in each flux site
NA_Value <- apply(dfAll_Sites, 1, function(x){any(is.na(x))})
dfAll_Sites<- dfAll_Sites[!NA_Value,]

# Remove high gap filled fraction
dfAll_Sites<- dfAll_Sites[dfAll_Sites$mean_Uncert>0.85,]

# Restructure dataframe
dfAll_Sites<-gather(dfAll_Sites, Type_Flux, values, -Annual_Preci, -year, -Ecosystem, -mean_Uncert, -Climate, -Disturbance, -Year_Disturbance, -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement)

#Reoder column
dfAll_Sites<- dfAll_Sites[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", "mean_Uncert", "Year_Disturbance", "Stand_Age", "Disturbance", "Climate", "Ecosystem")]

# Reclassify climate classification
dfAll_Sites$Climate<-ifelse((dfAll_Sites$Climate =="Af" | dfAll_Sites$Climate =="Am"), "Tropical",
                    ifelse((dfAll_Sites$Climate =="Cfa" | dfAll_Sites$Climate =="Cfb" | dfAll_Sites$Climate =="Cfc" | dfAll_Sites$Climate =="Csa" | dfAll_Sites$Climate =="Csb"), "Temperate",
                           ifelse((dfAll_Sites$Climate =="Dfb" | dfAll_Sites$Climate =="Dfc" | dfAll_Sites$Climate =="Dsc"), "Continental", NA)))

#.4.Function fit choice for ecosystem response
GPP<-dfAll_Sites[dfAll_Sites$Type_Flux %in% c("GPP"),]
GPP_High_Dist<-GPP[GPP$Int_Replacement %in% "High",]
Harvest_GPP<-GPP_High_Dist[GPP_High_Dist$Disturbance %in% c("Harvest"),]
Fire_GPP<-GPP_High_Dist[GPP_High_Dist$Disturbance %in% c("Wildfire"),]

# 4.1 Gamma function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=Plot_Flux_Fire)
f1= fitModel(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data=Plot_Flux_Fire, start = list(A=1000, B=0.170, k= -0.00295))
coef(f1)
plotFun(f1(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,200), add=T)

#Compute r2
fit<- lm(f1(Flux_Ratio$Stand_Age)~Flux_Ratio$values)
summary(fit)

#Compute function
Gamma_Harvest <- function(x){495.98872446 *(x^0.40451111)*(exp(-0.01201367*x))}
Gamma_Fire<- function(x){152.237501652*(x^0.482614951)*(exp(-0.003643999*x))}

# 4.2 Ricker function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=Flux_Ratio)
f2= fitModel(values~A*(Stand_Age*(exp((k*Stand_Age)/2))), data=Flux_Ratio, start = list(A=500, k= -0.3))
coef(f2)
plotFun(f2(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f2(Harvest_GPP$Year_Disturbance)~Harvest_GPP$values)
summary(fit)

#Compute function
Ricker_Harvest<- function(x){99.38607842*x*(exp(-0.04625329*x/2))}
Ricker_Fire<- function(x){29.91918974*x*(exp(-0.03419415*x/2))}

# 4.3 Second order polymonial function

# Conpute parameters of the function
plotPoints(values ~ Year_Disturbance, data=Fire_GPP)
f3= fitModel(values~A*Year_Disturbance^2+B*Year_Disturbance+C, data=Fire_GPP)
coef(f3)
plotFun(f3(Year_Disturbance)~Year_Disturbance, Year_Disturbance.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f3(Fire_GPP$Year_Disturbance)~Fire_GPP$values)
summary(fit)

#Compute function
Poly2_Harvest<- function(x){-0.5958727*x^2 + 56.0966845*x + 783.6275529}
Poly2_Fire<- function(x){-0.08703708*x^2 + 8.26272839*x + 409.20396672}

# 4.4 Third order polymonial function

# Conpute parameters of the function
plotPoints(values ~ Year_Disturbance, data=Fire_GPP)
f4= fitModel(values~A*Year_Disturbance^3+B*Year_Disturbance^2+C*Year_Disturbance+D, data=Fire_GPP)
coef(f4)
plotFun(f4(Year_Disturbance)~Year_Disturbance, Year_Disturbance.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f4(Fire_GPP$Year_Disturbance)~Fire_GPP$values)
summary(fit)

#Compute function
Poly3_Harvest<- function(x){0.009480182*x^3 -2.015796171*x^2 + 108.741923329*x + 523.052930475}
Poly3_Fire<- function(x){0.004709568*x^3 -0.657932473*x^2 + 24.532387387*x + 321.946134971}

# 4.5 Amiro function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=Flux_Ratio)
f6= fitModel(values~A*(1-exp(k*Stand_Age)), data=Flux_Ratio, start = list(A=1.23, k= -0.224))
coef(f6)
plotFun(f6(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f4(Fire_GPP_Reco$Year_Disturbance)~Fire_GPP_Reco$values)
summary(fit)

#Compute function
Amiro_Harvest <- function(x){1.2173626 *(1-exp(-0.5255131*x))}
Amiro_Fire<- function(x){1.0734594*(1-exp(-0.2010539 *x))}

# 4.6 Exponential function

# Conpute parameters of the function
f5= fitModel(values~A+B*exp(k*Year_Disturbance), data=Fire_GPP_Reco, start = list(A=5, B=-5, k=-0.005))
coef(f5)
plotFun(f5(Year_Disturbance)~Year_Disturbance, Year_Disturbance.lim=range(0,120), add=T)

#Compute function
Expo_Harvest <- function(x){1.23+exp(-0.224*x)}
Expo_Fire<- function(x){20.730+9.24*exp(1*x)}
plot(Expo_Fire, xlim=c(0,100))

# 5. Plot ecosystem response with fit function

# 5.1 Partition flux data per variable

#Subset data
Plot_Flux<- dfAll_Sites[dfAll_Sites$Type_Flux %in% c("GPP", "Respiration"),]
Plot_Flux_Harvest<- Plot_Flux[Plot_Flux$Disturbance %in% c("Harvest"),]
Plot_Flux_Fire<- Plot_Flux[Plot_Flux$Disturbance %in% c("Wildfire"),]

#Plot data
#Harvest
gg1<-ggplot(Plot_Flux_Harvest, aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 6)) +
  facet_wrap(~Type_Flux, ncol=1)+
  facet_grid(Type_Flux~Disturbance)+
  stat_function(fun=Gamma_Flux_Harvest, color="red")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.85,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")

#Wildfire
gg2<-ggplot(Plot_Flux_Fire, aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3,8 )) +
  facet_wrap(~Type_Flux, ncol=1)+
  facet_grid(Type_Flux~Disturbance)+
  stat_function(fun=Gamma_Flux_Fire, color="red")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.85,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")

# 5.2. Ratio GPP and Reco

#Subset data
Flux_Ratio<- dfAll_Sites[dfAll_Sites$Type_Flux %in% c("GPP_ER"),]
Flux_Ratio<- Flux_Ratio[Flux_Ratio$Int_Replacement %in% c("High"),]
Flux_Ratio_Harvest<- Flux_Ratio[Flux_Ratio$Disturbance %in% c("Harvest"),]
Flux_Ratio_Wildfire<- Flux_Ratio[Flux_Ratio$Disturbance %in% c("Wildfire"),]

#Plot data
#Harvest
gg3<-ggplot(Flux_Ratio_Harvest,aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 6)) +
  facet_wrap(~Disturbance, ncol=1)+
  stat_function(fun=Gamma_ratio_Harvest, color="red") +
  expand_limits(x = 0, y = 0)+
  geom_hline(yintercept=1, linetype=2, colour="grey", size=0.7)+
  xlab("Year since disturbance") + ylab("GPP/Reco")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.85,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")

#Wildfire
gg4<-ggplot(Flux_Ratio_Wildfire,aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 8)) +
  facet_wrap(~Disturbance, ncol=1)+
  stat_function(fun=Gamma_ratio_Fire, color="red") +
  expand_limits(x = 0, y = 0)+
  geom_hline(yintercept=1, linetype=2, colour="grey", size=0.7)+
  xlab("Year since disturbance") + ylab("GPP/Reco")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.85,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")

#Plot flux and ratio together
gA <- ggplotGrob(gg1)
gB <- ggplotGrob(gg2)
gC <- ggplotGrob(gg3)
gD <- ggplotGrob(gg4)

maxWidth = grid::unit.pmax(gA$widths[1:2], gB$widths[1:2], gC$widths[1:2], gD$widths[1:2])
gA$widths[1:2] <- as.list(maxWidth)
gB$widths[1:2] <- as.list(maxWidth)
gC$widths[1:2] <- as.list(maxWidth)
gD$widths[1:2] <- as.list(maxWidth)

gg5 <- arrangeGrob(
  gA, gC, gB, gD, nrow = 2, heights = c(0.5, 0.5))

print(gg5)

#6 Explain variabilty of the fluxes

#6.1. Anova analysis

#Subset dataset
dfAll_Sites$Year_Disturbance<-as.factor(dfAll_Sites$Year_Disturbance)
dfAll_Sites$Climate<-as.factor(dfAll_Sites$Climate)
dfAll_Sites<-na.omit(dfAll_Sites)
An_GPP<- dfAll_Sites[dfAll_Sites$Type_Flux %in% c("sum_GPP"),]

#Compute Anova for GPP
an<-aov(log(values)~Disturbance, An_GPP)
summary(an)

#Compute Tukey test
TukeyHSD(x = an)
par(mfrow=c(2,2))
plot(an)

#Compute linear model
fit = lm(values ~ Disturbance + Climate + Year_Disturbance, data=An_GPP)

m3 = fit
m2 = update(m3, ~ . - Ecosytem)
m1 = update(m2, ~ . - Climate)
m0=  update(m1, ~ . - Year_Disturbance)

af<-anova(m0,m1,m2,m3)
summary(af)
afss <- af$"Sum of Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

#Compute performance of the model
af<-step(fit, direction = c("forward"), steps = 2000)
summary(af)


#7. Spatial map with the flux sites

# Import site information
Site_Date<-read.csv("Input/Potential_Sites.csv", header = TRUE)

# read shapefile
wmap <- readOGR(dsn="Input", layer="ne_110m_land")

# convert to dataframe
wmap_df <- fortify(wmap)

#Transform dataframe
wmap_robin <- spTransform(wmap, CRS("+proj=igh +ellps=sphere +towgs84=0,0,0 +lon_0=100w +x_0=-11119487.43"))
wmap_df_robin <- as(wmap_robin, "data.frame")

# create a blank ggplot theme
theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank(),
                         panel.background = element_blank(),
                         plot.background = element_rect(fill=NULL),
                         panel.border = element_blank(),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         plot.title = element_text(size=22)))

# add graticule and bounding box (longlat)
grat <- readOGR("Input", layer="ne_110m_graticules_15") 
grat_df <- fortify(grat)

bbox <- readOGR("Input", layer="ne_110m_wgs84_bounding_box") 
bbox_df<- fortify(bbox)

# graticule and BBox (Robin)
grat_robin <- spTransform(grat, CRS("+proj=igh +ellps=sphere +towgs84=0,0,0 +lon_0=100w +x_0=-11119487.43"))  # reproject graticule
grat_df_robin <- fortify(grat_robin)
bbox_robin <- spTransform(bbox, CRS("+proj=igh +ellps=sphere +towgs84=0,0,0 +lon_0=100w +x_0=-11119487.43"))  # reproject bounding box
bbox_robin_df <- fortify(bbox_robin)

#Fluxnet
Fluxnet_Site <- SpatialPointsDataFrame(coords=Site_Date[,c("y","x")],data=data.frame(Site_Date),proj4string=CRS(proj="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
Fluxnet_Site <- spTransform(Fluxnet_Site, CRS("+proj=igh +ellps=sphere +towgs84=0,0,0 +lon_0=100w +x_0=-11119487.43"))
Fluxnet_Site <- as(Fluxnet_Site, "data.frame")

#Plot
ggplot(bbox_robin_df, aes(long,lat, group=group)) + 
  geom_polygon(fill="white") +
  geom_polygon(data=wmap_df_robin, aes(long,lat, group=group, fill=hole)) + 
  geom_path(data=grat_df_robin, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  geom_point(data=Fluxnet_Site, aes(x=y.1,y=x.1, group= NULL),  shape=4, color="red", size=3)+
  labs(title="Fluxnet sites location") + 
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("grey", "white"), guide="none") # change colors & remove legend

