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

#Remove outlier sites (US-Blo, CZ-Bk1, US-DK3, IT_Ro1, IT-Ro2, US_Bns, US-NC2)
Sum_Sd_Flux<- Sum_Sd_Flux[-c(28,36,37,49,50,51,53,63)]

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

#Subset dataframe for fitting process 
Flux_High<-dfAll_Sites[dfAll_Sites$Int_Replacement %in% c("High"),]
GPP_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco_High<-Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
GPP_High_Harvest<-GPP_High[GPP_High$Disturbance %in% c("Harvest"),]
Reco_High_Harvest<-Reco_High[Reco_High$Disturbance %in% c("Harvest"),]
GPP_High_Fire<-GPP_High[GPP_High$Disturbance %in% c("Wildfire"),]
Reco_High_Fire<-Reco_High[Reco_High$Disturbance %in% c("Wildfire"),]
Ratio_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]
Ratio_High_Harvest<-Ratio_High[Ratio_High$Disturbance %in% c("Harvest"),]
Ratio_High_Fire<-Ratio_High[Ratio_High$Disturbance %in% c("Wildfire"),]

# 4.1 Gamma function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=GPP_High_Fire)
f1= fitModel(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data=GPP_High_Fire, start = list(A=1000, B=0.170, k= -0.00295))
coef(f1)
plotFun(f1(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,200), add=T)

#Compute r2
fit<- lm(f1(GPP_High_Fire$Stand_Age)~GPP_High_Fire$values)
summary(fit)

# 4.2 Ricker function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=GPP_High_Fire)
f2= fitModel(values~A*(Stand_Age*(exp((k*Stand_Age)/2))), data=GPP_High_Fire, start = list(A=500, k= -0.3))
coef(f2)
plotFun(f2(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f2(GPP_High_Fire$Stand_Age)~GPP_High_Fire$values)
summary(fit)

# 4.3 Second order polymonial function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=GPP_High_Fire)
f3= fitModel(values~A*Stand_Age^2+B*Stand_Age+C, data=GPP_High_Fire)
coef(f3)
plotFun(f3(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f3(GPP_High_Fire$Stand_Age)~GPP_High_Fire$values)
summary(fit)

# 4.4 Third order polymonial function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=GPP_High_Fire)
f4= fitModel(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data=GPP_High_Fire)
coef(f4)
plotFun(f4(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f4(GPP_High_Fire$Stand_Age)~GPP_High_Fire$Stand_Age)
summary(fit)

# 4.5 Amiro function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=GPP_High_Fire)
f5= fitModel(values~A*(1-exp(k*Stand_Age)), data=GPP_High_Fire, start = list(A=1500, k= -0.224))
coef(f5)
plotFun(f5(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f5(GPP_High_Fire$Stand_Age)~GPP_High_Fire$values)
summary(fit)

# 5. Plot ecosystem response with fit function

# 5.1 Annual carbon flux

#Compute fit function with the best fit
Harvest_GPP_Fun<-function(x){-0.4392582*x^2 + 48.6588570*x + 260.8064506}
Harvest_Reco_Fun<-function(x){-0.2554387*x^2 +26.4671995*x + 604.9361162}
Fire_GPP_Fun<- function(x){13.83017087*(x^1.40517046)*(exp(-0.02274242*x))}
Fire_Reco_Fun<- function(x){39.82891365*x*(exp(-0.03179599*x/2))}

#Plot data
#GPP Harvest
gg1<-ggplot(GPP_High_Harvest, aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 6)) +
  facet_grid(Type_Flux~Disturbance, scales="free_x")+
  stat_function(fun=Harvest_GPP_Fun, color="red")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.80,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")+
  ylim(0,3000)

print(gg1)

#Reco Harvest
gg2<-ggplot(Reco_High_Harvest, aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 6)) +
  facet_grid(Type_Flux~Disturbance, scales="free_x")+
  stat_function(fun=Harvest_Reco_Fun, color="red")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.80,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")+
  ylim(0,3000)

print(gg2)

#GPP Wildfire
gg3<-ggplot(GPP_High_Fire, aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 8)) +
  facet_grid(Type_Flux~Disturbance, scales="free_x")+
  stat_function(fun=Fire_GPP_Fun, color="red")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.80,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")+
  ylim(0,3000)

print(gg3)

#Reco Wildfire
gg4<-ggplot(Reco_High_Fire, aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 8)) +
  facet_grid(Type_Flux~Disturbance, scales="free_x")+
  stat_function(fun=Fire_Reco_Fun, color="red")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.80,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")+
  ylim(0,3000)

print(gg4)

# 5.2. Ratio GPP and Reco

#Compute fit function with the best fit
Harvest_Ratio_Fun<-function(x){0.396598614*(x^0.387995395)*(exp(-0.007407526*x))}
Fire_Ratio_Fun<-function(x){0.194575390*(x^0.537502454)*(exp(-0.006993452*x))}

#Plot data
#Harvest
gg5<-ggplot(Ratio_High_Harvest,aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 6)) +
  facet_grid(Type_Flux~Disturbance, scales="free_x")+
  stat_function(fun=Harvest_Ratio_Fun, color="red") +
  expand_limits(x = 0, y = 0)+
  geom_hline(yintercept=1, linetype=2, colour="grey", size=0.7)+
  xlab("Year since disturbance") + ylab("GPP/Reco")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.80,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")+
  ylim(0,2.5)

print(gg5)

#Wildfire
gg6<-ggplot(Ratio_High_Fire,aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 8)) +
  facet_grid(Type_Flux~Disturbance, scales="free_x")+
  stat_function(fun=Fire_Ratio_Fun, color="red") +
  expand_limits(x = 0, y = 0)+
  geom_hline(yintercept=1, linetype=2, colour="grey", size=0.7)+
  xlab("Year since disturbance") + ylab("GPP/Reco")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.80,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")+
  ylim(0,2.5)

print(gg6)

#5.3.Plot carbon flux and GPP/Reco ratio together
gA <- ggplotGrob(gg1)
gB <- ggplotGrob(gg3)
gC <- ggplotGrob(gg2)
gD <- ggplotGrob(gg4)
gE <- ggplotGrob(gg5)
gF <- ggplotGrob(gg6)

maxWidth = grid::unit.pmax(gA$widths[1:2], gB$widths[1:2], gC$widths[1:2], gD$widths[1:2], gE$widths[1:2], gF$widths[1:2])
gA$widths[1:2] <- as.list(maxWidth)
gB$widths[1:2] <- as.list(maxWidth)
gC$widths[1:2] <- as.list(maxWidth)
gD$widths[1:2] <- as.list(maxWidth)
gE$widths[1:2] <- as.list(maxWidth)
gF$widths[1:2] <- as.list(maxWidth)

gg7 <- arrangeGrob(
  gA,gB,gC,gD,gE,gF, nrow = 3, heights = c(0.5, 0.5))

print(gg7)

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

