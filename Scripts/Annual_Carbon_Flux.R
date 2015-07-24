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

#1. Create a df for all fluxnet sites dataframe

#Open all files
dir <- file.path(path, 'Fluxnet_Data')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = fLoadFluxNCIntoDataframe(VarList.V.s=c('NEE_f','GPP_f','Reco','NEE_fqcOK'),
                                          FileName.s=list[i], NcPackage.s = 'RNetCDF')
}

# Remove data before disturbance from the list of dataframe
Site_Date<-read.csv("Input/Potential_Sites.csv", header = TRUE)
Site_Date$Date_Disturbance<- as.Date(Site_Date$Date_Disturbance)
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]$DateTime=as.Date(Fluxnet_Site[[i]]$DateTime)
}

for (i in seq_along(Fluxnet_Site)){
Fluxnet_Site[[i]]<-with(Fluxnet_Site[[i]], Fluxnet_Site[[i]][(Fluxnet_Site[[i]]$DateTime > Site_Date$Date_Disturbance[i]),])
}

#3. Analyse annual carbon flux

# Compute sum/mean annual +/- sd
Sum_Sd_Flux<-lapply(Fluxnet_Site, function (x) ddply(x, "year",
                    summarize,
                    sum_NEE= sum(NEE_f, na.rm=T),
                    sum_GPP= sum(GPP_f, na.rm=T),
                    sum_TER= sum(Reco,na.rm=T),
                    mean_NEE= mean(NEE_f, na.rm=T),
                    mean_GPP= mean(GPP_f, na.rm=T),
                    mean_TER= mean(Reco,na.rm=T),
                    sd_NEE= sd(NEE_f, na.rm=T),
                    sd_GPP= sd(GPP_f, na.rm=T),
                    sd_Reco= sd(Reco, na.rm=T),
                    mean_Uncert=mean(NEE_fqcOK, na.rm=T)))

#Add type of ecosystem, type of climate, type of disturbance and year after disturbance for each sites
for (i in seq_along(Sum_Sd_Flux)){
  Sum_Sd_Flux[[i]]$Ecosystem<- Site_Date$Ecosystem[i]
  Sum_Sd_Flux[[i]]$Climate<- Site_Date$Climate[i]
  Sum_Sd_Flux[[i]]$Disturbance<- Site_Date$Type_Disturbance[i]
  Sum_Sd_Flux[[i]]$Year_Disturbance<- Sum_Sd_Flux[[i]]$year- year(Site_Date$Date_Disturbance[i])
  Sum_Sd_Flux[[i]]$GPP_ER<- Sum_Sd_Flux[[i]]$sum_GPP / Sum_Sd_Flux[[i]]$sum_TER
  Sum_Sd_Flux[[i]]$Site_ID<- Site_Date$ID[i]
  Sum_Sd_Flux[[i]]$Stand_Replacement<- Site_Date$Stand_Replacement[i]
  Sum_Sd_Flux[[i]]$Int_Replacement<- Site_Date$Intensity_Replacement[i]
}

#Combine the flux sites in one dataframe
dfAll_Sites<- do.call("rbind", Sum_Sd_Flux)
dfAll_Sites<-dfAll_Sites[-c(122),]# the measurements seems to be an outlier

# Remove high gap filled fraction
dfAll_Sites<- dfAll_Sites[dfAll_Sites$mean_Uncert>0.80,]

# Restructure dataframe
dfAll_Sites<-gather(dfAll_Sites, Type_Flux, values, -year, -Ecosystem, -mean_Uncert, -Climate, -Disturbance, -Year_Disturbance, -Site_ID, -Stand_Replacement, -Int_Replacement)

#Reoder column
dfAll_Sites<- dfAll_Sites[c("Site_ID", "year", "Stand_Replacement","Int_Replacement","Type_Flux", "values", "mean_Uncert", "Year_Disturbance", "Disturbance", "Climate", "Ecosystem")]

# Reclassify climate classification
dfAll_Sites$Climate<-ifelse((dfAll_Sites$Climate =="Af" | dfAll_Sites$Climate =="Am"), "Tropical",
                    ifelse((dfAll_Sites$Climate =="Cfa" | dfAll_Sites$Climate =="Cfb" | dfAll_Sites$Climate =="Cfc" | dfAll_Sites$Climate =="Csa" | dfAll_Sites$Climate =="Csb"), "Temperate",
                           ifelse((dfAll_Sites$Climate =="Dfb" | dfAll_Sites$Climate =="Dfc" | dfAll_Sites$Climate =="Dsc"), "Continental", NA)))

#.4.Function fit choice for ecosystem response
GPP<-dfAll_Sites[dfAll_Sites$Type_Flux %in% c("sum_GPP"),]
GPP_High_Dist<-GPP[GPP$Stand_Replacement ==1,]
GPP_High_Dist<-GPP_High_Dist[GPP_High_Dist$Int_Replacement %in% "High",]
Harvest_GPP<-GPP_High_Dist[GPP_High_Dist$Disturbance %in% c("Harvest"),]
Fire_GPP<-GPP_High_Dist[GPP_High_Dist$Disturbance %in% c("Wildfire"),]

# 4.1 Gamma function

# Conpute parameters of the function
plotPoints(values ~ Year_Disturbance, data=Harvest_GPP)
f1= fitModel(values~A*(Year_Disturbance^B)*(exp(k*Year_Disturbance)), data=Fire_GPP, start = list(A=1000, B=0.170, k= -0.00295))
coef(f1)
plotFun(f1(Year_Disturbance)~Year_Disturbance, Year_Disturbance.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f1(Harvest_GPP$Year_Disturbance)~Harvest_GPP$values)
summary(fit)

#Compute function
Gamma_Harvest <- function(x){495.98872446 *(x^0.40451111)*(exp(-0.01201367*x))}
Gamma_Fire<- function(x){152.237501652*(x^0.482614951)*(exp(-0.003643999*x))}

# 4.2 Ricker function

# Conpute parameters of the function
plotPoints(values ~ Year_Disturbance, data=Harvest_GPP)
f2= fitModel(values~A*(Year_Disturbance*(exp((k*Year_Disturbance)/2))), data=Harvest_GPP, start = list(A=500, k= -0.3))
coef(f2)
plotFun(f2(Year_Disturbance)~Year_Disturbance, Year_Disturbance.lim=range(0,120), add=T)

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

# 4.5. Spline function
smooth.spline2 <- function(formula, data, ...) { 
  mat <- model.frame(formula, data) 
  smooth.spline(mat[, 2], mat[, 1]) 
} 
predictdf.smooth.spline <- function(model, xseq, se, level) { 
  pred <- predict(model, xseq) 
  data.frame(x = xseq, y = pred$y) 
}

# 4.6 Amiro function for ratio GPP and Reco
GPP_Reco<-dfAll_Sites[dfAll_Sites$Type_Flux %in% c("GPP_ER"),]
GPP_Reco<-na.omit(GPP_Reco)
Harvest_GPP_Reco<-GPP_Reco[GPP_Reco$Disturbance %in% c("Harvest"),]
Fire_GPP_Reco<-GPP_Reco[GPP_Reco$Disturbance %in% c("Fire"),]

# Conpute parameters of the function
plotPoints(values ~ Year_Disturbance, data=Harvest_GPP_Reco)
f6= fitModel(values~A*(1-exp(k*Year_Disturbance)), data=Harvest_GPP_Reco, start = list(A=1.23, k= -0.224))
coef(f6)
plotFun(f6(Year_Disturbance)~Year_Disturbance, Year_Disturbance.lim=range(0,120), add=T)

#Compute r2
fit<- lm(f4(Fire_GPP_Reco$Year_Disturbance)~Fire_GPP_Reco$values)
summary(fit)

#Compute function
Amiro_Harvest <- function(x){1.2173626 *(1-exp(-0.5255131*x))}
Amiro_Fire<- function(x){1.0734594*(1-exp(-0.2010539 *x))}

# 4.5 Exponential function for ratio GPP and Reco

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
Plot_Flux<- dfAll_Sites[dfAll_Sites$Type_Flux %in% c("sum_GPP", "sum_TER"),]
Plot_Flux_Harvest<- Plot_Flux[Plot_Flux$Disturbance %in% c("Harvest"),]
Plot_Flux_Storm<- Plot_Flux[Plot_Flux$Disturbance %in% c("Storm"),]
Plot_Flux_Fire<- Plot_Flux[Plot_Flux$Disturbance %in% c("Wildfire"),]
Plot_Flux_ENF<- Plot_Flux_Harvest[Plot_Flux_Harvest$Ecosystem %in% c("ENF"),]
Plot_Flux_DBF<- Plot_Flux_Harvest[Plot_Flux_Harvest$Ecosystem %in% c("DBF"),]

gg1<-ggplot(Plot_Flux_Fire, aes(Year_Disturbance, values, colour=mean_Uncert)) +
  geom_point(size=3) +
  facet_grid(Type_Flux~Disturbance, scales = "free")+
  # geom_smooth(method="smooth.spline2", se=F, color="red")+
  stat_function(fun=Ricker_Fire, color="red")+
  # geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("Year since Disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(name="Fraction of gap-filled data", low="#FF0000", high = "#00FF33")

# 5.2. Ratio GPP and Reco

#Subset data
gg2<-ggplot(Harvest_GPP_Reco,aes(Year_Disturbance, values, colour=mean_Uncert)) +
  geom_point (size=3.5)+
  facet_wrap(~Disturbance, ncol=1, scales="free")+
  stat_function(fun=Amiro_Harvest, color="red") +
  expand_limits(x = 0, y = 0)+
  geom_hline(yintercept=1, linetype=2, colour="grey", size=0.7)+
  # geom_errorbar(limits_TER, width=0.07, linetype=6)+
  xlab("Year since disturbance") + ylab("GPP/Reco")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(name="Fraction of gap-filled data", low="#FF0000", high = "#00FF33")

#Export plot
ggsave(gg7, filename = 'Latex/Figures/All_Sites/Dist_GPP_TER.eps', width = 14, height = 8)

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

