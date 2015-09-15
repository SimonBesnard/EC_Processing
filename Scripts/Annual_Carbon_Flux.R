## Script to analyse annual carbon flux
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library(RNetCDF)
library (ggplot2)
library(scales)
library(lubridate)
library(REddyProc)
library (dplyr)
library (plyr)
library(gtools)
library(tidyr)
library(manipulate)
library(MASS)
library(devtools)
library(gridExtra)
library (bootstrap)
library (foreach)
library (nlstools)
library(randomForest)
library(hexbin)
library(RColorBrewer)
library(plot3D)
library(ggRandomForests)
library(randomForestSRC)

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
                    mean_NEE= mean(NEE_f, na.rm=T),
                    mean_GPP= mean(GPP_f, na.rm=T),
                    mean_TER= mean(Reco,na.rm=T),
                    sd_NEE= sd(NEE_f, na.rm=T),
                    sd_GPP= sd(GPP_f, na.rm=T),
                    sd_Reco= sd(Reco, na.rm=T),
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
  Sum_Sd_Flux[[i]]$Year_Disturbance<- Sum_Sd_Flux[[i]]$year- year(as.Date(Site_Date$Date_Disturbance[i]))
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
dfAll_Sites<-gather(dfAll_Sites, Type_Flux, values, -Annual_Preci, -year, 
                    -Ecosystem, -mean_Uncert, -Climate, -Disturbance, 
                    -Year_Disturbance, -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
                    -Tair, -Tsoil, -SWC, -Rg)
names(dfAll_Sites)

#Reoder column
dfAll_Sites<- dfAll_Sites[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", 
                            "mean_Uncert", "Tair", "Tsoil","SWC","Rg","Year_Disturbance", "Stand_Age", "Disturbance", "Climate", "Ecosystem")]

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

# Conpute parameters of the function with bootsraping
Fun_Gamma <- nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data=GPP_High_Harvest, 
                            start = list(A=1000, B=0.170, k= -0.00295),algorithm="port")
Boo_Gamma <- nlsBoot(Fun_Gamma,niter=200)
plot(Boo_Gamma, type="boxplot")

# Store parameters combination in a new dataframe
Param_Gamma<-Boo_Gamma$coefboot

# 4.2 Ricker function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=Reco_High_Fire)
f2= fitModel(values~A*(Stand_Age*(exp((k*Stand_Age)/2))), data=Reco_High_Fire, start = list(A=500, k= -0.3))
coef(f2)
plotFun(f2(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,120), add=T)

# 4.3 Second order polymonial function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=Reco_High_Fire)
f3= fitModel(values~A*Stand_Age^2+B*Stand_Age+C, data=Reco_High_Fire)
coef(f3)
plotFun(f3(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,120), add=T)

# 4.4 Third order polymonial function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=Reco_High_Fire)
f4= fitModel(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data=Reco_High_Fire)
coef(f4)
plotFun(f4(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,170), add=T)

# 4.5 Amiro function

# Conpute parameters of the function
plotPoints(values ~ Stand_Age, data=Reco_High_Fire)
f5= fitModel(values~A*(1-exp(k*Stand_Age)), data=Reco_High_Fire, start = list(A=1500, k= -0.224))
coef(f5)
plotFun(f5(Stand_Age)~Stand_Age, Stand_Age.lim=range(0,120), add=T)

# 4.6. Compute statistical tests for the fitting functions

#R-squared
R2<- lm(f1(Reco_High_Fire$Stand_Age)~Reco_High_Fire$values)
summary(R2)

#Pearson test
Pearson<-cor.test(f1(Reco_High_Fire$Stand_Age),Reco_High_Fire$values)
print(Pearson)

#Spearman test
Spearman<-cor.test(f1(Reco_High_Fire$Stand_Age),Reco_High_Fire$values,method = "spearman")
print(Spearman)

#Kendall test
Kendall<-cor.test(f1(Reco_High_Fire$Stand_Age),Reco_High_Fire$values,method = "kendall")
print(Kendall)

#Normalized RMSE
NRMSE<-nrmse(f1(Reco_High_Fire$Stand_Age),Reco_High_Fire$values, norm="maxmin")# one can also use sd for the normalization
print(NRMSE)

#Variance Reco
VR<-var.test(f1(Reco_High_Fire$Stand_Age), Reco_High_Fire$values, Reco=1, conf.level=0.95)
print(VR)

#Modeling efficiency
MEF<-NSE(f1(Reco_High_Fire$Stand_Age), Reco_High_Fire$values)
print(MEF)

#Bootstrap (https://www.youtube.com/watch?v=cv2TSOHPp1o; https://cran.r-project.org/web/packages/broom/vignettes/bootstrapping.html; 
# https://rpubs.com/dgolicher/median_boot)

#4.7.Import df with the statistical test output and plot

#Import csv file
MEF_Out<-read.csv(file = "Output/MEF_Test.csv", header = T)

#Plot data
MEF_Out$Flux<-factor(MEF_Out$Flux,c("GPP", "Reco", "GPP/Reco"))

gg<-ggplot(MEF_Out, aes(x=Values, y=Function, label=Values)) +
  geom_point(size=3, shape=3) +
  geom_text(size=3, vjust=2)+
  facet_grid(Flux~Disturbance)+
  xlab("MEF") + ylab("Type of function")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlim(0,1)
print(gg)

# 5. Plot ecosystem response with the best fit function

# 5.1 Annual carbon flux

#Compute fit function with the best fit
Harvest_GPP_Fun<-function(x){-5.377303e-04*x^3 -3.558818e-01*x^2 + 4.543758e+01*x+2.786747e+02}
Harvest_Reco_Fun<-function(x){-0.001386029*x^3 -0.040531572*x^2 + 18.164191883*x+650.992435538}
Fire_GPP_Fun<- function(x){1.935848e-03*x^3 -5.929391e-01*x^2 + 5.042659e+01*x -1.842747e+02}
Fire_Reco_Fun<- function(x){8.987988e-04*x^3 -2.907249e-01*x^2 + 2.584671e+01*x + 2.138053e+02}

#Plot data
#GPP Harvest
gg1<-ggplot(GPP_High_Harvest, aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 6)) +
  facet_grid(Type_Flux~Disturbance, scales="free_x")+
  # geom_line(aes(y=predict(Fun_Gamma, colour="red")))+
  stat_function(fun=Harvest_GPP_Fun, color="red")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.80,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")+
  ylim(0,3000)+
  annotate("text", x = 27, y = 2800, label = "y=-5.38e-04x^3 -3.56e-01x^2 + 4.54e+01x + 2.79e+02", size=3.5)
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
  ylim(0,3000)+
  annotate("text", x = 20, y = 2800, label = "y=-0.001x^3 -0.04x^2 + 18.16x + 650.99", size=3.5)
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
  ylim(0,3000)+
  annotate("text", x = 45, y = 2800, label = "y=1.94e-03x^3 -5.93e-01x^2 + 5.04e+01x -1.84e+02", size=3.5)
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
  ylim(0,3000)+
  annotate("text", x = 45, y = 2800, label = "y=8.99e-04x^3 -2.91e-01x^2 + 2.58e+01x + 2.14e+02", size=3.5)
print(gg4)

# 5.2. Ratio GPP and Reco

#Compute fit function with the best fit
Harvest_Ratio_Fun<-function(x){0.396598614*(x^0.387995395)*(exp(-0.007407526*x))}
Fire_Ratio_Fun<-function(x){1.187271e-06*x^3 -3.611638e-04*x^2 + 3.212131e-02*x+3.052951e-01}

#Plot data
#Harvest
gg5<-ggplot(Ratio_High_Harvest,aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 6)) +
  facet_wrap(~Disturbance, scales="free_x")+
  stat_function(fun=Harvest_Ratio_Fun, color="red") +
  expand_limits(x = 0, y = 0)+
  geom_hline(yintercept=1, linetype=2, colour="grey", size=0.7)+
  xlab("Year since disturbance") + ylab("GPP/Reco")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.80,1))+
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")+
  ylim(0,2.5)+
  annotate("text", x = 13, y = 2.30, label = "y=0.40x^0.39 X exp(-0.007x)", size=3.5)
print(gg5)

#Wildfire
gg6<-ggplot(Ratio_High_Fire,aes(Stand_Age, values, colour=mean_Uncert)) +
  geom_point(aes(size = Annual_Preci), alpha = 0.5, position = "jitter") +
  scale_size(range = c(3, 8)) +
  facet_wrap(~Disturbance, scales="free_x")+
  stat_function(fun=Fire_Ratio_Fun, color="red") +
  expand_limits(x = 0, y = 0)+
  geom_hline(yintercept=1, linetype=2, colour="grey", size=0.7)+
  xlab("Year since disturbance") + ylab("GPP/Reco")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(low="#FF0000", high = "#00FF33", limits=c(0.80,1))+
  guides(fill=guide_legend(ncol=2)) +
  labs(colour="Fraction of gap-filled data", size="Annual precipitation (mm.y-1)")+
  ylim(0,2.5)+
  annotate("text", x = 40, y = 2.3, label = "y=1.19e-06x^3 -3.61e-04x^2 + 3.21e-02x+3.05e-01", size=3.5)
print(gg6)

#5.3.Plot carbon fluxes and GPP/Reco ratio together
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
Flux_High<-dfAll_Sites[dfAll_Sites$Int_Replacement %in% c("High"),]
GPP_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco_High<-Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
NEE_High<-Flux_High[Flux_High$Type_Flux %in% c("NEE"),]
GPP_High_Harvest<-GPP_High[GPP_High$Disturbance %in% c("Harvest"),]
Reco_High_Harvest<-Reco_High[Reco_High$Disturbance %in% c("Harvest"),]
GPP_High_Fire<-GPP_High[GPP_High$Disturbance %in% c("Wildfire"),]
Reco_High_Fire<-Reco_High[Reco_High$Disturbance %in% c("Wildfire"),]
Ratio_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]
Ratio_High_Harvest<-Ratio_High[Ratio_High$Disturbance %in% c("Harvest"),]
Ratio_High_Fire<-Ratio_High[Ratio_High$Disturbance %in% c("Wildfire"),]

#Compute linear model
fit = lm(values ~ Annual_Preci + Stand_Age +Tair + Tsoil + Rg, data=GPP_High)

m5 = fit
m4 = update(m5, ~ . - Tair)
m3 = update(m4, ~ . - Tsoil)
m2=  update(m3, ~ . - Rg)
m1=  update(m2, ~ . - Stand_Age)
m0=  update(m1, ~ . - Annual_Preci)

af<-anova(m0,m1,m2,m3,m4,m5)
summary(af)
afss <- af$"Sum of Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

#Compute performance of the model
af<-step(fit, direction = c("forward"), steps = 2000)
summary(af)

# 6.2 Principal Components Analysis

#Subset the data
Reco_High_PCA<-data.frame(cbind(Reco_High$values, Reco_High$Annual_Preci, Reco_High$Stand_Age,
                                   Reco_High$Tair, Reco_High$Tsoil,
                                   Reco_High$Rg))

colnames(Reco_High_PCA)<-c("Reco", "Precipitation", "Stand_Age", "Tair", "Tsoil", "Rg")

#Plot the pairwise scatterplots
gg7<-ggpairs(Reco_High_PCA,
        columns = c(1,2,3,4,5,6),
        legends=T,
        lower = list(continuous = "points"),
        diag = list(continuous = "density"),
        axisLabels = "none",
        title = "Correlation matrix - Reco")+ 
  theme_bw(base_size = 14, base_family = "Helvetica")
print(gg7)

# Perform the principal component analysis
arc.pca1<-prcomp(na.omit(Reco_High_PCA), scores=T, cor=T)
summary(arc.pca1, loadings=T)

# create data frame with scores
scores = as.data.frame(arc.pca1$x)

# plot of observations
gg8<-ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "tomato", alpha = 0.8, size = 4)
print(gg8)

# function to create a circle of correlations
circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(na.omit(Reco_High_PCA), arc.pca1$x))

# data frame with arrows coordinates
arrows = data.frame(x1 = c(0, 0, 0, 0, 0, 0), y1 = c(0, 0, 0, 0, 0, 0), x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
gg9<-ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "red") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs", y = "pc2 axis")+
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("PCA-Reco")
print(gg9)

# 7. Compute Random Forest on the annual flux

#7.1 GPP

#Subset the data
GPP_High<-data.frame(cbind(GPP_High$values, GPP_High$Annual_Preci, GPP_High$Stand_Age,
                           GPP_High$Tair, GPP_High$Tsoil,
                           GPP_High$Rg))

colnames(GPP_High)<-c("GPP", "Precipitation", "Stand_Age", "Tair", "Tsoil", "Rg")

# plot panels for each covariate colored by the logical chas variable.
dta <- melt(GPP_High, id.vars=c("GPP"))

ggplot(dta, aes(x=GPP, y=value))+
  geom_point(alpha=.4)+
  geom_rug(data=dta %>% filter(is.na(value)))+
  # labs(y="", x=st.labs["GPP"]) +
  scale_color_brewer(palette="Set2")+
  facet_wrap(~variable, scales="free_y", ncol=2)+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw(base_size = 12, base_family = "Helvetica") 

## Create a basic Random Forest

# Create Training & Test Sets
set.seed(42)
rownames(GPP_High)<-1:nrow(GPP_High)
rows <- sample(x=1:nrow(GPP_High),size=0.7 *  nrow(GPP_High))
train <- GPP_High[rows,]
test <-GPP_High[!rownames(GPP_High) %in% rows,]

#Train a Random Foret 
set.seed(42)
rf_GPP<-randomForest(GPP~., data=train, ntree = 1000, na.action = "na.omit", importance=T)

# Compute correlation
cor_rf_GPP<-(cor(predict(rf_GPP,train), train$GPP, use="pairwise"))^2

## Create a final Random forest
nobs <- nrow(na.omit(train))
crossclass <- sample(5,nobs,TRUE)
crossPredict <- rep(NA,nobs)
for (i in 1:5) {
  indtrain<-which(crossclass!=i)
  indvalidate<-setdiff(1:nobs,indtrain)
  rf_GPP_CV<-randomForest(GPP~., data=train[indtrain,], ntree = 1000, importance = T, na.action = "na.omit")
  crossPredict[indvalidate]<-predict(rf_GPP_CV, train[indvalidate,])
}

#Compute correlation
cor_rf_GPP_CV<-(cor(crossPredict, train$GPP, use='pairwise'))^2

# Plot the OOB errors  and variable importance importance .
imp<-importance(rf_GPP_CV)
vars<-dimnames(imp)[[1]]
imp<-data.frame(vars=vars,imp=as.numeric(imp[,1]))
imp<-imp[order(imp$imp,decreasing=T),]
par(mfrow=c(1,2))
varImpPlot(model0,main='Variable Importance Plot: Base Model')
plot(model0,main='Error vs No. of trees plot: Base Model')

## Model prediction

# Estimate MSE and RMSE
pred<-predict(object=rf_GPP_CV, newdata=test)
actual<-test$GPP
result<-data.frame(actual=actual,predicted=pred)
paste('Function Call: ', rf_GPP_CV$call)
paste('Mean Squared error: ',mean(rf_GPP_CV$mse))
paste('Root Mean Squared error: ',mean(sqrt(rf_GPP_CV$mse)))

#Plot observation vs. prediction
ggplot(result)+
  geom_point(aes(x=actual,y=predicted,color=predicted-actual),alpha=0.7)+
  ggtitle('Plotting Error')+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw(base_size = 12, base_family = "Helvetica") 

#7.2. NEE 
#Subset the data
NEE_High<-data.frame(cbind(NEE_High$values, NEE_High$Annual_Preci, NEE_High$Stand_Age,
                           NEE_High$Tair, NEE_High$Tsoil,
                           NEE_High$Rg))

colnames(NEE_High)<-c("NEE", "Precipitation", "Stand_Age", "Tair", "Tsoil", "Rg")
NEE_High<-na.omit(NEE_High)

# plot panels for each covariate colored by the logical chas variable.
dta <- melt(NEE_High, id.vars=c("NEE"))

ggplot(dta, aes(x=NEE, y=value))+
  geom_point(alpha=.4)+
  geom_rug(data=dta %>% filter(is.na(value)))+
  # labs(y="", x=st.labs["NEE"]) +
  scale_color_brewer(palette="Set2")+
  facet_wrap(~variable, scales="free_y", ncol=2)+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw(base_size = 12, base_family = "Helvetica") 

## Create a basic Random Forest

# Create Training & Test Sets
set.seed(42)
rownames(NEE_High)<-1:nrow(NEE_High)
rows <- sample(x=1:nrow(NEE_High),size=0.7 *  nrow(NEE_High))
train <- NEE_High[rows,]
test <-NEE_High[!rownames(NEE_High) %in% rows,]

#Train a Random Foret 
set.seed(42)
rf_NEE<-randomForest(NEE~., data=train, ntree = 1000, na.action = "na.omit", importance=T)

# Compute correlation
cor_rf_NEE<-(cor(predict(rf_NEE,train), train$NEE, use="pairwise"))^2

## Create a final Random forest
nobs <- nrow(na.omit(train))
crossclass <- sample(5,nobs,TRUE)
crossPredict <- rep(NA,nobs)
for (i in 1:5) {
  indtrain<-which(crossclass!=i)
  indvalidate<-setdiff(1:nobs,indtrain)
  rf_NEE_CV<-randomForest(NEE~., data=train[indtrain,], ntree = 1000, importance = T, na.action = "na.omit")
  crossPredict[indvalidate]<-predict(rf_NEE_CV, train[indvalidate,])
}

#Compute correlation
cor_rf_NEE_CV<-(cor(crossPredict, train$NEE, use='pairwise'))^2

# Plot the OOB errors  and variable importance importance .
imp<-importance(rf_NEE_CV)
vars<-dimnames(imp)[[1]]
imp<-data.frame(vars=vars,imp=as.numeric(imp[,1]))
imp<-imp[order(imp$imp,decreasing=T),]
par(mfrow=c(1,2))
varImpPlot(model0,main='Variable Importance Plot: Base Model')
plot(model0,main='Error vs No. of trees plot: Base Model')

## Model prediction

# Estimate MSE and RMSE
pred<-predict(object=rf_NEE_CV, newdata=test)
actual<-test$NEE
result<-data.frame(actual=actual,predicted=pred)
paste('Function Call: ', rf_NEE_CV$call)
paste('Mean Squared error: ',mean(rf_NEE_CV$mse))
paste('Root Mean Squared error: ',mean(sqrt(rf_NEE_CV$mse)))

#Plot observation vs. prediction
ggplot(result)+
  geom_point(aes(x=actual,y=predicted,color=predicted-actual),alpha=0.7)+
  ggtitle('Plotting Error')+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw(base_size = 12, base_family = "Helvetica") 

#7. Map with the location of the fluxnet sites

# Import site information
Site_Info<-read.csv("Input/Potential_Sites.csv", header = TRUE)
Site_Info<- Site_Info[-c(28,36,37,49,50,51,53,63),]
Site_Info<-Site_Info[Site_Info$Intensity_Replacement %in% c("High"),]

#7.1 Land shapefile

# read land shapefile
wmap <- readOGR(dsn="Input/Flux_Site_Map", layer="ne_110m_land")

#Change CRS
wmap_robinson <- spTransform(wmap, CRS("+proj=wintri"))

# convert Spdf to a dataframe
wmap_robinson_df <- fortify(wmap_robinson)

# 7.2 Graticule and bounding box (longlat)

#Read shapefile
grat <- readOGR("Input/Flux_Site_Map", layer="ne_110m_graticules_15") 

#Change CRS
grat_robinson <- spTransform(grat, CRS("+proj=wintri"))

# convert Spdf to a dataframe
grat_robinson_df <- fortify(grat_robinson)

# 7.3 Bouding box

#Read shapefile
bbox <- readOGR("Input/Flux_Site_Map", layer="ne_110m_wgs84_bounding_box")

#Change CRS
bbox_robinson <- spTransform(bbox, CRS("+proj=wintri"))

# convert Spdf to a dataframe
bbox_robinson_df <- fortify(bbox_robinson)

# 7.4. Import Fluxnet information
Site <- SpatialPointsDataFrame(coords=Site_Info[,c("y","x")],data=data.frame(Site_Info),proj4string=CRS(proj="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
Site_Robinson <- spTransform(Site, CRS("+proj=wintri"))
Site_Robinson_df <- as(Site_Robinson, "data.frame")

# 7.5 Plot fluxnet site

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


gg8<-ggplot(bbox_robinson_df, aes(long,lat, group=group)) + 
  geom_polygon(fill="white") +
  geom_polygon(data=wmap_robinson_df, aes(long,lat, group=group, fill=hole)) + 
  geom_path(data=grat_robinson_df, aes(long, lat, group=group, fill=NULL), linetype="dashed", color="grey50") +
  geom_point(data=Site_Robinson_df, aes(x=y.1,y=x.1, group= NULL),  shape=20, color="red", size=3)+
  # labs(title="Fluxnet sites location") + 
  coord_equal() + 
  theme_opts +
  scale_fill_manual(values=c("grey", "white"), guide="none") # change colors & remove legend
print(gg8)