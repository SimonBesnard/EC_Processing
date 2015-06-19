## Script to process EC data
## Author: Simon Besnard
## 15.05.2015
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

library(RNetCDF)
library (ggplot2)
library(scales)
library(lubridate)
library(gridExtra)
library(REddyProc)
library (dplyr)
library (plyr)
require(testthat)
library(gtools)
library(FactoMineR)
library (leaps)

#1. Create a df for all fluxnet sites dataframe

#Open all files
dir <- file.path(path, 'Fluxnet_Data/Data_4disturbed')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = try(fLoadFluxNCIntoDataframe(VarList.V.s=c('NEE_f','GPP_f','Reco','NEE_fqcok'),
                                          FileName.s=list[i], NcPackage.s = 'RNetCDF'))
}

# Remove the broken flux net sites from the list of dataframe
Fluxnet_Site<-Fluxnet_Site[sapply(Fluxnet_Site, function(x) !inherits(x, "try-error"))]

# Remove data before disturbance from the list of dataframe
Site_Date<-read.csv("Input/Potential_Sites.csv", header = TRUE)
Site_Date$Measure_Date<- as.Date(Site_Date$Measure_Date)
for (i in seq_along(Fluxnet_Site)){
  Fluxnet_Site[[i]]$DateTime=as.Date(Fluxnet_Site[[i]]$DateTime)
}

for (i in seq_along(Fluxnet_Site)){
Fluxnet_Site[[i]]<-with(Fluxnet_Site[[i]], Fluxnet_Site[[i]][(Fluxnet_Site[[i]]$DateTime > Site_Date$Measure_Date[i]),])
}

#3. Analyse inter-annual variability of the fluxex

# Compute mean annual +/- sd
Mean_Sd_Flux<-lapply(Fluxnet_Site, function (x) ddply(x, "year",
                    summarize,
                    mean_NEE= mean(NEE_f, na.rm=T),
                    mean_GPP= mean(GPP_f, na.rm=T),
                    mean_TER= mean(Reco,na.rm=T),
                    sd_NEE= sd(NEE_f, na.rm=T),
                    sd_GPP= sd(GPP_f, na.rm=T),
                    sd_Reco= sd(Reco, na.rm=T)))

#Add type of ecosystem, type of climate, type of disturbance and year after disturbance for each sites
for (i in seq_along(Mean_Sd_Flux)){
  Mean_Sd_Flux[[i]]$Ecosytem<- Site_Date$Ecosystem[i]
  Mean_Sd_Flux[[i]]$Climate<- Site_Date$Climate[i]
  Mean_Sd_Flux[[i]]$Disturbance<- Site_Date$Type_Disturbance[i]
  Mean_Sd_Flux[[i]]$Year_Disturbance<- Mean_Sd_Flux[[i]]$year- year(Site_Date$Measure_Date[i])
}

#Combine the flux sites in one dataframe
dfAll_Sites<- do.call("rbind", Mean_Sd_Flux)
dfAll_Sites<-dfAll_Sites[-c(53),]# the measurements seems to be an outlier

#Plot data mean/sd flux data
limits_NEE <- aes(ymax = mean_NEE + sd_NEE, ymin=mean_NEE - sd_NEE)
limits_GPP <- aes(ymax = mean_GPP + sd_GPP, ymin=mean_GPP - sd_GPP)
limits_TER <- aes(ymax = mean_TER + sd_Reco, ymin=mean_TER - sd_Reco)

#NEE
gg1<-ggplot(dfAll_Sites, aes(dfAll_Sites$Year_Disturbance, dfAll_Sites$mean_NEE)) +
  facet_wrap(~Disturbance, ncol=1)+
  geom_point(size=2, shape=3) +
  # geom_path()+
  # geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("Year since Disturbance") + ylab("Annual Mean-NEE (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#GPP
gg2<-ggplot(dfAll_Sites, aes(dfAll_Sites$Year_Disturbance, dfAll_Sites$mean_GPP)) +
  facet_wrap(~Disturbance, ncol=1)+
  geom_point(size=2, shape=3) +
  # geom_path()+
  # geom_errorbar(limits_GPP, width=0.07, linetype=6)+
  xlab("Year since disturbance") + ylab("Annual Mean-GPP (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#TER
gg3<-ggplot(dfAll_Sites, aes(dfAll_Sites$Year_Disturbance, dfAll_Sites$mean_TER)) +
  geom_point(size=2, shape=3) +
  facet_wrap(~Disturbance, ncol=1)+
  # geom_path()+
  # geom_errorbar(limits_TER, width=0.07, linetype=6)+
  xlab("Year since disturbance") + ylab("Annual Mean-TER (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave(gg3, filename = 'Latex/Figures/Dist_TER_Mean.eps', width = 10, height = 8)

#4 Explain variabilty of the fluxes

#4.1. Fluxes

#NEE
fit = lm(mean_NEE ~ Disturbance + Climate + Ecosytem, data=dfAll_Sites)

#Compute Anova
summary(fit)
af<-anova(fit)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

m3 = fit
m2 = update(m3, ~ . - Ecosytem)
m1 = update(m2, ~ . - Climate)
af<-anova(m1,m2,m3)
afss <- af$"Sum of Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

#Compute performance of the model
step(fit, direction = c("both"), steps = 2000)

#GPP
fit = lm(mean_GPP ~ Climate + Ecosytem + Disturbance, data=dfAll_Sites)
summary(fit)

m3 = fit
m2 = update(m3, ~ . - Disturbance)
m1 = update(m2, ~ . - Ecosytem)

anova(m1,m2,m3)

#TER
fit = lm(mean_TER ~ Climate + Ecosytem + Disturbance, data=dfAll_Sites)
summary(fit)

m3 = fit
m2 = update(m3, ~ . - Disturbance)
m1 = update(m2, ~ . - Ecosytem)

anova(m1,m2,m3)
