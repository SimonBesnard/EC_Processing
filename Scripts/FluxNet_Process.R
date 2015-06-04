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


library(RNetCDF)
library (ggplot2)
library(scales)
library(lubridate)
library(gridExtra)
library(REddyProc)
library (dplyr)

#Load the nc file
AU_Tum_Site<-"BR-Sa2.2001.2002.daily.nc"

#Load EC variables 
df <- fLoadFluxNCIntoDataframe(VarList.V.s=c('year', 'NEE_f','GPP_f','Reco','NEE_fqcok'), 
                               FileName.s=AU_Tum_Site, 
                               Dir.s=file.path(path, "Fluxnet_Data/Data_4disturbed/"), 
                               NcPackage.s = 'RNetCDF')
# Remove data before disturbance
df<-df[-c(1:730),]

# Compute mean annual +/- sd
Mean_Sd_Flux<-ddply(df, "year",
  summarize,
  mean_NEE= mean(NEE_f, na.rm=T),
  mean_GPP= mean(GPP_f, na.rm=T),
  mean_TER= mean(Reco,na.rm=T),
  sd_NEE= sd(NEE_f, na.rm=T),
  sd_GPP= sd(GPP_f, na.rm=T),
  sd_Reco= sd(Reco, na.rm=T))

#Plot data raw flux data
# NEE
p1<-ggplot(df, aes(df$DateTime, df$NEE_f)) + geom_line(size=0.2) + 
  xlab("") + ylab("NEE (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#GPP
p3<-ggplot(df, aes(df$DateTime, df$GPP_f)) + geom_line(size=0.2) + 
  xlab("") + ylab("GPP (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#TER
p4<-ggplot(df, aes(df$DateTime, df$Reco)) + geom_line(size=0.2) + 
  xlab("") + ylab("Reco (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p1, filename = 'Latex/Figures/BR_Sa2/BR_SA2_NEE.eps', width = 12, height = 4)

#Plot data mean/sd flux data

#NEE
limits_NEE <- aes(ymax = mean_NEE + sd_NEE, ymin=mean_NEE - sd_NEE)
limits_GPP <- aes(ymax = mean_GPP + sd_GPP, ymin=mean_GPP - sd_GPP)
limits_TER <- aes(ymax = mean_TER + sd_Reco, ymin=mean_TER - sd_Reco)

gg1<-ggplot(Mean_Sd_Flux, aes(Mean_Sd_Flux$year, Mean_Sd_Flux$mean_NEE)) + geom_point(size=3) +
  geom_path()+
  geom_errorbar(limits_NEE, width=0.25)+
  xlab("") + ylab("Annual Mean-NEE (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg2<-ggplot(Mean_Sd_Flux, aes(Mean_Sd_Flux$year, Mean_Sd_Flux$mean_GPP)) + geom_point(size=3) +
  geom_path()+
  geom_errorbar(limits_GPP, width=0.25)+
  xlab("") + ylab("Annual Mean-GPP (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg3<-ggplot(Mean_Sd_Flux, aes(Mean_Sd_Flux$year, Mean_Sd_Flux$mean_TER)) + geom_point(size=3) +
  geom_path()+
  geom_errorbar(limits_TER, width=0.25)+
  xlab("") + ylab("Annual Mean-TER (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(gg3, filename = 'Latex/Figures/Au_Tum_TER_Mean.eps', width = 9, height = 4)

