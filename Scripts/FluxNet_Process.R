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

#1. Create a df for all fluxnet sites dataframe

#Open all files
dir <- file.path(path, 'Fluxnet_Data/Data_4disturbed')
list <- list.files(dir, pattern=glob2rx('*.nc'), full.names=TRUE)

#Loop over the list of files
Fluxnet_Site <- list()
for(i in seq_along(list)) {
  Fluxnet_Site[[i]] = try(fLoadFluxNCIntoDataframe(VarList.V.s=c('year', 'NEE_f','GPP_f','Reco','NEE_fqcok'),
                                          FileName.s=list[i], NcPackage.s = 'RNetCDF'))
}

#function that adds a part to the colum names of a data frame for identification,txtn = identifier POS=position to start
source('Function/dfColRename.R')
dfAll_Sites_NEE=smartbind(Au_Tum_df$NEE_f, BR_Sa3_df$NEE_f, BR_Sa2_df$NEE_f, BR_Ji1_df$NEE_f, CA_Ca2_df$NEE_f, CA_Ca1_df$NEE_f, fill=NA)

# Remove data before disturbance
# df<-df[-c(1:730),]

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
p2<-ggplot(df, aes(df$DateTime, df$GPP_f)) + geom_line(size=0.2) + 
  xlab("") + ylab("GPP (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#TER
p3<-ggplot(df, aes(df$DateTime, df$Reco)) + geom_line(size=0.2) + 
  xlab("") + ylab("Reco (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p3, filename = 'Latex/Figures/BR_Sa3/BR_Sa3_TER.png', width = 12, height = 4)

#Plot data mean/sd flux data

#NEE
limits_NEE <- aes(ymax = mean_NEE + sd_NEE, ymin=mean_NEE - sd_NEE)
limits_GPP <- aes(ymax = mean_GPP + sd_GPP, ymin=mean_GPP - sd_GPP)
limits_TER <- aes(ymax = mean_TER + sd_Reco, ymin=mean_TER - sd_Reco)

gg1<-ggplot(Mean_Sd_Flux, aes(Mean_Sd_Flux$year, Mean_Sd_Flux$mean_NEE)) + geom_point(size=3) +
  geom_path()+
  geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("") + ylab("Annual Mean-NEE (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg2<-ggplot(Mean_Sd_Flux, aes(Mean_Sd_Flux$year, Mean_Sd_Flux$mean_GPP)) + geom_point(size=3) +
  geom_path()+
  geom_errorbar(limits_GPP, width=0.07, linetype=6)+
  xlab("") + ylab("Annual Mean-GPP (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

gg3<-ggplot(Mean_Sd_Flux, aes(Mean_Sd_Flux$year, Mean_Sd_Flux$mean_TER)) + geom_point(size=3) +
  geom_path()+
  geom_errorbar(limits_TER, width=0.07, linetype=6)+
  xlab("") + ylab("Annual Mean-TER (g.m-2.day-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(gg3, filename = 'Latex/Figures/BR_Sa2/BR_Sa2_TER_Mean.png', width = 10, height = 6)


