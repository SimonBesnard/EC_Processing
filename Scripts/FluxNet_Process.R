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
install.packages("gamm4")

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
library (RColorBrewer)
library(gamm4)

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

# Remove uncomplete year after disturbance
Fluxnet_Site[[1]]$DateTime

#3. Analyse inter-annual variability of the fluxex

# Compute sum annual +/- sd
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
                    sd_Reco= sd(Reco, na.rm=T)))

#Add type of ecosystem, type of climate, type of disturbance and year after disturbance for each sites
for (i in seq_along(Sum_Sd_Flux)){
  Sum_Sd_Flux[[i]]$Ecosytem<- Site_Date$Ecosystem[i]
  Sum_Sd_Flux[[i]]$Climate<- Site_Date$Climate[i]
  Sum_Sd_Flux[[i]]$Disturbance<- Site_Date$Type_Disturbance[i]
  Sum_Sd_Flux[[i]]$Year_Disturbance<- Sum_Sd_Flux[[i]]$year- year(Site_Date$Measure_Date[i])
  Sum_Sd_Flux[[i]]$GPP_ER<- Sum_Sd_Flux[[i]]$sum_GPP / Sum_Sd_Flux[[i]]$sum_TER
  Sum_Sd_Flux[[i]]$Site_ID<- Site_Date$ID[i]
}

#Combine the flux sites in one dataframe
dfAll_Sites<- do.call("rbind", Sum_Sd_Flux)
dfAll_Sites<-dfAll_Sites[-c(88),]# the measurements seems to be an outlier

#Fitting non-linear model to ecosystem response
fit.mean <- function(x){1.23*(1-exp(-0.224*x))}

# 4. Plot all sites together

#Keep harvest and fire disturbance
dfAll_Sites<- dfAll_Sites[!dfAll_Sites$Disturbance %in% c("Insect_Outbreaks", "Thinning"),]

dfAll_Sites<-read.csv("Output/dfAll_Sites.csv", sep=",")

#Compute error
limits_NEE <- aes(ymax = mean_NEE + sd_NEE, ymin=mean_NEE - sd_NEE)
limits_GPP <- aes(ymax = mean_GPP + sd_GPP, ymin=mean_GPP - sd_GPP)
limits_TER <- aes(ymax = mean_TER + sd_Reco, ymin=mean_TER - sd_Reco)

#4.1. Without partitioning variable
colourCount = length(unique(dfAll_Sites$Site_ID))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

#NEE
gg1<-ggplot(dfAll_Sites, aes(Year_Disturbance, sum_NEE, shape=Disturbance)) +
  geom_point(size=3) +
  # geom_path()+
  # geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("Year since Disturbance") + ylab("NEE (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  # scale_colour_manual(name="Site ID", values=getPalette(colourCount))+
  guides(colour = guide_legend(title.position ="top", title.hjust =0.5, override.aes = list(size=3), ncol=2))

#GPP
gg2<-ggplot(dfAll_Sites, aes(Year_Disturbance, sum_GPP, shape=Disturbance)) +
  geom_point(size=3) +
  # geom_path()+
  # geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("Year since Disturbance") + ylab("GPP (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  # scale_colour_manual(name="Site ID", values=getPalette(colourCount))+
  guides(colour = guide_legend(title.position ="top", title.hjust =0.5, override.aes = list(size=3), ncol=2))

#TER
gg3<-ggplot(dfAll_Sites, aes(Year_Disturbance, sum_TER, shape=Disturbance, colour=Site_ID)) +
  geom_point(size=3) +
  # geom_path()+
  # geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("Year since Disturbance") + ylab("TER (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_manual(name="Site ID", values=getPalette(colourCount))+
  guides(colour = guide_legend(title.position ="top", title.hjust =0.5, override.aes = list(size=3), ncol=2))

# 4.2. Partition flux data per variable

#NEE
gg4<-ggplot(dfAll_Sites, aes(Year_Disturbance, Value, shape=Disturbance)) +
  facet_grid(Flux_Type~Disturbance)+
  geom_point(size=3, shape=3) +
  geom_smooth()+
  # geom_path()+
  # geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("Year since Disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_manual(name="Site ID", values=getPalette(colourCount))+
  guides(colour = guide_legend(title.position ="top", title.hjust =0.5, override.aes = list(size=3), ncol=2))

#GPP
gg5<-ggplot(dfAll_Sites, aes(Year_Disturbance, sum_GPP)) +
  facet_wrap(~Disturbance, ncol = 1)+
  geom_point(size=3, shape=3) +
  geom_smooth()+
  # geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("Year since Disturbance") + ylab("GPP (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  # scale_colour_manual(name="Site ID", values=getPalette(colourCount))+
  guides(colour = guide_legend(title.position ="top", title.hjust =0.5, override.aes = list(size=3), ncol=2))

#TER
gg6<-ggplot(dfAll_Sites, aes(Year_Disturbance, sum_TER)) +
  facet_wrap(~Disturbance, ncol = 1)+
  geom_point(size=3, shape=3) +
  geom_smooth()+
  # geom_path()+
  # geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("Year since Disturbance") + ylab("TER (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  # scale_colour_manual(name="Site ID", values=getPalette(colourCount))+
  guides(colour = guide_legend(title.position ="top", title.hjust =0.5, override.aes = list(size=3), ncol=2))

# 4.3. Ratio GPP and Reco
gg7<-ggplot(dfAll_Sites,aes(Year_Disturbance, GPP_ER)) +
  geom_point (shape=3, size=3)+
  facet_wrap(~Disturbance, ncol=1, scales="free_y")+
  stat_function(fun=fit.mean, color="red") +
  expand_limits(x = 0, y = 0)+
  # scale_y_continuous(breaks = round(seq(min(dfAll_Sites$GPP_ER), max(dfAll_Sites$GPP_ER), by = 0.5),1))+
  geom_hline(yintercept=1, linetype=2, colour="grey", size=0.7)+
  # geom_errorbar(limits_TER, width=0.07, linetype=6)+
  xlab("Year since disturbance") + ylab("GPP/Reco")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Export plot
ggsave(gg7, filename = 'Latex/Figures/All_Sites/Dist_GPP_TER.eps', width = 14, height = 8)

#4 Explain variabilty of the fluxes

#4.1. Fluxes

#Compute linear model
fit = lm(mean_GPP ~ Disturbance + Climate + Ecosytem, data=dfAll_Sites)

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
af<-step(fit, direction = c("both"), steps = 2000)
summary(af)
