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
# install.packages("gamm4")
# install.packages("tidyr")
# install.packages("manipulate")
# install.packages("fitdistrplus")
# install.packages("evd")
# install.packages("flexmix")


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
library(tidyr)
library(manipulate)
library(MASS)
library(fitdistrplus)
library(evd)
library(flexmix)

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
                    sd_Reco= sd(Reco, na.rm=T),
                    mean_Uncert=mean(NEE_fqcOK, na.rm=T)))

#Add type of ecosystem, type of climate, type of disturbance and year after disturbance for each sites
for (i in seq_along(Sum_Sd_Flux)){
  Sum_Sd_Flux[[i]]$Ecosytem<- Site_Date$Ecosystem[i]
  Sum_Sd_Flux[[i]]$Climate<- Site_Date$Climate[i]
  Sum_Sd_Flux[[i]]$Disturbance<- Site_Date$Type_Disturbance[i]
  Sum_Sd_Flux[[i]]$Year_Disturbance<- Sum_Sd_Flux[[i]]$year- year(Site_Date$Measure_Date[i])
  Sum_Sd_Flux[[i]]$Species<- Site_Date$Species[i]
  Sum_Sd_Flux[[i]]$GPP_ER<- Sum_Sd_Flux[[i]]$sum_GPP / Sum_Sd_Flux[[i]]$sum_TER
  Sum_Sd_Flux[[i]]$Site_ID<- Site_Date$ID[i]
}

#Combine the flux sites in one dataframe
dfAll_Sites<- do.call("rbind", Sum_Sd_Flux)
dfAll_Sites<-dfAll_Sites[-c(88),]# the measurements seems to be an outlier

# Restructure dataframe
dfAll_Sites<-gather(dfAll_Sites, Type_Flux, values, -year, -Species, -Ecosytem, -mean_Uncert, -Climate, -Disturbance, -Year_Disturbance, -Site_ID)

#Reoder column
dfAll_Sites<- dfAll_Sites[c("Site_ID", "year", "Type_Flux", "values", "mean_Uncert", "Year_Disturbance", "Disturbance", "Climate", "Ecosytem", "Species")]

# Reclassify climate classification
dfAll_Sites$Climate<-ifelse((dfAll_Sites$Climate =="Af" | dfAll_Sites$Climate =="Am"), "Tropical",
                    ifelse((dfAll_Sites$Climate =="Cfa" | dfAll_Sites$Climate =="Cfb" | dfAll_Sites$Climate =="Cfc" | dfAll_Sites$Climate =="Csa" | dfAll_Sites$Climate =="Csb"), "Temperate",
                           ifelse((dfAll_Sites$Climate =="Dfb" | dfAll_Sites$Climate =="Dfc" | dfAll_Sites$Climate =="Dsc"), "Continental", NA)))

# 4. Plot all sites together

#Fitting non-linear model to ecosystem response - Amiro et al. (2010)
Amiro_Fun <- function(x){1.23*(1-exp(-0.224*x))}

#Fitting non-linear model to ecosystem response - Gamma function - Tang et al. (2014)
Gamma_Fun <- function(x){975.14*(x^0.050)*(exp(-0.00295*x))}

#Fitting non-linear model to ecosystem response - Ricker function
Ricker_Fun<- function(x){500*x*(exp(-0.3*x/2))}

#Fitting non-linear model to ecosystem response - Spline function
smooth.spline2 <- function(formula, data, ...) { 
  mat <- model.frame(formula, data) 
  smooth.spline(mat[, 2], mat[, 1]) 
} 
predictdf.smooth.spline <- function(model, xseq, se, level) { 
  pred <- predict(model, xseq) 
  data.frame(x = xseq, y = pred$y) 
}

#Fitting non-linear model to ecosystem response - Beta function
fit <- fgev(Plot_Flux$values)
hist(Plot_Flux$values,prob=T,col="gray", xlab="", ylab="", main="")
legend("topright", 
       legend=c("density", "fgev", "flexmix"), 
       fill=c("darkgreen", "blue", "darkred")
)
xval <- seq(from=0, to=max(Plot_Flux$values))

# density
fit1 <- density(Plot_Flux$values)
lines(fit1, col="darkgreen", lwd=2)

# generalized extreme value distribution
fit2 <- fgev(Plot_Flux$values)
param2 <- fit2$estimate
loc <- param2[["loc"]]
scal <- param2[["scale"]]
shape <- param2[["shape"]]
lines(xval, dgev(xval, loc=loc, scale=scal, shape=shape), col="blue", lwd=2)

# mixture of two Gamma distributions
# http://r.789695.n4.nabble.com/Gamma-mixture-models-with-flexmix-tt3328929.html#none
fit3 <- flexmix(values~1, data=subset(Plot_Flux, values>0), k=2, 
                model = list(FLXMRglm(family = "Gamma"), FLXMRglm(family = "Gamma"))
)
param3 <- parameters(fit3)[[1]] # don't know why this is a list
interc <- param3[1,]
shape <- param3[2,]
lambda <- prior(fit3)
yval <- lambda[[1]]*dgamma(xval, shape=shape[[1]], rate=interc[[1]]*shape[[1]]) + 
  lambda[[2]]*dgamma(xval, shape=shape[[2]], rate=interc[[2]]*shape[[2]])
lines(xval, yval, col="darkred", lwd=2)

#Compute error bars
limits_NEE <- aes(ymax = mean_NEE + sd_NEE, ymin=mean_NEE - sd_NEE)
limits_GPP <- aes(ymax = mean_GPP + sd_GPP, ymin=mean_GPP - sd_GPP)
limits_TER <- aes(ymax = mean_TER + sd_Reco, ymin=mean_TER - sd_Reco)

x<-Plot_Flux[, c("values", "Year_Disturbance")] 
dput(x)


# 4.1. Partition flux data per variable

#Subset data
Plot_Flux<- dfAll_Sites[dfAll_Sites$Type_Flux %in% c("sum_NEE", "sum_GPP", "sum_TER"),]
Plot_Flux<- Plot_Flux[Plot_Flux$Disturbance %in% c("Fire"),]
Plot_Flux<- Plot_Flux[!Plot_Flux$values ==0,]

gg4<-ggplot(Plot_Flux, aes(Year_Disturbance, values, colour=mean_Uncert)) +
  geom_point(size=3) +
  facet_grid(Type_Flux~Disturbance, scales = "free")+
  # stat_smooth(fun="lm", color="red", formula = y ~ x + I(x^2), size = 1)+
  # geom_smooth(method = "smooth.spline2", se= F)+
  stat_function(fun=Gamma_Fun, color="red")+
  # geom_errorbar(limits_NEE, width=0.07, linetype=6)+
  xlab("Year since Disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(name="Uncertainty", low="#FF0000", high = "#00FF33")

# 4.2. Ratio GPP and Reco

#Subset data
Plot_Ratio<- dfAll_Sites[dfAll_Sites$Type_Flux %in% c("GPP_ER"),]
Plot_Ratio<- Plot_Ratio[Plot_Ratio$Species %in% c("Douglas pine", "Black spruce", "Jack pine"),]

gg7<-ggplot(Plot_Ratio,aes(Year_Disturbance, values, shape=Disturbance, colour=mean_Uncert)) +
  geom_point (size=3.5)+
  facet_wrap(~Disturbance, ncol=1, scales="free")+
  stat_function(fun=Amiro_Fun, color="red") +
  expand_limits(x = 0, y = 0)+
  geom_hline(yintercept=1, linetype=2, colour="grey", size=0.7)+
  # geom_errorbar(limits_TER, width=0.07, linetype=6)+
  xlab("Year since disturbance") + ylab("GPP/Reco")+ 
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_colour_gradient(name="Uncertainty", low="#FF0000", high = "#00FF33")

#Export plot
ggsave(gg7, filename = 'Latex/Figures/All_Sites/Dist_GPP_TER.eps', width = 14, height = 8)

#5 Explain variabilty of the fluxes

#5.1. Anova analysis

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
