## Script to compute statistical analysis annual carbon flux
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library (ggplot2)
library(scales)
library (dplyr)
library (plyr)
library(tidyr)
library(gridExtra)
library (car)
library (lme4)

#1 Explain variabilty of the fluxes

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")

#1.1. Subset dataset
Flux_High<-dfAll_Sites[dfAll_Sites$Int_Replacement %in% c("High"),]

#Subset type of flux

#NEP
NEP_High<-Flux_High[Flux_High$Type_Flux %in% c("NEP"),]

#GPP
GPP_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP"),]

#Reco
Reco_High<-Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]

#Ratio GPP-Reco
Ratio_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]

# 1.2 Compute ANOVA and forward stepwise regresion for climate variables type of disturbance and stand age

# 1.2.1 NEP

#Perform Anova
lm1 <- lm(values~ Annual_Preci + Tair + Stand_Age + Disturbance, data=NEP_High) 
Anova_NEP<-aov(lm1)
summary(Anova_NEP)

# Forward stepwise regresion 
lm1 <- lm(values~ Stand_Age+ Tair + Annual_Preci + Disturbance, data=NEP_High) 
fwd.NEP<-step(lm1, direction="both")
summary(fwd.NEP)
summary(fwd.NEP)$r.squared

# 1.2.2 GPP

#Perform Anova
lm1 <- lm(values~ Annual_Preci + Tair + Stand_Age + Disturbance, data=GPP_High) 
Anova_GPP<-aov(lm1)
summary(Anova_GPP)

# Forward stepwise regresion 
lm1 <- lm(values~ Stand_Age+ Tair + Annual_Preci + Disturbance, data=GPP_High) 
fwd.GPP<-step(lm1, direction="both")
summary(fwd.GPP)
summary(fwd.GPP)$r.squared

# 1.2.3 Respiration

#Perform Anova
lm1 <- lm(values~ Annual_Preci + Tair + Stand_Age + Disturbance, data=Reco_High) 
Anova_Reco<-aov(lm1)
summary(Anova_Reco)

# Forward stepwise regresion 
lm1 <- lm(values~ Stand_Age+ Tair + Annual_Preci + Disturbance, data=Reco_High) 
fwd.Reco<-step(lm1, direction="both")
summary(fwd.Reco)
summary(fwd.Reco)$r.squared

# 1.2.4. Ratio GPP-Reco

#Perform Anova
lm1 <- lm(values~ Annual_Preci + Tair + Stand_Age + Disturbance, data=Ratio_High) 
Anova_Ratio<-aov(lm1)
summary(Anova_Ratio)

# Forward stepwise regresion 
lm1 <- lm(values~ Stand_Age+ Tair + Annual_Preci + Disturbance, data=Ratio_High) 
fwd.Ratio<-step(lm1, direction="both")
summary(fwd.Ratio)
summary(fwd.Ratio)$r.squared


