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

# 1. 2. Compute ANOVA for type of disturbance

# 1.2.1 NEP

#Boxplot and interactions plots
attach(NEP_High)
boxplot(values~Disturbance)
interaction.plot(Stand_Age, Disturbance, values)

#Perform Anova
lm1 <- lm(values~ Disturbance, data=NEP_High) 
Anova_NEP<-aov(lm1)
summary(Anova_NEP)

# 1.2.2 GPP

#Boxplot and interactions plots
attach(GPP_High)
boxplot(values~Disturbance)
interaction.plot(Stand_Age, Disturbance, values)

#Perform Anova
lm1 <- lm(values~ Disturbance*Stand_Age, data=GPP_High) 
Anova_GPP<-aov(lm1)
summary(Anova_GPP)

# 1.2.3 Respiration

#Boxplot and interactions plots
attach(Reco_High)
boxplot(values~Disturbance)
interaction.plot(Stand_Age, Disturbance, values)

#Perform Anova
lm1 <- lm(values~ Disturbance*Stand_Age, data=Reco_High) 
Anova_Reco<-aov(lm1)
summary(Anova_Reco)

# 1.2.4. Ratio GPP-Reco

#Boxplot and interactions plots
attach(Ratio_High)
boxplot(values~Disturbance)
interaction.plot(Stand_Age, Disturbance, values)

#Perform Anova
lm1 <- lm(values~ Disturbance*Stand_Age, data=Ratio_High) 
Anova_Ratio<-aov(lm1)
summary(Anova_Ratio)

# 1.3 Compute ANOVA for climate variables

# 1.3.1 NEP

#Perform Anova
lm1 <- lm(values~ Annual_Preci + Tair + Stand_Age, data=NEP_High) 
Anova_NEP<-aov(lm1)
summary(Anova_NEP)

# 1.3.2 GPP

#Perform Anova
lm1 <- lm(values~ Disturbance*Stand_Age, data=GPP_High) 
Anova_GPP<-aov(lm1)
summary(Anova_GPP)

# 1.3.3 Respiration

#Perform Anova
lm1 <- lm(values~ Disturbance*Stand_Age, data=Reco_High) 
Anova_Reco<-aov(lm1)
summary(Anova_Reco)

# 1.3.4. Ratio GPP-Reco

#Perform Anova
nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = Ratio_High, 
    start = list(A = 100, B = 0.170, k = -0.00295))

f1<- function(x) {0.384394*(x^0.341429)*(exp(-0.004749 *x))}
lm1 <- lm(values~ Annual_Preci + Tair + Stand_Age, data=Ratio_High)
Anova_Ratio<-aov(lm1)
summary(Anova_Ratio)

df<-Ratio_High
df<-df[c("values", "Annual_Preci", "Tair", "Stand_Age")]
colnames(df)<- c("Y", "P", "TA", "A")
df<-df[c(1:10),]

dput(df)


