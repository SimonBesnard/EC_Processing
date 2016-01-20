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
library (minpack.lm)
library (randomForest)
library (reshape)
library (hydroGOF)
library (relaimpo)

#1 Explain variabilty of the fluxes using linear regression analysis

# 1.1 Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
NEP<-readRDS("Output/NEP.rds")
NEP_Mean_Site<-readRDS("Output/NEP_Mean_Site.rds")
GPP<-readRDS("Output/GPP.rds")
GPP_Mean_Site<-readRDS("Output/GPP_Mean_Site.rds")
Reco<-readRDS("Output/Reco.rds")
Reco_Mean_Site<-readRDS("Output/Reco_Mean_Site.rds")
Ratio_GPP_Reco<-readRDS("Output/Ratio_GPP_Reco.rds")
Ratio_GPP_Reco_Mean_Site<-readRDS("Output/Ratio_GPP_Reco_Mean_Site.rds")
Ratio_NEP_GPP<-readRDS("Output/Ratio_NEP_GPP.rds")
Ratio_NEP_GPP_Mean_Site<- readRDS("Output/Ratio_NEP_GPP_Mean_Site.rds")
Ratio_NEP_GPPmax<- readRDS("Output/Ratio_NEP_GPPmax.rds")
Ratio_NEP_GPPmax_Mean_Site<- readRDS("Output/Ratio_NEP_GPPmax_Mean_Site.rds")

# 1.2 Analysis for NEP

## All years per site
#-----------------------------------------------------------------
#Compute transform function

# Age
Fun_NEP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = NEP,
               start = list(A= 2.379e+02, B= -3.295e-03, C=-7.879e+02, D=-1.519e-01), control = list(maxiter = 500))
coef(Fun_NEP)
f_Age_NEP<- function (x) {3.092250e+02*(exp(-4.718321e-03*x)) -1.037645e+03*(exp(-1.752753e-01*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_MAT.R")
stat_Temp(NEP)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = NEP, 
                start = list(A=-0.02505, B=-0.11885, C= 28.03883, D=4.61699), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_NEP<- function (x) {-0.06796338*x^3 +1.20942716*x^2 +  20.89095063*x+9.85633575}

#Precipitation
# Select the best function and implement it
source("Function/NEP_P.R")
stat_P(NEP)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = NEP,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_NEP<- function (x) {2.612202e+02*(1-exp(-1.931099e-03*x))}

# Append transform climate variables and stand age
NEP$f_P<- f_P_NEP(NEP$Annual_Preci)
NEP$f_Tair<- f_Tair_NEP(NEP$Tair)
NEP$f_Age<- f_Age_NEP(NEP$Stand_Age)

# Stepwise regression
lm.NEP<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=NEP)
step.NEP<- step(lm.NEP, direction = "backward")
print(step.NEP)
summary(step.NEP)

# Compute Importance variable
bootswiss <- boot.relimp(values~f_P + f_Tair + f_Age + f_P:f_Age + f_Tair:f_Age, 
                         data= NEP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
NEP<- NEP[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
NEP$prediction <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.NEP<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.NEP<- step(lm.NEP, direction = "backward")
  NEP.pred = predict(object = step.NEP, newdata = test.df)
  NEP$prediction[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$prediction, NEP$values)^2
RMSE_NEP <- rmse(NEP$prediction, NEP$values)
NSE_NEP<-NSE(NEP$prediction, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$prediction, NEP$values)

## Mean site
#-----------------------------------------------------------------
#Compute transform function
  
# Age
Fun_NEP_Mean_Site<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = NEP_Mean_Site,
                 start = list(A= 2.379e+02, B= -3.295e-03, C=-7.879e+02, D=-1.519e-01), control = list(maxiter = 500))
coef(Fun_NEP_Mean_Site)
f_Age_NEP_Mean_Site<- function (x) {1.999622e+02*(exp(3.885497e-04*x)) -9.927851e+02*(exp(-2.548488e-01*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_MAT.R")
stat_Temp(NEP_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = NEP_Mean_Site, 
                start = list(A=-0.02505, B=-0.11885, C= 28.03883, D=4.61699), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_NEP_Mean_Site<- function (x) {-0.1060114*x^3+2.8695169*x^2+5.9529467*x+3.8572301}

#Precipitation
# Select the best function and implement it
source("Function/NEP_P.R")
stat_P(NEP_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = NEP_Mean_Site,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_NEP_Mean_Site<- function (x) {2.608070e+02*(1-exp(-1.516355e-03*x))}

# Append transform climate variables and stand age
NEP_Mean_Site$f_P<- f_P_NEP_Mean_Site(NEP_Mean_Site$Annual_Preci)
NEP_Mean_Site$f_Tair<- f_Tair_NEP_Mean_Site(NEP_Mean_Site$Tair)
NEP_Mean_Site$f_Age<- f_Age_NEP_Mean_Site(NEP_Mean_Site$Stand_Age)

# Stepwise regression
lm.NEP_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=NEP_Mean_Site)
step.NEP_Mean_Site<- step(lm.NEP_Mean_Site, direction = "backward")
print(step.NEP_Mean_Site)
summary(step.NEP_Mean_Site)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Age, 
                         data= NEP_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
NEP_Mean_Site<- NEP_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
NEP_Mean_Site$prediction <- NA
for(id in unique(NEP_Mean_Site$Site_ID)){
  train.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID != id,]
  test.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.NEP_Mean_Site<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.NEP_Mean_Site<- step(lm.NEP_Mean_Site, direction = "backward")
  NEP_Mean_Site.pred = predict(object = step.NEP_Mean_Site, newdata = test.df)
  NEP_Mean_Site$prediction[NEP_Mean_Site$Site_ID == id] <- NEP_Mean_Site.pred
}

R2_NEP_Mean_Site<- cor(NEP_Mean_Site$prediction, NEP_Mean_Site$values)^2
RMSE_NEP_Mean_Site <- rmse(NEP_Mean_Site$prediction, NEP_Mean_Site$values)
NSE_NEP_Mean_Site<-NSE(NEP_Mean_Site$prediction, NEP_Mean_Site$values, na.rm=TRUE)
Bias_NEP_Mean_Site<-pbias(NEP_Mean_Site$prediction, NEP_Mean_Site$values)

# 1.3 Analysis for GPP

## All years per site
#-----------------------------------------------------------------#
# Compute transform function

# Age
Fun_GPP<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = GPP, 
               start = list(A=1287.1816, k= -0.1344), control = list(maxiter = 500))
coef(Fun_GPP)
f_Age_GPP<- function (x) {1350.7327181*(1-exp(-0.1322037*x))}

#Tair
# Select the best function and implement it
source("Function/GPP_MAT.R")
stat_GPP_Temp(GPP)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = GPP, 
                start = list(A=0.1914, B=74.6523, C= 691.5736), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_GPP<- function (x) {0.4075516*x^2+66.0585525*x+737.0781211}

#Precipitation
# Select the best function and implement it
source("Function/GPP_P.R")
stat_GPP_P(GPP)
Fun_Preci<-nlsLM(values~A*Annual_Preci^3+B*Annual_Preci^2+C*Annual_Preci+D, data = GPP, 
                 start = list(A=2.078e-07, B=-9.925e-04, C= 2.106e+00, D=1.038e+02), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_GPP<- function (x) {2.286300e-07*x^3-1.091184e-03*x^2+2.228790e+00*x+1.357104e+02}

# Append transform climate variables and stand age
GPP$f_P<- f_P_GPP(GPP$Annual_Preci)
GPP$f_Tair<- f_Tair_GPP(GPP$Tair)
GPP$f_Age<- f_Age_GPP(GPP$Stand_Age)

# Stepwise regression
lm.GPP<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=GPP)
step.GPP<- stepAIC(lm.GPP, direction = "backward")
print(step.GPP)
summary(step.GPP)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair, 
                         data= GPP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
GPP<- GPP[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
GPP$prediction <- NA
for(id in unique(GPP$Site_ID)){
  train.df <- GPP[GPP$Site_ID != id,]
  test.df <- GPP[GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.GPP<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.GPP<- step(lm.GPP, direction = "backward")
  GPP.pred = predict(object = step.GPP, newdata = test.df)
  GPP$prediction[GPP$Site_ID == id] <- GPP.pred
}

R2_GPP<- cor(GPP$prediction, GPP$values)^2
RMSE_GPP <- (sum((GPP$prediction-GPP$values)^2)/length(GPP$values))^(1/2)
NSE_GPP<-NSE(GPP$prediction, GPP$values, na.rm=TRUE)
Bias_GPP<-pbias(GPP$prediction, GPP$values)

## Mean site
-----------------------------------------------------------------
# Compute transform function
  
# Age
Fun_GPP_Mean_Site<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = GPP_Mean_Site, 
                 start = list(A=1287.1816, k= -0.1344), control = list(maxiter = 500))
coef(Fun_GPP_Mean_Site)
f_Age_GPP_Mean_Site<- function (x) {1318.1458225*(1-exp(-0.1346331*x))}

#Tair
# Select the best function and implement it
source("Function/GPP_MAT.R")
stat_GPP_Temp(GPP_Mean_Site)
Fun_Tair<-nlsLM(values~A/(1+exp(B-C*Tair)), data=GPP_Mean_Site,
                start = list(A= 2.239e+04, B=3.292e+00, C=5.549e-02), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_GPP_Mean_Site<- function (x) {1.445531e+06/(1+exp(7.624249e+00-5.974631e-02*x))}

#Precipitation
# Select the best function and implement it
source("Function/GPP_P.R")
stat_GPP_P(GPP_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = GPP_Mean_Site,
                 start = list(A= 3.168e+03, B=6.259e-04), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_GPP_Mean_Site<- function (x) {4.183496e+03*(1-exp(-4.453655e-04*x))}

# Append transform climate variables and stand age
GPP_Mean_Site$f_P<- f_P_GPP_Mean_Site(GPP_Mean_Site$Annual_Preci)
GPP_Mean_Site$f_Tair<- f_Tair_GPP_Mean_Site(GPP_Mean_Site$Tair)
GPP_Mean_Site$f_Age<- f_Age_GPP_Mean_Site(GPP_Mean_Site$Stand_Age)

# Stepwise regression
lm.GPP_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=GPP_Mean_Site)
step.GPP_Mean_Site<- stepAIC(lm.GPP_Mean_Site, direction = "backward")
print(step.GPP_Mean_Site)
summary(step.GPP_Mean_Site)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age, 
                         data= GPP_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
GPP_Mean_Site<- GPP_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
GPP_Mean_Site$prediction <- NA
for(id in unique(GPP_Mean_Site$Site_ID)){
  train.df <- GPP_Mean_Site[GPP_Mean_Site$Site_ID != id,]
  test.df <- GPP_Mean_Site[GPP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.GPP_Mean_Site<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.GPP_Mean_Site<- step(lm.GPP_Mean_Site, direction = "backward")
  GPP_Mean_Site.pred = predict(object = step.GPP_Mean_Site, newdata = test.df)
  GPP_Mean_Site$prediction[GPP_Mean_Site$Site_ID == id] <- GPP_Mean_Site.pred
}

R2_GPP_Mean_Site<- cor(GPP_Mean_Site$prediction, GPP_Mean_Site$values)^2
RMSE_GPP_Mean_Site <- (sum((GPP_Mean_Site$prediction-GPP_Mean_Site$values)^2)/length(GPP_Mean_Site$values))^(1/2)
NSE_GPP_Mean_Site<-NSE(GPP_Mean_Site$prediction, GPP_Mean_Site$values, na.rm=TRUE)
Bias_GPP_Mean_Site<-pbias(GPP_Mean_Site$prediction, GPP_Mean_Site$values)  

# 1.4 Analysis for Reco

## All years per site
#-----------------------------------------------------------------#

# Compute transform function
#Age
Fun_Reco<-nlsLM(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = Reco, 
                start = list(A = 631.614933, B = 0.154252, k = -0.001269), control = list(maxiter = 500))
coef(Fun_Reco)
f_Age_Reco<- function (x) {621.240665526*(x^0.183866864)*(exp(-0.002275885*x))}

#Tair
# Select the best function and implement it
source('Function/Reco_MAT.R')
stat_Reco_Temp(Reco)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Reco, 
                start = list(A= 0.4456, B=51.8351, C= 680.5344), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Reco<- function (x) {1.340493*x^2+31.803942*x+724.341953}

#Precipitation
# Select the best function and implement it
source("Function/Reco_P.R")
stat_Reco_P(Reco)
Fun_Preci<-nlsLM(values~A*Annual_Preci^2+B*Annual_Preci+C, data = Reco, 
                 start = list(A=3.583e-05, B=6.324e-01, C= 5.028e+02), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Reco<- function (x) {3.867042e-05*x^2+7.173393e-01*x+4.743947e+02}

# Append transform climate variables and stand age
Reco$f_P<- f_P_Reco(Reco$Annual_Preci)
Reco$f_Tair<- f_Tair_Reco(Reco$Tair)
Reco$f_Age<- f_Age_Reco(Reco$Stand_Age)

# Stepwise regression
lm.Reco<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Reco)
step.Reco<- stepAIC(lm.Reco, direction = "backward")
print(step.Reco)
summary(step.Reco)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair, 
                         data= Reco, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Reco<- Reco[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Reco$prediction <- NA
for(id in unique(Reco$Site_ID)){
  train.df <- Reco[Reco$Site_ID != id,]
  test.df <- Reco[Reco$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Reco<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Reco<- step(lm.Reco, direction = "backward")
  Reco.pred = predict(object = step.Reco, newdata = test.df)
  Reco$prediction[Reco$Site_ID == id] <- Reco.pred
}

R2_Reco<- cor(Reco$prediction, Reco$values)^2
RMSE_Reco <- (sum((Reco$prediction-Reco$values)^2)/length(Reco$values))^(1/2)
NSE_Reco<-NSE(Reco$prediction, Reco$values, na.rm=TRUE)
Bias_Reco<-pbias(Reco$prediction, Reco$values)  

## Mean site
-----------------------------------------------------------------
  
# Compute transform function
#Age
Fun_Reco_Mean_Site<-nlsLM(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = Reco_Mean_Site, 
                  start = list(A = 631.614933, B = 0.154252, k = -0.001269), control = list(maxiter = 500))
coef(Fun_Reco_Mean_Site)
f_Age_Reco_Mean_Site<- function (x) {595.699417750*(x^0.201937763)*(exp(-0.002877207*x))}

#Tair
# Select the best function and implement it
source('Function/Reco_MAT.R')
stat_Reco_Temp(Reco_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = Reco_Mean_Site, 
                start = list(A=0.4033, B=-12.0200, C= 125.5088, D=712.5760), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Reco_Mean_Site<- function (x) {0.3329494*x^3-8.1626636*x^2+80.7048041*x+729.3687097}

#Precipitation
# Select the best function and implement it
source("Function/Reco_P.R")
stat_Reco_P(Reco_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Reco_Mean_Site,
                 start = list(A= 2.723e+03, B=6.334e-04), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Reco_Mean_Site<- function (x) {4.683706e+03*(1-exp(-3.225884e-04*x))}

# Append transform climate variables and stand age
Reco_Mean_Site$f_P<- f_P_Reco_Mean_Site(Reco_Mean_Site$Annual_Preci)
Reco_Mean_Site$f_Tair<- f_Tair_Reco_Mean_Site(Reco_Mean_Site$Tair)
Reco_Mean_Site$f_Age<- f_Age_Reco_Mean_Site(Reco_Mean_Site$Stand_Age)

# Stepwise regression
lm.Reco_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Reco_Mean_Site)
step.Reco_Mean_Site<- stepAIC(lm.Reco_Mean_Site, direction = "backward")
print(step.Reco_Mean_Site)
summary(step.Reco_Mean_Site)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair, 
                         data= Reco_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Reco_Mean_Site<- Reco_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Reco_Mean_Site$prediction <- NA
for(id in unique(Reco_Mean_Site$Site_ID)){
  train.df <- Reco_Mean_Site[Reco_Mean_Site$Site_ID != id,]
  test.df <- Reco_Mean_Site[Reco_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Reco_Mean_Site<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Reco_Mean_Site<- step(lm.Reco_Mean_Site, direction = "backward")
  Reco_Mean_Site.pred = predict(object = step.Reco_Mean_Site, newdata = test.df)
  Reco_Mean_Site$prediction[Reco_Mean_Site$Site_ID == id] <- Reco_Mean_Site.pred
}

R2_Reco_Mean_Site<- cor(Reco_Mean_Site$prediction, Reco_Mean_Site$values)^2
RMSE_Reco_Mean_Site <- (sum((Reco_Mean_Site$prediction-Reco_Mean_Site$values)^2)/length(Reco_Mean_Site$values))^(1/2)
NSE_Reco_Mean_Site<-NSE(Reco_Mean_Site$prediction, Reco_Mean_Site$values, na.rm=TRUE)
Bias_Reco_Mean_Site<-pbias(Reco_Mean_Site$prediction, Reco_Mean_Site$values)  

# 1.5 Analysis for Ratio GPP-Reco

## All years per site
#-----------------------------------------------------------------
# Compute transform function

#Age
Fun_Ratio_GPP_Reco<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = Ratio_GPP_Reco, 
                         start = list(A= 1.1582, k= -0.2312), control = list(maxiter = 500))
coef(Fun_Ratio_GPP_Reco)
f_Age_Ratio_GPP_Reco<- function (x) {1.2013623*(1-exp(-0.2315584*x))}

#Tair
# Select the best function and implement it
source("Function/GPP_ER_MAT.R")
stat_GPP_ER_Temp(Ratio_GPP_Reco)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_GPP_Reco, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_GPP_Reco<- function (x) {-0.001131294*x^2+0.037115995*x+ 0.994279812}

#Precipitation
# Select the best function and implement it
source ("Function/GPP_ER_P.R")
stat_GPP_ER_P(Ratio_GPP_Reco)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data =Ratio_GPP_Reco,
                  start = list(A= 1.14235, B=0.00931), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_GPP_Reco<- function (x) {1.187800212*(1-exp(-0.009106726*x))}

# Append transform climate variables and stand age
Ratio_GPP_Reco$f_P<- f_P_Ratio_GPP_Reco(Ratio_GPP_Reco$Annual_Preci)
Ratio_GPP_Reco$f_Tair<- f_Tair_Ratio_GPP_Reco(Ratio_GPP_Reco$Tair)
Ratio_GPP_Reco$f_Age<- f_Age_Ratio_GPP_Reco(Ratio_GPP_Reco$Stand_Age)

# Stepwise regression
lm.Ratio_GPP_Reco<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_GPP_Reco)
step.Ratio_GPP_Reco<- stepAIC(lm.Ratio_GPP_Reco, direction = "backward")
print(step.Ratio_GPP_Reco)
summary(step.Ratio_GPP_Reco)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair, 
                         data= Ratio_GPP_Reco, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_GPP_Reco<- Ratio_GPP_Reco[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_GPP_Reco$prediction <- NA
for(id in unique(Ratio_GPP_Reco$Site_ID)){
  train.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID != id,]
  test.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_GPP_Reco<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Ratio_GPP_Reco<- step(lm.Ratio_GPP_Reco, direction = "backward")
  Ratio_GPP_Reco.pred = predict(object = step.Ratio_GPP_Reco, newdata = test.df)
  Ratio_GPP_Reco$prediction[Ratio_GPP_Reco$Site_ID == id] <- Ratio_GPP_Reco.pred
}

R2_Ratio_GPP_Reco<- cor(Ratio_GPP_Reco$prediction, Ratio_GPP_Reco$values)^2
RMSE_Ratio_GPP_Reco <- (sum((Ratio_GPP_Reco$prediction-Ratio_GPP_Reco$values)^2)/length(Ratio_GPP_Reco$values))^(1/2)
NSE_Ratio_GPP_Reco<-NSE(Ratio_GPP_Reco$prediction, Ratio_GPP_Reco$values, na.rm=TRUE)
Bias_Ratio_GPP_Reco<-pbias(Ratio_GPP_Reco$prediction, Ratio_GPP_Reco$values)  

## Mean site
-----------------------------------------------------------------
# Compute transform function
  
#Age
Fun_Ratio_GPP_Reco_Mean_Site<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = Ratio_GPP_Reco_Mean_Site, 
                            start = list(A= 1.1582, k= -0.2312), control = list(maxiter = 500))
coef(Fun_Ratio_GPP_Reco_Mean_Site)
f_Age_Ratio_GPP_Reco_Mean_Site<- function (x) {1.1980160*(1-exp(-0.2468009*x))}

#Tair
# Select the best function and implement it
source("Function/GPP_ER_MAT.R")
stat_GPP_ER_Temp(Ratio_GPP_Reco_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = Ratio_GPP_Reco_Mean_Site, 
                start = list(A=-1.553e-05, B=-7.161e-04, C= 3.382e-02, D=9.782e-01), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_GPP_Reco_Mean_Site<- function (x) {-0.0001560855*x^3+ 0.0042596229*x^2-0.0024131144*x +0.9909759452}

#Precipitation
# Select the best function and implement it
source ("Function/GPP_ER_P.R")
stat_GPP_ER_P(Ratio_GPP_Reco_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data =Ratio_GPP_Reco_Mean_Site,
                 start = list(A= 1.14235, B=0.00931), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_GPP_Reco_Mean_Site<- function (x) {1.19799534*(1-exp(-0.00605618*x))}

# Append transform climate variables and stand age
Ratio_GPP_Reco_Mean_Site$f_P<- f_P_Ratio_GPP_Reco_Mean_Site(Ratio_GPP_Reco_Mean_Site$Annual_Preci)
Ratio_GPP_Reco_Mean_Site$f_Tair<- f_Tair_Ratio_GPP_Reco_Mean_Site(Ratio_GPP_Reco_Mean_Site$Tair)
Ratio_GPP_Reco_Mean_Site$f_Age<- f_Age_Ratio_GPP_Reco_Mean_Site(Ratio_GPP_Reco_Mean_Site$Stand_Age)

# Stepwise regression
lm.Ratio_GPP_Reco_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_GPP_Reco_Mean_Site)
step.Ratio_GPP_Reco_Mean_Site<- stepAIC(lm.Ratio_GPP_Reco_Mean_Site, direction = "backward")
print(step.Ratio_GPP_Reco_Mean_Site)
summary(step.Ratio_GPP_Reco_Mean_Site)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair, 
                         data= Ratio_GPP_Reco_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_GPP_Reco_Mean_Site<- Ratio_GPP_Reco_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_GPP_Reco_Mean_Site$prediction <- NA
for(id in unique(Ratio_GPP_Reco_Mean_Site$Site_ID)){
  train.df <- Ratio_GPP_Reco_Mean_Site[Ratio_GPP_Reco_Mean_Site$Site_ID != id,]
  test.df <- Ratio_GPP_Reco_Mean_Site[Ratio_GPP_Reco_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_GPP_Reco_Mean_Site<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Ratio_GPP_Reco_Mean_Site<- step(lm.Ratio_GPP_Reco_Mean_Site, direction = "backward")
  Ratio_GPP_Reco_Mean_Site.pred = predict(object = step.Ratio_GPP_Reco_Mean_Site, newdata = test.df)
  Ratio_GPP_Reco_Mean_Site$prediction[Ratio_GPP_Reco_Mean_Site$Site_ID == id] <- Ratio_GPP_Reco_Mean_Site.pred
}

R2_Ratio_GPP_Reco_Mean_Site<- cor(Ratio_GPP_Reco_Mean_Site$prediction, Ratio_GPP_Reco_Mean_Site$values)^2
RMSE_Ratio_GPP_Reco_Mean_Site <- (sum((Ratio_GPP_Reco_Mean_Site$prediction-Ratio_GPP_Reco_Mean_Site$values)^2)/length(Ratio_GPP_Reco_Mean_Site$values))^(1/2)
NSE_Ratio_GPP_Reco_Mean_Site<-NSE(Ratio_GPP_Reco_Mean_Site$prediction, Ratio_GPP_Reco_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_GPP_Reco_Mean_Site<-pbias(Ratio_GPP_Reco_Mean_Site$prediction, Ratio_GPP_Reco_Mean_Site$values)  

# 1.6 Analysis for Ratio NEP-GPP

## All years per site
#-----------------------------------------------------------------
  
# Compute transform function

#Age
Fun_Ratio_NEP_GPP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPP,
                         start = list(A=0.165450, B= -0.003772, C=-1.319022, D=-0.148503), control = list(maxiter = 500))
coef(Fun_Ratio_NEP_GPP)
f_Age_Ratio_NEP_GPP<- function (x) {0.208727985*(exp(-0.004670439*x)) -1.530989643*(exp(-0.176894536*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_GPP_MAT.R")
stat_NEP_GPP_Temp(Ratio_NEP_GPP)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPP, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPP<- function (x) {-0.001054757*x^2+0.034763686*x-0.069793114}

#Precipitation
# Select the best function and implement it
source("Function/NEP_GPP_P.R")
stat_NEP_GPP_P(Ratio_NEP_GPP)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Ratio_NEP_GPP,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPP<- function (x) {0.10654012*(1-exp(-0.00215687*x))}

# Append transform climate variables and stand age
Ratio_NEP_GPP$f_P<- f_P_Ratio_NEP_GPP(Ratio_NEP_GPP$Annual_Preci)
Ratio_NEP_GPP$f_Tair<- f_Tair_Ratio_NEP_GPP(Ratio_NEP_GPP$Tair)
Ratio_NEP_GPP$f_Age<- f_Age_Ratio_NEP_GPP(Ratio_NEP_GPP$Stand_Age)

# Stepwise regression
lm.Ratio_NEP_GPP<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_NEP_GPP)
step.Ratio_NEP_GPP<- stepAIC(lm.Ratio_NEP_GPP, direction = "backward")
print(step.Ratio_NEP_GPP)
summary(step.Ratio_NEP_GPP)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair, 
                         data= Ratio_NEP_GPP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPP<- Ratio_NEP_GPP[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPP$prediction <- NA
for(id in unique(Ratio_NEP_GPP$Site_ID)){
  train.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != id,]
  test.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_NEP_GPP<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Ratio_NEP_GPP<- step(lm.Ratio_NEP_GPP, direction = "backward")
  Ratio_NEP_GPP.pred = predict(object = step.Ratio_NEP_GPP, newdata = test.df)
  Ratio_NEP_GPP$prediction[Ratio_NEP_GPP$Site_ID == id] <- Ratio_NEP_GPP.pred
}

R2_Ratio_NEP_GPP<- cor(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values)^2
RMSE_Ratio_NEP_GPP <- (sum((Ratio_NEP_GPP$prediction-Ratio_NEP_GPP$values)^2)/length(Ratio_NEP_GPP$values))^(1/2)
NSE_Ratio_NEP_GPP<-NSE(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values, na.rm=TRUE)
Bias_Ratio_NEP_GPP<-pbias(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values)  

## Mean site
-----------------------------------------------------------------
  
# Compute transform function
  
#Age
Fun_Ratio_NEP_GPP_Mean_Site<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPP_Mean_Site,
                           start = list(A=0.165450, B= -0.003772, C=-1.319022, D=-0.148503), control = list(maxiter = 500))
coef(Fun_Ratio_NEP_GPP_Mean_Site)
f_Age_Ratio_NEP_GPP_Mean_Site<- function (x) {0.1295690954*(exp(0.0008833587*x)) -1.4131439463*(exp(-0.1923803358*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_GPP_MAT.R")
stat_NEP_GPP_Temp(Ratio_NEP_GPP_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = Ratio_NEP_GPP_Mean_Site, 
                start = list(A=-1.553e-05, B=-7.161e-04, C= 3.382e-02, D=9.782e-01), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPP_Mean_Site<- function (x) {-0.0001115377*x^3+0.0029222585*x^2+0.0051844796*x-0.0932437869}

#Precipitation
# Select the best function and implement it
source("Function/NEP_GPP_P.R")
stat_NEP_GPP_P(Ratio_NEP_GPP_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Ratio_NEP_GPP_Mean_Site,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPP_Mean_Site<- function (x) {0.113883713*(1-exp(-0.001191419*x))}

# Append transform climate variables and stand age
Ratio_NEP_GPP_Mean_Site$f_P<- f_P_Ratio_NEP_GPP_Mean_Site(Ratio_NEP_GPP_Mean_Site$Annual_Preci)
Ratio_NEP_GPP_Mean_Site$f_Tair<- f_Tair_Ratio_NEP_GPP_Mean_Site(Ratio_NEP_GPP_Mean_Site$Tair)
Ratio_NEP_GPP_Mean_Site$f_Age<- f_Age_Ratio_NEP_GPP_Mean_Site(Ratio_NEP_GPP_Mean_Site$Stand_Age)

# Stepwise regression
lm.Ratio_NEP_GPP_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_NEP_GPP_Mean_Site)
step.Ratio_NEP_GPP_Mean_Site<- stepAIC(lm.Ratio_NEP_GPP_Mean_Site, direction = "backward")
print(step.Ratio_NEP_GPP_Mean_Site)
summary(step.Ratio_NEP_GPP_Mean_Site)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair + f_Tair:f_Age, 
                         data= Ratio_NEP_GPP_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPP_Mean_Site<- Ratio_NEP_GPP_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPP_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPP_Mean_Site$Site_ID)){
  train.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != id,]
  test.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_NEP_GPP_Mean_Site<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Ratio_NEP_GPP_Mean_Site<- step(lm.Ratio_NEP_GPP_Mean_Site, direction = "backward")
  Ratio_NEP_GPP_Mean_Site.pred = predict(object = step.Ratio_NEP_GPP_Mean_Site, newdata = test.df)
  Ratio_NEP_GPP_Mean_Site$prediction[Ratio_NEP_GPP_Mean_Site$Site_ID == id] <- Ratio_NEP_GPP_Mean_Site.pred
}

R2_Ratio_NEP_GPP_Mean_Site<- cor(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values)^2
RMSE_Ratio_NEP_GPP_Mean_Site <- (sum((Ratio_NEP_GPP_Mean_Site$prediction-Ratio_NEP_GPP_Mean_Site$values)^2)/length(Ratio_NEP_GPP_Mean_Site$values))^(1/2)
NSE_Ratio_NEP_GPP_Mean_Site<-NSE(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_NEP_GPP_Mean_Site<-pbias(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values)  

# 1.7 Analysis for Ratio NEP-GPPmax

## All years per site
#-----------------------------------------------------------------

# Compute transform function

#Age
Fun_Ratio_NEP_GPPmax<-Fun_Ratio_NEP_GPPmax<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPPmax,
                                                  start = list(A=-0.776705, B= -0.161076, C=0.189838, D=-0.002193), control = list(maxiter = 500))
coef(Fun_Ratio_NEP_GPPmax)
f_Age_Ratio_NEP_GPPmax<- function (x) {-0.848039211*(exp(-0.139083664*x)) + 0.276009334*(exp(-0.005064603*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_GPPmax_MAT.R")
stat_NEP_GPPmax_Temp(Ratio_NEP_GPPmax)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPPmax, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPPmax<- function (x) {-0.001049839*x^2+ 0.034868203*x-0.008556744}

#Precipitation
# Select the best function and implement it
source("Function/NEP_GPPmax_P.R")
stat_NEP_GPPmax_P(Ratio_NEP_GPPmax)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Ratio_NEP_GPPmax,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPPmax<- function (x) {0.169259704*(1-exp(-0.004557485*x))}

# Append transform climate variables and stand age
Ratio_NEP_GPPmax$f_P<- f_P_Ratio_NEP_GPPmax(Ratio_NEP_GPPmax$Annual_Preci)
Ratio_NEP_GPPmax$f_Tair<- f_Tair_Ratio_NEP_GPPmax(Ratio_NEP_GPPmax$Tair)
Ratio_NEP_GPPmax$f_Age<- f_Age_Ratio_NEP_GPPmax(Ratio_NEP_GPPmax$Stand_Age)

# Stepwise regression
lm.Ratio_NEP_GPPmax<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_NEP_GPPmax)
step.Ratio_NEP_GPPmax<- stepAIC(lm.Ratio_NEP_GPPmax, direction = "backward")
print(step.Ratio_NEP_GPPmax)
summary(step.Ratio_NEP_GPPmax)

# Compute Importance variable
bootswiss <- boot.relimp(values~ (f_P + f_Tair + f_Age)^2, 
                         data= Ratio_NEP_GPPmax, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPPmax<- Ratio_NEP_GPPmax[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPPmax$prediction <- NA
for(id in unique(Ratio_NEP_GPPmax$Site_ID)){
  train.df <- Ratio_NEP_GPPmax[Ratio_NEP_GPPmax$Site_ID != id,]
  test.df <- Ratio_NEP_GPPmax[Ratio_NEP_GPPmax$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_NEP_GPPmax<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Ratio_NEP_GPPmax<- step(lm.Ratio_NEP_GPPmax, direction = "backward")
  Ratio_NEP_GPPmax.pred = predict(object = step.Ratio_NEP_GPPmax, newdata = test.df)
  Ratio_NEP_GPPmax$prediction[Ratio_NEP_GPPmax$Site_ID == id] <- Ratio_NEP_GPPmax.pred
}

R2_Ratio_NEP_GPPmax<- cor(Ratio_NEP_GPPmax$prediction, Ratio_NEP_GPPmax$values)^2
RMSE_Ratio_NEP_GPPmax <- (sum((Ratio_NEP_GPPmax$prediction-Ratio_NEP_GPPmax$values)^2)/length(Ratio_NEP_GPPmax$values))^(1/2)
NSE_Ratio_NEP_GPPmax<-NSE(Ratio_NEP_GPPmax$prediction, Ratio_NEP_GPPmax$values, na.rm=TRUE)
Bias_Ratio_NEP_GPPmax<-pbias(Ratio_NEP_GPPmax$prediction, Ratio_NEP_GPPmax$values)

## Mean site
-----------------------------------------------------------------
  
# Compute transform function
  
#Age
Fun_Ratio_NEP_GPPmax_Mean_Site<-Fun_Ratio_NEP_GPPmax_Mean_Site<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPPmax_Mean_Site,
                                                    start = list(A=-0.776705, B= -0.161076, C=0.189838, D=-0.002193), control = list(maxiter = 500))
coef(Fun_Ratio_NEP_GPPmax_Mean_Site)
f_Age_Ratio_NEP_GPPmax_Mean_Site<- function (x) {-0.7628673769*(exp(-0.1905213298*x)) + 0.1659909077*(exp(0.0002685997*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_GPPmax_MAT.R")
stat_NEP_GPPmax_Temp(Ratio_NEP_GPPmax_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = Ratio_NEP_GPPmax, 
                start = list(A=-1.384e-05, B=-7.172e-04, C= 3.258e-02, D=-1.939e-02), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPPmax_Mean_Site<- function (x) {-6.116178e-05*x^3+ 8.995759e-04*x^2+2.224402e-02*x-1.055564e-02}

#Precipitation
# Select the best function and implement it
source("Function/NEP_GPPmax_P.R")
stat_NEP_GPPmax_P(Ratio_NEP_GPPmax_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Ratio_NEP_GPPmax_Mean_Site,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPPmax_Mean_Site<- function (x) {0.177883852*(1-exp(-0.002138833*x))}

# Append transform climate variables and stand age
Ratio_NEP_GPPmax_Mean_Site$f_P<- f_P_Ratio_NEP_GPPmax_Mean_Site(Ratio_NEP_GPPmax_Mean_Site$Annual_Preci)
Ratio_NEP_GPPmax_Mean_Site$f_Tair<- f_Tair_Ratio_NEP_GPPmax_Mean_Site(Ratio_NEP_GPPmax_Mean_Site$Tair)
Ratio_NEP_GPPmax_Mean_Site$f_Age<- f_Age_Ratio_NEP_GPPmax_Mean_Site(Ratio_NEP_GPPmax_Mean_Site$Stand_Age)

# Stepwise regression
lm.Ratio_NEP_GPPmax_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_NEP_GPPmax_Mean_Site)
step.Ratio_NEP_GPPmax_Mean_Site<- stepAIC(lm.Ratio_NEP_GPPmax_Mean_Site, direction = "backward")
print(step.Ratio_NEP_GPPmax_Mean_Site)
summary(step.Ratio_NEP_GPPmax_Mean_Site)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair, 
                         data= Ratio_NEP_GPPmax_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPPmax_Mean_Site<- Ratio_NEP_GPPmax_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPPmax_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPPmax_Mean_Site$Site_ID)){
  train.df <- Ratio_NEP_GPPmax_Mean_Site[Ratio_NEP_GPPmax_Mean_Site$Site_ID != id,]
  test.df <- Ratio_NEP_GPPmax_Mean_Site[Ratio_NEP_GPPmax_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_NEP_GPPmax_Mean_Site<- lm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Ratio_NEP_GPPmax_Mean_Site<- step(lm.Ratio_NEP_GPPmax_Mean_Site, direction = "backward")
  Ratio_NEP_GPPmax_Mean_Site.pred = predict(object = step.Ratio_NEP_GPPmax_Mean_Site, newdata = test.df)
  Ratio_NEP_GPPmax_Mean_Site$prediction[Ratio_NEP_GPPmax_Mean_Site$Site_ID == id] <- Ratio_NEP_GPPmax_Mean_Site.pred
}

R2_Ratio_NEP_GPPmax_Mean_Site<- cor(Ratio_NEP_GPPmax_Mean_Site$prediction, Ratio_NEP_GPPmax_Mean_Site$values)^2
RMSE_Ratio_NEP_GPPmax_Mean_Site <- (sum((Ratio_NEP_GPPmax_Mean_Site$prediction-Ratio_NEP_GPPmax_Mean_Site$values)^2)/length(Ratio_NEP_GPPmax_Mean_Site$values))^(1/2)
NSE_Ratio_NEP_GPPmax_Mean_Site<-NSE(Ratio_NEP_GPPmax_Mean_Site$prediction, Ratio_NEP_GPPmax_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_NEP_GPPmax_Mean_Site<-pbias(Ratio_NEP_GPPmax_Mean_Site$prediction, Ratio_NEP_GPPmax_Mean_Site$values)

# 2. Plot observed vs. actual for the different flux

# 2.1 All years per site

# Prepare dataset for plotting
pred_NEP<- NEP[c("Type_Flux", "values", "prediction")]
pred_GPP<- GPP[c("Type_Flux", "values", "prediction")]
pred_Reco<- Reco[c("Type_Flux", "values", "prediction")]
pred_Ratio_GPP_Reco<- Ratio_GPP_Reco[c("Type_Flux", "values", "prediction")]
pred_Ratio_NEP_GPP<- Ratio_NEP_GPP[c("Type_Flux", "values", "prediction")]
pred_Ratio_NEP_GPPmax<- Ratio_NEP_GPPmax[c("Type_Flux", "values", "prediction")]
pred_All<- rbind(pred_NEP, pred_GPP, pred_Reco, pred_Ratio_GPP_Reco, pred_Ratio_NEP_GPP, pred_Ratio_NEP_GPPmax)
levels(pred_All$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-ER", "Ratio NEP-GPP", 
                                "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "Ratio NEP-GPPclimax")

df_NEP<-pred_All[pred_All$Type_Flux %in% c("NEP"),]
df_GPP<-pred_All[pred_All$Type_Flux %in% c("GPP"),]
df_Reco<-pred_All[pred_All$Type_Flux %in% c("Respiration"),]
df_NEP_GPP<-pred_All[pred_All$Type_Flux %in% c("Ratio NEP-GPP"),]
df_GPP_Reco<-pred_All[pred_All$Type_Flux %in% c("Ratio GPP-ER"),]
df_NEP_GPPmax<-pred_All[pred_All$Type_Flux %in% c("Ratio NEP-GPPclimax"),]

#Plot prediction vs. observation

# NEP
gg1<- ggplot(df_NEP, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(-700, 800)+
  ylim(-700, 800)

#GPP
gg2<- ggplot(df_GPP, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(0, 4000)+
  ylim(0, 4000)

# Respiration
gg3<- ggplot(df_Reco, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(0, 4000)+
  ylim(0, 4000)

#Ratio GPP-Reco
gg4<- ggplot(df_GPP_Reco, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(0, 2)+
  ylim(0, 2)

# Ratio NEP-GPP
gg5<- ggplot(df_NEP_GPP, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(-1.5, 1)+
  ylim(-1.5, 1)

# Ratio NEP-GPPclimax
gg6<- ggplot(df_NEP_GPPmax, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(-1, 1.3)+
  ylim(-1, 1.3)

#Plot all plots together
pdf("Latex/Figures/Pred_Flux_All_Site.eps", width = 15, height = 10) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file

# 2.2 Average site

# Prepare dataset for plotting
pred_NEP_Mean_Site<- NEP_Mean_Site[c("Type_Flux", "values", "prediction")]
pred_GPP_Mean_Site<- GPP_Mean_Site[c("Type_Flux", "values", "prediction")]
pred_Reco_Mean_Site<- Reco_Mean_Site[c("Type_Flux", "values", "prediction")]
pred_Ratio_GPP_Reco_Mean_Site<- Ratio_GPP_Reco_Mean_Site[c("Type_Flux", "values", "prediction")]
pred_Ratio_NEP_GPP_Mean_Site<- Ratio_NEP_GPP_Mean_Site[c("Type_Flux", "values", "prediction")]
pred_Ratio_NEP_GPPmax_Mean_Site<- Ratio_NEP_GPPmax_Mean_Site[c("Type_Flux", "values", "prediction")]
pred_All<- rbind(pred_NEP_Mean_Site, pred_GPP_Mean_Site, pred_Reco_Mean_Site, pred_Ratio_GPP_Reco_Mean_Site, 
                 pred_Ratio_NEP_GPP_Mean_Site, pred_Ratio_NEP_GPPmax_Mean_Site)
levels(pred_All$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-ER", "Ratio NEP-GPP", 
                                "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4", "Ratio NEP-GPPclimax")

df_NEP_Mean_Site<-pred_All[pred_All$Type_Flux %in% c("NEP"),]
df_GPP_Mean_Site<-pred_All[pred_All$Type_Flux %in% c("GPP"),]
df_Reco_Mean_Site<-pred_All[pred_All$Type_Flux %in% c("Respiration"),]
df_NEP_GPP_Mean_Site<-pred_All[pred_All$Type_Flux %in% c("Ratio NEP-GPP"),]
df_GPP_Reco_Mean_site<-pred_All[pred_All$Type_Flux %in% c("Ratio GPP-ER"),]
df_NEP_GPPmax_Mean_Site<-pred_All[pred_All$Type_Flux %in% c("Ratio NEP-GPPclimax"),]

#Plot prediction vs. observation

# NEP
gg1<- ggplot(df_NEP_Mean_Site, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(-700, 800)+
  ylim(-700, 800)

#GPP
gg2<- ggplot(df_GPP_Mean_Site, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(0, 4000)+
  ylim(0, 4000)

# Respiration
gg3<- ggplot(df_Reco_Mean_Site, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(0, 4000)+
  ylim(0, 4000)

#Ratio GPP-Reco
gg4<- ggplot(df_GPP_Reco_Mean_site, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(0, 2)+
  ylim(0, 2)

# Ratio NEP-GPP
gg5<- ggplot(df_NEP_GPP_Mean_Site, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(-1.5, 1)+
  ylim(-1.5, 1)

# Ratio NEP-GPPclimax
gg6<- ggplot(df_NEP_GPPmax_Mean_Site, aes(x=prediction, y=values))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("Predicted")+
  ylab("Observed")+
  xlim(-1, 1.3)+
  ylim(-1, 1.3)

#Plot all plots together
pdf("Latex/Figures/Pred_Flux_Mean_Site.eps", width = 15, height = 10) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file