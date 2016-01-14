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

#1 Explain variabilty of the fluxes using random forest approach

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")

# 1.1. Subset dataframe for fitting process 

#Subset data set

Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]
NEP<- Flux_High[Flux_High$Type_Flux %in% c("NEP"),]
GPP<- Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco<- Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
Ratio_NEP_GPP<- Flux_High[Flux_High$Type_Flux %in% c("NEP_GPP"),]
Ratio_GPP_Reco<- Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]

#Compute GPPmax based on the lieth model
#Min model
params <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$Tair))))
GPP$GPPp<-with(as.list(params), GPP1000*((1-exp(-k*GPP$Annual_Preci))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))
#Multiplicative model
# params <- c(a=-0.119, b=1.315, c=-0.000664)
# GPP$GPPmat<-with(as.list(params), 1/(1+exp(a*GPP$Tair+b))) 
# GPP$GPPp<-with(as.list(params),(1-exp(c*GPP$Annual_Preci)))
# GPP<-transform(GPP, GPPmax = 3000*GPPmat*GPPp)
GPP$GPP_GPPmax<- NEP$values/GPP$GPPmax
Ratio_NEP_GPPmax<- GPP
Ratio_NEP_GPPmax<-Ratio_NEP_GPPmax[, !(colnames(Ratio_NEP_GPPmax) %in% c("values", "Type_Flux", "GPPmat", "GPPmax", "GPPp"))]
Ratio_NEP_GPPmax<-gather(Ratio_NEP_GPPmax, variable, values, -Annual_Preci, -year, 
                         -Ecosystem, -Climate, -Disturbance,
                         -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
                         -Tair, -Rg, -Study, -Lat, -Long)
Ratio_NEP_GPPmax<- Ratio_NEP_GPPmax[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "values", "Annual_Preci", 
                                      "Tair","Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")]
colnames(Ratio_NEP_GPPmax)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", 
                              "Tair", "Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")

# 1.2 Analysis for NEP

## All years per site
-----------------------------------------------------------------
#Compute transform function

# Age
Fun_NEP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = NEP,
               start = list(A= 2.379e+02, B= -3.295e-03, C=-7.879e+02, D=-1.519e-01), control = list(maxiter = 500))
coef(Fun_NEP)
f_Age_NEP<- function (x) {2.866217e+02*(exp(-4.939159e-03*x)) -1.013988e+03*(exp(-1.802970e-01*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_MAT.R")
stat_Temp(NEP)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = NEP, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_NEP<- function (x) {-0.9499793*x^2 + 32.6173090*x+18.9896462}

#Precipitation
# Select the best function and implement it
source("Function/NEP_P.R")
stat_P(NEP)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = NEP,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_NEP<- function (x) {2.560257e+02*(1-exp(-1.575503e-03*x))}

#Compute random forest with climate variables and stand age
NEPyoung<- subset(NEP, Stand_Age < 50)
NEPold<- subset(NEP, Stand_Age > 50)
NEP$f_P<- f_P_NEP(NEP$Annual_Preci)
NEP$f_Tair<- f_Tair_NEP(NEP$Tair)
NEP$f_Age<- f_Age_NEP(NEP$Stand_Age)
NEP<- NEP[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
NEP$prediction <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  NEP.ruf <- randomForest(values~Annual_Preci*Tair*Stand_Age, 
                                data = train.df,
                                importance = T, 
                                ntree=2000)
  NEP.ruf.pred = predict(object = NEP.ruf, newdata = test.df)
  NEP$prediction[NEP$Site_ID == id] <- NEP.ruf.pred
}

importance(NEP.ruf)
R2_NEP<- cor(NEP$prediction, NEP$values)^2
RMSE_NEP <- (sum((NEP$prediction-NEP$values)^2)/length(NEP$values))^(1/2)
NSE_NEP<-NSE(NEP$prediction, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$prediction, NEP$values)  

## Mean site
-----------------------------------------------------------------
#Compute transform function
  
# Age
Fun_NEP_Mean_Site<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = NEP_Mean_Site,
                 start = list(A= 2.379e+02, B= -3.295e-03, C=-7.879e+02, D=-1.519e-01), control = list(maxiter = 500))
coef(Fun_NEP_Mean_Site)
f_Age_NEP_Mean_Site<- function (x) {1.713517e+02*(exp(1.300205e-03*x)) -9.377685e+02*(exp(-2.526771e-01*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_MAT.R")
stat_Temp(NEP_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = NEP_Mean_Site, 
                start = list(A=-0.02505, B=-0.11885, C= 28.03883, D=4.61699), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_NEP_Mean_Site<- function (x) {-0.09597789*x^3+2.57706732*x^2+6.59768444*x+7.92810680}

#Precipitation
# Select the best function and implement it
source("Function/NEP_P.R")
stat_P(NEP_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = NEP_Mean_Site,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_NEP_Mean_Site<- function (x) {2.636420e+02*(1-exp(-1.222983e-03*x))}

#Compute random forest with climate variables and stand age
NEP_Mean_Siteyoung<- subset(NEP_Mean_Site, Stand_Age < 50)
NEP_Mean_Siteold<- subset(NEP_Mean_Site, Stand_Age > 50)
NEP_Mean_Site$f_P<- f_P_NEP_Mean_Site(NEP_Mean_Site$Annual_Preci)
NEP_Mean_Site$f_Tair<- f_Tair_NEP_Mean_Site(NEP_Mean_Site$Tair)
NEP_Mean_Site$f_Age<- f_Age_NEP_Mean_Site(NEP_Mean_Site$Stand_Age)
NEP_Mean_Site<- NEP_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
NEP_Mean_Site$prediction <- NA
for(id in unique(NEP_Mean_Site$Site_ID)){
  train.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID != id,]
  test.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  NEP_Mean_Site.ruf <- randomForest(values~Annual_Preci*Tair*Stand_Age, 
                          data = train.df,
                          importance = T, 
                          ntree=2000)
  NEP_Mean_Site.ruf.pred = predict(object = NEP_Mean_Site.ruf, newdata = test.df)
  NEP_Mean_Site$prediction[NEP_Mean_Site$Site_ID == id] <- NEP_Mean_Site.ruf.pred
}

importance(NEP_Mean_Site.ruf)
R2_NEP_Mean_Site<- cor(NEP_Mean_Site$prediction, NEP_Mean_Site$values)^2
RMSE_NEP_Mean_Site <- (sum((NEP_Mean_Site$prediction-NEP_Mean_Site$values)^2)/length(NEP_Mean_Site$values))^(1/2)
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
f_Age_GPP<- function (x) {1328.6538012*(1-exp(-0.1301034*x))}

#Tair
# Select the best function and implement it
source("Function/GPP_MAT.R")
stat_GPP_Temp(GPP)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = GPP, 
                start = list(A=0.1914, B=74.6523, C= 691.5736), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_GPP<- function (x) {0.526569*x^2+70.305293*x+719.849481}

#Precipitation
# Select the best function and implement it
source("Function/GPP_P.R")
stat_GPP_P(GPP)
Fun_Preci<-nlsLM(values~A*Annual_Preci^3+B*Annual_Preci^2+C*Annual_Preci+D, data = GPP, 
                 start = list(A=2.078e-07, B=-9.925e-04, C= 2.106e+00, D=1.038e+02), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_GPP<- function (x) {2.037852e-07*x^3-9.761853e-04*x^2+2.114104e+00*x+1.268553e+02}

# 1.3.2 Compute random forest with climate variables and stand age
GPPyoung<- subset(GPP, Stand_Age < 50)
GPPold<- subset(GPP, Stand_Age > 50)
GPP$f_P<- f_P_GPP(GPP$Annual_Preci)
GPP$f_Tair<- f_Tair_GPP(GPP$Tair)
GPP$f_Age<- f_Age_GPP(GPP$Stand_Age)
GPP<- GPP[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
GPP$prediction <- NA
for(id in unique(GPP$Site_ID)){
  train.df <- GPP[GPP$Site_ID != id,]
  test.df <- GPP[GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  GPP.ruf <- randomForest(values~ Annual_Preci + Tair + Stand_Age, 
                                 data = train.df,
                                 importance = T, 
                                 ntree=2000)
  GPP.ruf.pred = predict(object = GPP.ruf, newdata = test.df)
  GPP$prediction[GPP$Site_ID == id] <- GPP.ruf.pred
}

importance(GPP.ruf)
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
f_Age_GPP_Mean_Site<- function (x) {1296.9544039*(1-exp(-0.1275932*x))}

#Tair
# Select the best function and implement it
source("Function/GPP_MAT.R")
stat_GPP_Temp(GPP_Mean_Site)
Fun_Tair<-nlsLM(values~A/(1+exp(B-C*Tair)), data=GPP_Mean_Site,
                start = list(A= 2.239e+04, B=3.292e+00, C=5.549e-02), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_GPP_Mean_Site<- function (x) {1.174518e+06/(1+exp(7.417708e+00-6.025796e-02*x))}

#Precipitation
# Select the best function and implement it
source("Function/GPP_P.R")
stat_GPP_P(GPP_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = GPP_Mean_Site,
                 start = list(A= 3.168e+03, B=6.259e-04), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_GPP_Mean_Site<- function (x) {4.545596e+03*(1-exp(-3.875732e-04*x))}

# 1.3.2 Compute random forest with climate variables and stand age
GPP_Mean_Siteyoung<- subset(GPP_Mean_Site, Stand_Age < 50)
GPP_Mean_Siteold<- subset(GPP_Mean_Site, Stand_Age > 50)
GPP_Mean_Site$f_P<- f_P_GPP_Mean_Site(GPP_Mean_Site$Annual_Preci)
GPP_Mean_Site$f_Tair<- f_Tair_GPP_Mean_Site(GPP_Mean_Site$Tair)
GPP_Mean_Site$f_Age<- f_Age_GPP_Mean_Site(GPP_Mean_Site$Stand_Age)
GPP_Mean_Site<- GPP_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
GPP_Mean_Site$prediction <- NA
for(id in unique(GPP_Mean_Site$Site_ID)){
  train.df <- GPP_Mean_Site[GPP_Mean_Site$Site_ID != id,]
  test.df <- GPP_Mean_Site[GPP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  GPP_Mean_Site.ruf <- randomForest(values~ Annual_Preci + Tair + Stand_Age, 
                          data = train.df,
                          importance = T, 
                          ntree=2000)
  GPP_Mean_Site.ruf.pred = predict(object = GPP_Mean_Site.ruf, newdata = test.df)
  GPP_Mean_Site$prediction[GPP_Mean_Site$Site_ID == id] <- GPP_Mean_Site.ruf.pred
}

importance(GPP_Mean_Site.ruf)
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
f_Age_Reco<- function (x) {601.18948133*(x^0.18755369)*(exp(-0.00209698*x))}

#Tair
# Select the best function and implement it
source('Function/Reco_MAT.R')
stat_Reco_Temp(Reco)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Reco, 
                start = list(A= 0.4456, B=51.8351, C= 680.5344), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Reco<- function (x) {1.469048*x^2+37.822558*x+700.984723}

#Precipitation
# Select the best function and implement it
source("Function/Reco_P.R")
stat_Reco_P(Reco)
Fun_Preci<-nlsLM(values~A*Annual_Preci^2+B*Annual_Preci+C, data = Reco, 
                 start = list(A=3.583e-05, B=6.324e-01, C= 5.028e+02), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Reco<- function (x) {4.087312e-05*x^2+7.222884e-01*x+4.514801e+02}

# Compute random forest with climate variables and stand age
Recoyoung<- subset(Reco, Stand_Age < 50)
Recoold<- subset(Reco, Stand_Age > 50)
Reco$f_P<- f_P_Reco(Reco$Annual_Preci)
Reco$f_Tair<- f_Tair_Reco(Reco$Tair)
Reco$f_Age<- f_Age_Reco(Reco$Stand_Age)
Reco<- Reco[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Reco$prediction <- NA
for(id in unique(Reco$Site_ID)){
  train.df <- Reco[Reco$Site_ID != id,]
  test.df <- Reco[Reco$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  Reco.ruf <- randomForest(values~Annual_Preci + Tair + Stand_Age, 
                                 data = train.df,
                                 importance = T, 
                                 ntree=2000)
  Reco.ruf.pred = predict(object = Reco.ruf, newdata = test.df)
  Reco$prediction[Reco$Site_ID == id] <- Reco.ruf.pred
}

importance(Reco.ruf)
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
f_Age_Reco_Mean_Site<- function (x) {581.076641656*(x^0.203116629)*(exp(-0.002721904*x))}

#Tair
# Select the best function and implement it
source('Function/Reco_MAT.R')
stat_Reco_Temp(Reco_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = Reco_Mean_Site, 
                start = list(A=0.4033, B=-12.0200, C= 125.5088, D=712.5760), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Reco_Mean_Site<- function (x) {0.3202669*x^3-7.6982750*x^2+78.0542317*x+722.2485193}

#Precipitation
# Select the best function and implement it
source("Function/Reco_P.R")
stat_Reco_P(Reco_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Reco_Mean_Site,
                 start = list(A= 2.723e+03, B=6.334e-04), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Reco_Mean_Site<- function (x) {5.047789e+03*(1-exp(-2.887937e-04*x))}

# Compute random forest with climate variables and stand age
Reco_Mean_Siteyoung<- subset(Reco_Mean_Site, Stand_Age < 50)
Reco_Mean_Siteold<- subset(Reco_Mean_Site, Stand_Age > 50)
Reco_Mean_Site$f_P<- f_P_Reco_Mean_Site(Reco_Mean_Site$Annual_Preci)
Reco_Mean_Site$f_Tair<- f_Tair_Reco_Mean_Site(Reco_Mean_Site$Tair)
Reco_Mean_Site$f_Age<- f_Age_Reco_Mean_Site(Reco_Mean_Site$Stand_Age)
Reco_Mean_Site<- Reco_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Reco_Mean_Site$prediction <- NA
for(id in unique(Reco_Mean_Site$Site_ID)){
  train.df <- Reco_Mean_Site[Reco_Mean_Site$Site_ID != id,]
  test.df <- Reco_Mean_Site[Reco_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  Reco_Mean_Site.ruf <- randomForest(values~Annual_Preci + Tair + Stand_Age, 
                           data = train.df,
                           importance = T, 
                           ntree=2000)
  Reco_Mean_Site.ruf.pred = predict(object = Reco_Mean_Site.ruf, newdata = test.df)
  Reco_Mean_Site$prediction[Reco_Mean_Site$Site_ID == id] <- Reco_Mean_Site.ruf.pred
}

importance(Reco_Mean_Site.ruf)
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
f_Age_Ratio_GPP_Reco<- function (x) {1.1849425*(1-exp(-0.2357734*x))}

#Tair
# Select the best function and implement it
source("Function/GPP_ER_MAT.R")
stat_GPP_ER_Temp(Ratio_GPP_Reco)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_GPP_Reco, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_GPP_Reco<- function (x) {-0.00115323*x^2+ 0.03545035*x+ 1.00065137}

#Precipitation
# Select the best function and implement it
source ("Function/GPP_ER_P.R")
stat_GPP_ER_P(Ratio_GPP_Reco)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data =Ratio_GPP_Reco,
                  start = list(A= 1.14235, B=0.00931), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_GPP_Reco<- function (x) {1.168676201*(1-exp(-0.009425726*x))}

# Compute random forest with climate variables and stand age
Ratio_GPP_Recoyoung<- subset(Ratio_GPP_Reco, Stand_Age < 50)
Ratio_GPP_Recoold<- subset(Ratio_GPP_Reco, Stand_Age > 50)
Ratio_GPP_Reco$f_P<- f_P_Ratio_GPP_Reco(Ratio_GPP_Reco$Annual_Preci)
Ratio_GPP_Reco$f_Tair<- f_Tair_Ratio_GPP_Reco(Ratio_GPP_Reco$Tair)
Ratio_GPP_Reco$f_Age<- f_Age_Ratio_GPP_Reco(Ratio_GPP_Reco$Stand_Age)
Ratio_GPP_Reco<- Ratio_GPP_Reco[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_GPP_Reco$prediction <- NA
for(id in unique(Ratio_GPP_Reco$Site_ID)){
  train.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID != id,]
  test.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  Ratio_GPP_Reco.ruf <- randomForest(values~Annual_Preci + Tair + Stand_Age, 
                                 data = train.df,
                                 importance = T, 
                                 ntree=2000)
  Ratio_GPP_Reco.ruf.pred = predict(object = Ratio_GPP_Reco.ruf, newdata = test.df)
  Ratio_GPP_Reco$prediction[Ratio_GPP_Reco$Site_ID == id] <- Ratio_GPP_Reco.ruf.pred
}

importance(Ratio_GPP_Reco.ruf)
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
f_Age_Ratio_GPP_Reco_Mean_Site<- function (x) {1.1866492*(1-exp(-0.2433028*x))}

#Tair
# Select the best function and implement it
source("Function/GPP_ER_MAT.R")
stat_GPP_ER_Temp(Ratio_GPP_Reco_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = Ratio_GPP_Reco_Mean_Site, 
                start = list(A=-1.553e-05, B=-7.161e-04, C= 3.382e-02, D=9.782e-01), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_GPP_Reco_Mean_Site<- function (x) {-0.0001550743*x^3+ 0.0042187286*x^2-0.0021057427*x + 0.9916571099}

#Precipitation
# Select the best function and implement it
source ("Function/GPP_ER_P.R")
stat_GPP_ER_P(Ratio_GPP_Reco_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data =Ratio_GPP_Reco_Mean_Site,
                 start = list(A= 1.14235, B=0.00931), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_GPP_Reco_Mean_Site<- function (x) {1.184776664*(1-exp(-.006074198*x))}

# Compute random forest with climate variables and stand age
Ratio_GPP_Reco_Mean_Siteyoung<- subset(Ratio_GPP_Reco_Mean_Site, Stand_Age < 50)
Ratio_GPP_Reco_Mean_Siteold<- subset(Ratio_GPP_Reco_Mean_Site, Stand_Age > 50)
Ratio_GPP_Reco_Mean_Site$f_P<- f_P_Ratio_GPP_Reco_Mean_Site(Ratio_GPP_Reco_Mean_Site$Annual_Preci)
Ratio_GPP_Reco_Mean_Site$f_Tair<- f_Tair_Ratio_GPP_Reco_Mean_Site(Ratio_GPP_Reco_Mean_Site$Tair)
Ratio_GPP_Reco_Mean_Site$f_Age<- f_Age_Ratio_GPP_Reco_Mean_Site(Ratio_GPP_Reco_Mean_Site$Stand_Age)
Ratio_GPP_Reco_Mean_Site<- Ratio_GPP_Reco_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_GPP_Reco_Mean_Site$prediction <- NA
for(id in unique(Ratio_GPP_Reco_Mean_Site$Site_ID)){
  train.df <- Ratio_GPP_Reco_Mean_Site[Ratio_GPP_Reco_Mean_Site$Site_ID != id,]
  test.df <- Ratio_GPP_Reco_Mean_Site[Ratio_GPP_Reco_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  Ratio_GPP_Reco_Mean_Site.ruf <- randomForest(values~Annual_Preci + Tair + Stand_Age, 
                                     data = train.df,
                                     importance = T, 
                                     ntree=2000)
  Ratio_GPP_Reco_Mean_Site.ruf.pred = predict(object = Ratio_GPP_Reco_Mean_Site.ruf, newdata = test.df)
  Ratio_GPP_Reco_Mean_Site$prediction[Ratio_GPP_Reco_Mean_Site$Site_ID == id] <- Ratio_GPP_Reco_Mean_Site.ruf.pred
}

importance(Ratio_GPP_Reco_Mean_Site.ruf)
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
f_Age_Ratio_NEP_GPP<- function (x) {0.19945826*(exp(-0.00517618*x)) -1.51419788*(exp(-0.17645323*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_GPP_MAT.R")
stat_NEP_GPP_Temp(Ratio_NEP_GPP)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPP, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPP<- function (x) {-0.001058435*x^2+0.033521466*x-0.065887553}

#Precipitation
# Select the best function and implement it
source("Function/NEP_GPP_P.R")
stat_NEP_GPP_P(Ratio_NEP_GPP)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Ratio_NEP_GPP,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPP<- function (x) {0.10654012*(1-exp(-0.00215687*x))}

# Compute random forest with climate variables and stand age
Ratio_NEP_GPPyoung<- subset(Ratio_NEP_GPP, Stand_Age < 50)
Ratio_NEP_GPPold<- subset(Ratio_NEP_GPP, Stand_Age > 50)
Ratio_NEP_GPP$f_P<- f_P_Ratio_NEP_GPP(Ratio_NEP_GPP$Annual_Preci)
Ratio_NEP_GPP$f_Tair<- f_Tair_Ratio_NEP_GPP(Ratio_NEP_GPP$Tair)
Ratio_NEP_GPP$f_Age<- f_Age_Ratio_NEP_GPP(Ratio_NEP_GPP$Stand_Age)
Ratio_NEP_GPP<- Ratio_NEP_GPP[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPP$prediction <- NA
for(id in unique(Ratio_NEP_GPP$Site_ID)){
  train.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != id,]
  test.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  Ratio_NEP_GPP.ruf <- randomForest(values~Annual_Preci + Tair + Stand_Age, 
                                 data = train.df,
                                 importance = T, 
                                 ntree=2000)
  Ratio_NEP_GPP.ruf.pred = predict(object = Ratio_NEP_GPP.ruf, newdata = test.df)
  Ratio_NEP_GPP$prediction[Ratio_NEP_GPP$Site_ID == id] <- Ratio_NEP_GPP.ruf.pred
}

importance(Ratio_NEP_GPP.ruf)
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
f_Age_Ratio_NEP_GPP_Mean_Site<- function (x) {0.114421517*(exp(0.001525113*x)) -1.365016664*(exp(-0.184229780*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_GPP_MAT.R")
stat_NEP_GPP_Temp(Ratio_NEP_GPP_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = Ratio_NEP_GPP_Mean_Site, 
                start = list(A=-1.553e-05, B=-7.161e-04, C= 3.382e-02, D=9.782e-01), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPP_Mean_Site<- function (x) {-0.0001097155*x^3+0.0028769820*x^2+0.0051210453*x-0.0928790839}

#Precipitation
# Select the best function and implement it
source("Function/NEP_GPP_P.R")
stat_NEP_GPP_P(Ratio_NEP_GPP_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Ratio_NEP_GPP_Mean_Site,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPP_Mean_Site<- function (x) {0.1273282875*(1-exp(-0.0006899331*x))}

# Compute random forest with climate variables and stand age
Ratio_NEP_GPP_Mean_Siteyoung<- subset(Ratio_NEP_GPP_Mean_Site, Stand_Age < 50)
Ratio_NEP_GPP_Mean_Siteold<- subset(Ratio_NEP_GPP_Mean_Site, Stand_Age > 50)
Ratio_NEP_GPP_Mean_Site$f_P<- f_P_Ratio_NEP_GPP_Mean_Site(Ratio_NEP_GPP_Mean_Site$Annual_Preci)
Ratio_NEP_GPP_Mean_Site$f_Tair<- f_Tair_Ratio_NEP_GPP_Mean_Site(Ratio_NEP_GPP_Mean_Site$Tair)
Ratio_NEP_GPP_Mean_Site$f_Age<- f_Age_Ratio_NEP_GPP_Mean_Site(Ratio_NEP_GPP_Mean_Site$Stand_Age)
Ratio_NEP_GPP_Mean_Site<- Ratio_NEP_GPP_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPP_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPP_Mean_Site$Site_ID)){
  train.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != id,]
  test.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  Ratio_NEP_GPP_Mean_Site.ruf <- randomForest(values~Annual_Preci + Tair + Stand_Age, 
                                    data = train.df,
                                    importance = T, 
                                    ntree=2000)
  Ratio_NEP_GPP_Mean_Site.ruf.pred = predict(object = Ratio_NEP_GPP_Mean_Site.ruf, newdata = test.df)
  Ratio_NEP_GPP_Mean_Site$prediction[Ratio_NEP_GPP_Mean_Site$Site_ID == id] <- Ratio_NEP_GPP_Mean_Site.ruf.pred
}

importance(Ratio_NEP_GPP_Mean_Site.ruf)
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
f_Age_Ratio_NEP_GPPmax<- function (x) {-0.837843308*(exp(-0.137152089*x)) + 0.271401155*(exp(-0.006167279*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_GPPmax_MAT.R")
stat_NEP_GPPmax_Temp(Ratio_NEP_GPPmax)
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPPmax, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPPmax<- function (x) {-0.001077101*x^2+ 0.033039865*x-0.001576806}

#Precipitation
# Select the best function and implement it
source("Function/NEP_GPPmax_P.R")
stat_NEP_GPPmax_P(Ratio_NEP_GPPmax)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Ratio_NEP_GPPmax,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPPmax<- function (x) {-5.248299e+02*(1-exp(2.287424e-07*x))}

# 1.7.2 Compute random forest with climate variables and stand age
Ratio_NEP_GPPmaxyoung<- subset(Ratio_NEP_GPPmax, Stand_Age < 50)
Ratio_NEP_GPPmaxold<- subset(Ratio_NEP_GPPmax, Stand_Age > 50)
Ratio_NEP_GPPmax$f_P<- f_P_Ratio_NEP_GPPmax(Ratio_NEP_GPPmax$Annual_Preci)
Ratio_NEP_GPPmax$f_Tair<- f_Tair_Ratio_NEP_GPPmax(Ratio_NEP_GPPmax$Tair)
Ratio_NEP_GPPmax$f_Age<- f_Age_Ratio_NEP_GPPmax(Ratio_NEP_GPPmax$Stand_Age)
Ratio_NEP_GPPmax<- Ratio_NEP_GPPmax[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPPmax$prediction <- NA
for(id in unique(Ratio_NEP_GPPmax$Site_ID)){
  train.df <- Ratio_NEP_GPPmax[Ratio_NEP_GPPmax$Site_ID != id,]
  test.df <- Ratio_NEP_GPPmax[Ratio_NEP_GPPmax$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  Ratio_NEP_GPPmax.ruf <- randomForest(values~Annual_Preci + Tair + Stand_Age, 
                                 data = train.df,
                                 importance = T, 
                                 ntree=2000)
  Ratio_NEP_GPPmax.ruf.pred = predict(object = Ratio_NEP_GPPmax.ruf, newdata = test.df)
  Ratio_NEP_GPPmax$prediction[Ratio_NEP_GPPmax$Site_ID == id] <- Ratio_NEP_GPPmax.ruf.pred
}

importance(Ratio_NEP_GPPmax.ruf)
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
f_Age_Ratio_NEP_GPPmax_Mean_Site<- function (x) {-0.7055358303*(exp(-0.1742680533*x)) + 0.1447636669*(exp(0.0009525487*x))}

#Tair
# Select the best function and implement it
source("Function/NEP_GPPmax_MAT.R")
stat_NEP_GPPmax_Temp(Ratio_NEP_GPPmax_Mean_Site)
Fun_Tair<-nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = Ratio_NEP_GPPmax, 
                start = list(A=-1.384e-05, B=-7.172e-04, C= 3.258e-02, D=-1.939e-02), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPPmax_Mean_Site<- function (x) {-3.803532e-05*x^3+ 1.351748e-04*x^2+2.570041e-02*x-4.324706e-03}

#Precipitation
# Select the best function and implement it
source("Function/NEP_GPPmax_P.R")
stat_NEP_GPPmax_P(Ratio_NEP_GPPmax_Mean_Site)
Fun_Preci<-nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = Ratio_NEP_GPPmax_Mean_Site,
                 start = list(A= 300, B=0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPPmax_Mean_Site<- function (x) {0.169380076*(1-exp(-0.001807818*x))}

# 1.7.2 Compute random forest with climate variables and stand age
Ratio_NEP_GPPmax_Mean_Siteyoung<- subset(Ratio_NEP_GPPmax_Mean_Site, Stand_Age < 50)
Ratio_NEP_GPPmax_Mean_Siteold<- subset(Ratio_NEP_GPPmax_Mean_Site, Stand_Age > 50)
Ratio_NEP_GPPmax_Mean_Site$f_P<- f_P_Ratio_NEP_GPPmax_Mean_Site(Ratio_NEP_GPPmax_Mean_Site$Annual_Preci)
Ratio_NEP_GPPmax_Mean_Site$f_Tair<- f_Tair_Ratio_NEP_GPPmax_Mean_Site(Ratio_NEP_GPPmax_Mean_Site$Tair)
Ratio_NEP_GPPmax_Mean_Site$f_Age<- f_Age_Ratio_NEP_GPPmax_Mean_Site(Ratio_NEP_GPPmax_Mean_Site$Stand_Age)
Ratio_NEP_GPPmax_Mean_Site<- Ratio_NEP_GPPmax_Mean_Site[c("Site_ID",  "Type_Flux", "Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPPmax_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPPmax_Mean_Site$Site_ID)){
  train.df <- Ratio_NEP_GPPmax_Mean_Site[Ratio_NEP_GPPmax_Mean_Site$Site_ID != id,]
  test.df <- Ratio_NEP_GPPmax_Mean_Site[Ratio_NEP_GPPmax_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  Ratio_NEP_GPPmax_Mean_Site.ruf <- randomForest(values~Annual_Preci + Tair + Stand_Age, 
                                       data = train.df,
                                       importance = T, 
                                       ntree=2000)
  Ratio_NEP_GPPmax_Mean_Site.ruf.pred = predict(object = Ratio_NEP_GPPmax_Mean_Site.ruf, newdata = test.df)
  Ratio_NEP_GPPmax_Mean_Site$prediction[Ratio_NEP_GPPmax_Mean_Site$Site_ID == id] <- Ratio_NEP_GPPmax_Mean_Site.ruf.pred
}

importance(Ratio_NEP_GPPmax_Mean_Site.ruf)
R2_Ratio_NEP_GPPmax_Mean_Site<- cor(Ratio_NEP_GPPmax_Mean_Site$prediction, Ratio_NEP_GPPmax_Mean_Site$values)^2
RMSE_Ratio_NEP_GPPmax_Mean_Site <- (sum((Ratio_NEP_GPPmax_Mean_Site$prediction-Ratio_NEP_GPPmax_Mean_Site$values)^2)/length(Ratio_NEP_GPPmax_Mean_Site$values))^(1/2)
NSE_Ratio_NEP_GPPmax_Mean_Site<-NSE(Ratio_NEP_GPPmax_Mean_Site$prediction, Ratio_NEP_GPPmax_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_NEP_GPPmax_Mean_Site<-pbias(Ratio_NEP_GPPmax_Mean_Site$prediction, Ratio_NEP_GPPmax_Mean_Site$values)

# 2. Explain variabilty of the fluxes using linear regression analysis

#. 2.1 NEP

## All years per site
#-----------------------------------------------------------------

# Compute Importance variable
# Stepwise regression
lm.NEP<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=NEP)
step.NEP<- step(lm.NEP, direction = "backward")
print(step.NEP)
summary(step.NEP)

# VarImp estimation
bootswiss <- boot.relimp(values~f_P + f_Tair + f_Age + f_Tair:f_Age, 
                         data= NEP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
NEP$prediction <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.NEP<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.NEP<- step(lm.NEP, direction = "backward")
  NEP.pred = predict(object = step.NEP, newdata = test.df)
  NEP$prediction[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$prediction, NEP$values)^2
RMSE_NEP <- rmse(NEP$prediction, NEP$values)
NSE_NEP<-NSE(NEP$prediction, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$prediction, NEP$values)

## Mean site
-----------------------------------------------------------------
  
# Compute Importance variable
# Stepwise regression
lm.NEP_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=NEP_Mean_Site)
step.NEP_Mean_Site<- step(lm.NEP_Mean_Site, direction = "backward")
print(step.NEP_Mean_Site)
summary(step.NEP_Mean_Site)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Age, 
                         data= NEP_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
NEP_Mean_Site$prediction <- NA
for(id in unique(NEP_Mean_Site$Site_ID)){
  train.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID != id,]
  test.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.NEP_Mean_Site<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.NEP_Mean_Site<- step(lm.NEP_Mean_Site, direction = "backward")
  NEP_Mean_Site.pred = predict(object = step.NEP_Mean_Site, newdata = test.df)
  NEP_Mean_Site$prediction[NEP_Mean_Site$Site_ID == id] <- NEP_Mean_Site.pred
}

R2_NEP_Mean_Site<- cor(NEP_Mean_Site$prediction, NEP_Mean_Site$values)^2
RMSE_NEP_Mean_Site <- rmse(NEP_Mean_Site$prediction, NEP_Mean_Site$values)
NSE_NEP_Mean_Site<-NSE(NEP_Mean_Site$prediction, NEP_Mean_Site$values, na.rm=TRUE)
Bias_NEP_Mean_Site<-pbias(NEP_Mean_Site$prediction, NEP_Mean_Site$values)

#. 2.2 GPP

## All years per site
-----------------------------------------------------------------
  
# Compute Importance variable

# Stepwise regression
lm.GPP<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=GPP)
step.GPP<- stepAIC(lm.GPP, direction = "backward")
print(step.GPP)
summary(step.GPP)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_Tair:f_Age, 
                         data= GPP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
GPP$prediction <- NA
for(id in unique(GPP$Site_ID)){
  train.df <- GPP[GPP$Site_ID != id,]
  test.df <- GPP[GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.GPP<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
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
  
# Compute Importance variable
  
# Stepwise regression
lm.GPP_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=GPP_Mean_Site)
step.GPP_Mean_Site<- stepAIC(lm.GPP_Mean_Site, direction = "backward")
print(step.GPP_Mean_Site)
summary(step.GPP_Mean_Site)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age, 
                         data= GPP_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
GPP_Mean_Site$prediction <- NA
for(id in unique(GPP_Mean_Site$Site_ID)){
  train.df <- GPP_Mean_Site[GPP_Mean_Site$Site_ID != id,]
  test.df <- GPP_Mean_Site[GPP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.GPP_Mean_Site<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.GPP_Mean_Site<- step(lm.GPP_Mean_Site, direction = "backward")
  GPP_Mean_Site.pred = predict(object = step.GPP_Mean_Site, newdata = test.df)
  GPP_Mean_Site$prediction[GPP_Mean_Site$Site_ID == id] <- GPP_Mean_Site.pred
}

R2_GPP_Mean_Site<- cor(GPP_Mean_Site$prediction, GPP_Mean_Site$values)^2
RMSE_GPP_Mean_Site <- (sum((GPP_Mean_Site$prediction-GPP_Mean_Site$values)^2)/length(GPP_Mean_Site$values))^(1/2)
NSE_GPP_Mean_Site<-NSE(GPP_Mean_Site$prediction, GPP_Mean_Site$values, na.rm=TRUE)
Bias_GPP_Mean_Site<-pbias(GPP_Mean_Site$prediction, GPP_Mean_Site$values)  

#. 2.3 Respiration

## All years per site
-----------------------------------------------------------------
  
# Compute Importance variable

# Stepwise regression
lm.Reco<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Reco)
step.Reco<- stepAIC(lm.Reco, direction = "backward")
print(step.Reco)
summary(step.Reco)

# VarImp estimation
bootswiss <- boot.relimp(values~ (f_P + f_Tair + f_Age)^2, 
                           data= Reco, 
                           b = 100,  
                           type = "lmg",
                           rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
Reco$prediction <- NA
for(id in unique(Reco$Site_ID)){
  train.df <- Reco[Reco$Site_ID != id,]
  test.df <- Reco[Reco$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Reco<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
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
  
# Compute Importance variable
  
# Stepwise regression
lm.Reco_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Reco_Mean_Site)
step.Reco_Mean_Site<- stepAIC(lm.Reco_Mean_Site, direction = "backward")
print(step.Reco_Mean_Site)
summary(step.Reco_Mean_Site)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair, 
                         data= Reco_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
Reco_Mean_Site$prediction <- NA
for(id in unique(Reco_Mean_Site$Site_ID)){
  train.df <- Reco_Mean_Site[Reco_Mean_Site$Site_ID != id,]
  test.df <- Reco_Mean_Site[Reco_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Reco_Mean_Site<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Reco_Mean_Site<- step(lm.Reco_Mean_Site, direction = "backward")
  Reco_Mean_Site.pred = predict(object = step.Reco_Mean_Site, newdata = test.df)
  Reco_Mean_Site$prediction[Reco_Mean_Site$Site_ID == id] <- Reco_Mean_Site.pred
}

R2_Reco_Mean_Site<- cor(Reco_Mean_Site$prediction, Reco_Mean_Site$values)^2
RMSE_Reco_Mean_Site <- (sum((Reco_Mean_Site$prediction-Reco_Mean_Site$values)^2)/length(Reco_Mean_Site$values))^(1/2)
NSE_Reco_Mean_Site<-NSE(Reco_Mean_Site$prediction, Reco_Mean_Site$values, na.rm=TRUE)
Bias_Reco_Mean_Site<-pbias(Reco_Mean_Site$prediction, Reco_Mean_Site$values)  

#. 2.4 Ratio GPP-ER

## All years per site
-----------------------------------------------------------------
  
# Compute Importance variable

# Stepwise regression
lm.Ratio_GPP_Reco<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_GPP_Reco)
step.Ratio_GPP_Reco<- stepAIC(lm.Ratio_GPP_Reco, direction = "backward")
print(step.Ratio_GPP_Reco)
summary(step.Ratio_GPP_Reco)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair, 
                           data= Ratio_GPP_Reco, 
                           b = 100,  
                           type = "lmg",
                           rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_GPP_Reco$prediction <- NA
for(id in unique(Ratio_GPP_Reco$Site_ID)){
  train.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID != id,]
  test.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_GPP_Reco<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
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
  
# Compute Importance variable
  
# Stepwise regression
lm.Ratio_GPP_Reco_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_GPP_Reco_Mean_Site)
step.Ratio_GPP_Reco_Mean_Site<- stepAIC(lm.Ratio_GPP_Reco_Mean_Site, direction = "backward")
print(step.Ratio_GPP_Reco_Mean_Site)
summary(step.Ratio_GPP_Reco_Mean_Site)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair, 
                         data= Ratio_GPP_Reco_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_GPP_Reco_Mean_Site$prediction <- NA
for(id in unique(Ratio_GPP_Reco_Mean_Site$Site_ID)){
  train.df <- Ratio_GPP_Reco_Mean_Site[Ratio_GPP_Reco_Mean_Site$Site_ID != id,]
  test.df <- Ratio_GPP_Reco_Mean_Site[Ratio_GPP_Reco_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_GPP_Reco_Mean_Site<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Ratio_GPP_Reco_Mean_Site<- step(lm.Ratio_GPP_Reco_Mean_Site, direction = "backward")
  Ratio_GPP_Reco_Mean_Site.pred = predict(object = step.Ratio_GPP_Reco_Mean_Site, newdata = test.df)
  Ratio_GPP_Reco_Mean_Site$prediction[Ratio_GPP_Reco_Mean_Site$Site_ID == id] <- Ratio_GPP_Reco_Mean_Site.pred
}

R2_Ratio_GPP_Reco_Mean_Site<- cor(Ratio_GPP_Reco_Mean_Site$prediction, Ratio_GPP_Reco_Mean_Site$values)^2
RMSE_Ratio_GPP_Reco_Mean_Site <- (sum((Ratio_GPP_Reco_Mean_Site$prediction-Ratio_GPP_Reco_Mean_Site$values)^2)/length(Ratio_GPP_Reco_Mean_Site$values))^(1/2)
NSE_Ratio_GPP_Reco_Mean_Site<-NSE(Ratio_GPP_Reco_Mean_Site$prediction, Ratio_GPP_Reco_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_GPP_Reco_Mean_Site<-pbias(Ratio_GPP_Reco_Mean_Site$prediction, Ratio_GPP_Reco_Mean_Site$values)  

#. 2.5 Ratio NEP-GPP

## All years per site
-----------------------------------------------------------------

# Compute Importance variable

# Stepwise regression
lm.Ratio_NEP_GPP<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_NEP_GPP)
step.Ratio_NEP_GPP<- stepAIC(lm.Ratio_NEP_GPP, direction = "backward")
print(step.Ratio_NEP_GPP)
summary(step.Ratio_NEP_GPP)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Age, 
                         data= Ratio_NEP_GPP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPP$prediction <- NA
for(id in unique(Ratio_NEP_GPP$Site_ID)){
  train.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != id,]
  test.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_NEP_GPP<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
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
  
# Compute Importance variable
  
# Stepwise regression
lm.Ratio_NEP_GPP_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_NEP_GPP_Mean_Site)
step.Ratio_NEP_GPP_Mean_Site<- stepAIC(lm.Ratio_NEP_GPP_Mean_Site, direction = "backward")
print(step.Ratio_NEP_GPP_Mean_Site)
summary(step.Ratio_NEP_GPP_Mean_Site)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair + f_Tair:f_Age, 
                         data= Ratio_NEP_GPP_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPP_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPP_Mean_Site$Site_ID)){
  train.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != id,]
  test.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_NEP_GPP_Mean_Site<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Ratio_NEP_GPP_Mean_Site<- step(lm.Ratio_NEP_GPP_Mean_Site, direction = "backward")
  Ratio_NEP_GPP_Mean_Site.pred = predict(object = step.Ratio_NEP_GPP_Mean_Site, newdata = test.df)
  Ratio_NEP_GPP_Mean_Site$prediction[Ratio_NEP_GPP_Mean_Site$Site_ID == id] <- Ratio_NEP_GPP_Mean_Site.pred
}

R2_Ratio_NEP_GPP_Mean_Site<- cor(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values)^2
RMSE_Ratio_NEP_GPP_Mean_Site <- (sum((Ratio_NEP_GPP_Mean_Site$prediction-Ratio_NEP_GPP_Mean_Site$values)^2)/length(Ratio_NEP_GPP_Mean_Site$values))^(1/2)
NSE_Ratio_NEP_GPP_Mean_Site<-NSE(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_NEP_GPP_Mean_Site<-pbias(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values)  

#. 2.6 Ratio NEP-GPPmax

## All years per site
-----------------------------------------------------------------
  
# Compute Importance variable

# Stepwise regression
lm.Ratio_NEP_GPPmax<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_NEP_GPPmax)
step.Ratio_NEP_GPPmax<- stepAIC(lm.Ratio_NEP_GPPmax, direction = "backward")
print(step.Ratio_NEP_GPPmax)
summary(step.Ratio_NEP_GPPmax)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair + f_Tair:f_Age, 
                         data= Ratio_NEP_GPPmax, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPPmax$prediction <- NA
for(id in unique(Ratio_NEP_GPPmax$Site_ID)){
  train.df <- Ratio_NEP_GPPmax[Ratio_NEP_GPPmax$Site_ID != id,]
  test.df <- Ratio_NEP_GPPmax[Ratio_NEP_GPPmax$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_NEP_GPPmax<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
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
  
# Compute Importance variable
  
# Stepwise regression
lm.Ratio_NEP_GPPmax_Mean_Site<-lm(values ~ (f_P + f_Tair + f_Age)^2, data=Ratio_NEP_GPPmax_Mean_Site)
step.Ratio_NEP_GPPmax_Mean_Site<- stepAIC(lm.Ratio_NEP_GPPmax_Mean_Site, direction = "backward")
print(step.Ratio_NEP_GPPmax_Mean_Site)
summary(step.Ratio_NEP_GPPmax_Mean_Site)

# VarImp estimation
bootswiss <- boot.relimp(values~ f_P + f_Tair + f_Age + f_P:f_Tair, 
                         data= Ratio_NEP_GPPmax_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPPmax_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPPmax_Mean_Site$Site_ID)){
  train.df <- Ratio_NEP_GPPmax_Mean_Site[Ratio_NEP_GPPmax_Mean_Site$Site_ID != id,]
  test.df <- Ratio_NEP_GPPmax_Mean_Site[Ratio_NEP_GPPmax_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age")]
  lm.Ratio_NEP_GPPmax_Mean_Site<- glm(values ~ (f_P + f_Tair + f_Age)^2, data=train.df)
  step.Ratio_NEP_GPPmax_Mean_Site<- step(lm.Ratio_NEP_GPPmax_Mean_Site, direction = "backward")
  Ratio_NEP_GPPmax_Mean_Site.pred = predict(object = step.Ratio_NEP_GPPmax_Mean_Site, newdata = test.df)
  Ratio_NEP_GPPmax_Mean_Site$prediction[Ratio_NEP_GPPmax_Mean_Site$Site_ID == id] <- Ratio_NEP_GPPmax_Mean_Site.pred
}

R2_Ratio_NEP_GPPmax_Mean_Site<- cor(Ratio_NEP_GPPmax_Mean_Site$prediction, Ratio_NEP_GPPmax_Mean_Site$values)^2
RMSE_Ratio_NEP_GPPmax_Mean_Site <- (sum((Ratio_NEP_GPPmax_Mean_Site$prediction-Ratio_NEP_GPPmax_Mean_Site$values)^2)/length(Ratio_NEP_GPPmax_Mean_Site$values))^(1/2)
NSE_Ratio_NEP_GPPmax_Mean_Site<-NSE(Ratio_NEP_GPPmax_Mean_Site$prediction, Ratio_NEP_GPPmax_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_NEP_GPPmax_Mean_Site<-pbias(Ratio_NEP_GPPmax_Mean_Site$prediction, Ratio_NEP_GPPmax_Mean_Site$values)

# 3. Plot observed vs. actual for the different flux
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
pdf("Latex/Figures/Pred_Flux.eps", width = 15, height = 10) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file

# 4. Compute residual of fitted model vs. climate variables

# 4.1 NEP

# Compute residual
NEP$Res<- residuals(Fun_NEP)

#Plot residuals
R2_Tair<-(cor(residuals(Fun_NEP), NEP$Tair, use="pairwise.complete.obs"))^2
R2_Preci<-(cor(residuals(Fun_NEP), NEP$Annual_Preci, use="pairwise.complete.obs"))^2

gg1 <- ggplot(data = NEP, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.11", x = 23, y = -325) +
  facet_wrap(~Type_Flux)

gg2 <- ggplot(data = NEP, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.04", x = 3200, y = -325) +
  facet_wrap(~Type_Flux)

# 4.2 GPP

# Compute residual
GPP$Res<- residuals(Fun_GPP)

#Plot residuals
R2_Tair<-(cor(residuals(Fun_GPP), GPP$Tair, use="pairwise.complete.obs"))^2
R2_Preci<-(cor(residuals(Fun_GPP), GPP$Annual_Preci, use="pairwise.complete.obs"))^2

gg3 <- ggplot(data = GPP, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.51", x = 23, y = -800) +
  facet_wrap(~Type_Flux)

gg4 <- ggplot(data = GPP, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.53", x = 3200, y = -800) +
  facet_wrap(~Type_Flux)

# 4.3 Reco

# Compute residual
Reco$Res<- residuals(Fun_Reco)

#Plot residuals
R2_Tair<-(cor(residuals(Fun_Reco), Reco$Tair, use="pairwise.complete.obs"))^2
R2_Preci<-(cor(residuals(Fun_Reco), Reco$Annual_Preci, use="pairwise.complete.obs"))^2

gg5 <- ggplot(data = Reco, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.44", x = 23, y = -700) +
  facet_wrap(~Type_Flux)

gg6 <- ggplot(data = Reco, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.54", x = 3200, y = -700) +
  facet_wrap(~Type_Flux)

# 4.4 Ratio GPP-Reco
levels(Ratio_GPP_Reco$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-Reco", "Ratio NEP-GPP", "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4" )

# Compute residual
Ratio_GPP_Reco$Res<- residuals(Fun_Ratio_GPP_Reco)

#Plot residuals
R2_Tair<-lm(residuals(Fun_Ratio_GPP_Reco)~ Ratio_GPP_Reco$Tair + I(Ratio_GPP_Reco$Tair^2))
summary(R2_Tair)
R2_Preci<-lm(residuals(Fun_Ratio_GPP_Reco)~ Ratio_GPP_Reco$Annual_Preci + I(Ratio_GPP_Reco$Annual_Preci^2))
summary(R2_Preci)

gg7 <- ggplot(data = Ratio_GPP_Reco, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.22", x = 23, y = -0.4) +
  facet_wrap(~Type_Flux)

gg8 <- ggplot(data = Ratio_GPP_Reco, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.06", x = 3200, y = -0.4) +
  facet_wrap(~Type_Flux)

# 4.5 Ratio NEP-GPP
levels(Ratio_NEP_GPP$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-Reco", "Ratio NEP-GPP", "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4")

# Compute residual
Ratio_NEP_GPP$Res<- residuals(Fun_Ratio_NEP_GPP)

#Plot residuals
R2_Tair<-lm(residuals(Fun_Ratio_NEP_GPP)~ Ratio_NEP_GPP$Tair + I(Ratio_NEP_GPP$Tair^2))
summary(R2_Tair)
R2_Preci<-lm(residuals(Fun_Ratio_NEP_GPP)~ Ratio_NEP_GPP$Annual_Preci + I(Ratio_NEP_GPP$Annual_Preci^2))
summary(R2_Preci)

gg9 <- ggplot(data = Ratio_NEP_GPP, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.22", x = 23, y = -0.55) +
  facet_wrap(~Type_Flux)

gg10 <- ggplot(data = Ratio_NEP_GPP, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.06", x = 3200, y = -0.55) +
  facet_wrap(~Type_Flux)

# 4.6 Ratio NEP-GPPclimax

# Compute residual
Ratio_NEP_GPPmax$Res<- residuals(Fun_Ratio_NEP_GPPmax)

#Plot residuals
R2_Tair<-lm(residuals(Fun_Ratio_NEP_GPPmax)~ Ratio_NEP_GPPmax$Tair + I(Ratio_NEP_GPPmax$Tair^2))
summary(R2_Tair)
R2_Preci<-lm(residuals(Fun_Ratio_NEP_GPPmax)~ Ratio_NEP_GPPmax$Annual_Preci + I(Ratio_NEP_GPPmax$Annual_Preci^2))
summary(R2_Preci)
levels(Ratio_NEP_GPPmax$Type_Flux) <- c("Ratio NEP-GPPmax")

gg11 <- ggplot(data = Ratio_NEP_GPPmax, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.17", x = 23, y = -0.5)+ 
  facet_wrap(~Type_Flux)

gg12 <- ggplot(data = Ratio_NEP_GPPmax, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.03", x = 3200, y = -0.5) +
  facet_wrap(~Type_Flux)

# 4.7 Plot all graphs together
pdf("Latex/Figures/Residual_Flux.eps", width = 17, height = 14) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, gg10, gg11, gg12, nrow=6) # Write the grid.arrange in the file
dev.off() # Close the file
