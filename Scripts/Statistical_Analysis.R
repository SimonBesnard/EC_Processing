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
library(relaimpo)
library (caret)
library(party)
library (minpack.lm)
library (randomUniformForest)
library (reshape)

#1 Explain variabilty of the fluxes

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

# Remove climate variable outlier
NEP[c(37, 45, 136),"Annual_Preci"]<-NA
GPP[c(37,45,136),"Annual_Preci"]<-NA
Reco[c(37, 45, 136),"Annual_Preci"]<-NA
Ratio_GPP_Reco[c(37, 45, 136),"Annual_Preci"]<-NA
Ratio_NEP_GPP[c(37, 45, 136),"Annual_Preci"]<-NA

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
GPP$GPPmax[GPP$GPPmax == 0] <- NA
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

# Create dataset mean per site
NEP_Mean_Site<-ddply(NEP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Annual_Preci=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))

GPP_Mean_Site<-ddply(GPP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Annual_Preci=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))

Reco_Mean_Site<-ddply(Reco, .(Site_ID, Type_Flux),
                      summarise,
                      Stand_Age= mean(Stand_Age, na.rm=T),
                      values=mean(values, na.rm=T),
                      Annual_Preci=mean(Annual_Preci, na.rm=T),
                      Tair=mean(Tair, na.rm=T))

Ratio_NEP_GPP_Mean_Site<-ddply(Ratio_NEP_GPP, .(Site_ID, Type_Flux),
                               summarise,
                               Stand_Age= mean(Stand_Age, na.rm=T),
                               values=mean(values, na.rm=T),
                               Annual_Preci=mean(Annual_Preci, na.rm=T),
                               Tair=mean(Tair, na.rm=T))

Ratio_GPP_Reco_Mean_Site<-ddply(Ratio_GPP_Reco, .(Site_ID, Type_Flux),
                                summarise,
                                Stand_Age= mean(Stand_Age, na.rm=T),
                                values=mean(values, na.rm=T),
                                Annual_Preci=mean(Annual_Preci, na.rm=T),
                                Tair=mean(Tair, na.rm=T))

Ratio_NEP_GPPmax_Mean_Site<-ddply(Ratio_NEP_GPPmax, .(Site_ID, Type_Flux),
                                  summarise,
                                  Stand_Age= median(Stand_Age, na.rm=T),
                                  values=mean(values, na.rm=T),
                                  Annual_Preci=mean(Annual_Preci, na.rm=T),
                                  Tair=mean(Tair, na.rm=T))

# 1.2 Analysis for NEP

#1.2.1 Compute transform function

# Age
Fun_NEP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = NEP,
               start = list(A= 2.379e+02, B= -3.295e-03, C=-7.879e+02, D=-1.519e-01), control = list(maxiter = 500))
coef(Fun_NEP)
f_Age_NEP<- function (x) {2.378735e+02*(exp(-3.295241e-03*x)) -7.878536e+02*(exp(-1.518790e-01*x))}

#Tair
Fun_Tair<-nlsLM(values~A/(1+exp(1.1315-0.119*Tair)), data = NEP,
                start = list(A= 200), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_NEP<- function (x) {409.3021 /(1+exp(1.1315-0.119*x))}

#Precipitation
NEP$Annual_Preci[NEP$Annual_Preci == 0] <- NA
Fun_Preci<-nlsLM(values~A*(1-exp(-0.000664*Annual_Preci)), data = NEP,
                 start = list(A= 300), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_NEP<- function (x) {394.0764*(1-exp(-0.000664*x))}

# 1.2.2 Compute random forest with climate variables and stand age
NEPyoung<- subset(NEP, Stand_Age < 50)
NEPold<- subset(NEP, Stand_Age > 50)
NEP$f_P<- f_P_NEP(NEP$Annual_Preci)
NEP$f_Tair<- f_Tair_NEP(NEP$Tair)
NEP$f_Age<- f_Age_NEP(NEP$Stand_Age)
NEP<- NEP[c("Site_ID",  "Annual_Preci", "Tair","Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
NEP<- na.omit(NEP)
NEP$prediction <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "Rg", "f_P", "f_Tair", "f_Age")]
  NEP.ruf <- randomUniformForest(values~Annual_Preci+ Tair+ Stand_Age+ Rg + f_P + f_Tair + f_Age, 
                                data = train.df,
                                importance = T,
                                ntree = 2000)
  NEP.ruf.pred = predict(object = NEP.ruf, X = test.df)
  NEP$prediction[NEP$Site_ID == id] <- NEP.ruf.pred
}

summary(NEP.ruf)
statsModel = model.stats(NEP$prediction, NEP$values, regression = TRUE)
RMSE_NEP <- (sum((NEP$prediction-NEP$values)^2)/length(NEP$values))^(1/2)

# 1.3 Analysis for GPP

# 1.3.1 Compute transform function
# Age
Fun_GPP<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = GPP, 
               start = list(A=1287.1816, k= -0.1344), control = list(maxiter = 500))
coef(Fun_GPP)
f_Age_GPP<- function (x) {1287.1921959 *(1-exp(-0.1343273*x))}

#Tair
Fun_Tair<-nlsLM(values~A/(1+exp(1.1315-0.119*Tair)), data = GPP,
                start = list(A= 1500), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_GPP<- function (x) {2973.416/(1+exp(1.1315-0.119*x))}

#Precipitation
GPP$Annual_Preci[GPP$Annual_Preci == 0] <- NA
Fun_Preci<-nlsLM(values~A*(1-exp(-0.000664*Annual_Preci)), data = GPP,
                 start = list(A= 1500), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_GPP<- function (x) {3188.806 *(1-exp(-0.000664*x))}

# 1.3.2 Compute random forest with climate variables and stand age
GPPyoung<- subset(GPP, Stand_Age < 50)
GPPold<- subset(GPP, Stand_Age > 50)
GPP$f_P<- f_P_GPP(GPP$Annual_Preci)
GPP$f_Tair<- f_Tair_GPP(GPP$Tair)
GPP$f_Age<- f_Age_GPP(GPP$Stand_Age)
GPP<- GPP[ccontrols=cforest_unbiased()("Site_ID",  "Annual_Preci", "Tair","Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
GPP<- na.omit(GPP)
GPP$prediction <- NA
for(id in unique(GPP$Site_ID)){
  train.df <- GPP[GPP$Site_ID != id,]
  test.df <- GPP[GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "Rg", "f_P", "f_Tair", "f_Age")]
  GPP.ruf <- randomUniformForest(values~Annual_Preci+ Tair+ Stand_Age+ Rg + f_P + f_Tair + f_Age, 
                                 data = train.df,
                                 importance = T,
                                 ntree = 2000)
  GPP.ruf.pred = predict(object = GPP.ruf, X = test.df)
  GPP$prediction[GPP$Site_ID == id] <- GPP.ruf.pred
}

summary(GPP.ruf)
statsModel = model.stats(GPP$prediction, GPP$values, regression = TRUE)
RMSE_GPP <- (sum((GPP$prediction-GPP$values)^2)/length(GPP$values))^(1/2)

# 1.4 Analysis for Reco

# 1.4.1 Compute transform function
#Age
Fun_Reco<-nlsLM(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = Reco, 
                start = list(A = 631.614933, B = 0.154252, k = -0.001269), control = list(maxiter = 500))
coef(Fun_Reco)
f_Age_Reco<- function (x) {631.617082164*(x^0.154251100)*(exp(-0.001269377 *x))}

#Tair
Fun_Tair<-nlsLM(values~A/(1+exp(1.1315-0.119*Tair)), data = Reco,
                start = list(A= 1500), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Reco<- function (x) {2591.7/(1+exp(1.1315-0.119*x))}

#Precipitation
Reco$Annual_Preci[Reco$Annual_Preci == 0] <- NA
Fun_Preci<-nlsLM(values~A*(1-exp(-0.000664*Annual_Preci)), data = Reco,
                 start = list(A= 1500), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Reco<- function (x) {2799.827 *(1-exp(-0.000664*x))}

# 1.4.2 Compute random forest with climate variables and stand age
Recoyoung<- subset(Reco, Stand_Age < 50)
Recoold<- subset(Reco, Stand_Age > 50)
Reco$f_P<- f_P_Reco(Reco$Annual_Preci)
Reco$f_Tair<- f_Tair_Reco(Reco$Tair)
Reco$f_Age<- f_Age_Reco(Reco$Stand_Age)
Reco<- Reco[c("Site_ID",  "Annual_Preci", "Tair","Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Reco<- na.omit(Reco)
Reco$prediction <- NA
for(id in unique(Reco$Site_ID)){
  train.df <- Reco[Reco$Site_ID != id,]
  test.df <- Reco[Reco$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "Rg", "f_P", "f_Tair", "f_Age")]
  Reco.ruf <- randomUniformForest(values~Annual_Preci+ Tair+ Stand_Age+ Rg + f_P + f_Tair + f_Age, 
                                 data = train.df,
                                 importance = T,
                                 ntree = 2000, na.action="omit")
  Reco.ruf.pred = predict(object = Reco.ruf, X = test.df)
  Reco$prediction[Reco$Site_ID == id] <- Reco.ruf.pred
}

summary(Reco.ruf)
statsModel = model.stats(Reco$prediction, Reco$values, regression = TRUE)
RMSE_Reco <- (sum((Reco$prediction-Reco$values)^2)/length(Reco$values))^(1/2)

# 1.5 Analysis for Ratio GPP-Reco

#1.5.1 Compute transform function

#Age
Fun_Ratio_GPP_Reco<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = Ratio_GPP_Reco, 
                         start = list(A= 1.1582, k= -0.2312), control = list(maxiter = 500))
coef(Fun_Ratio_GPP_Reco)
f_Age_Ratio_GPP_Reco<- function (x) {1.1582127*(1-exp(-0.2312178*x))}

#Tair
Fun_Tair<-nlsLM(values~A/(1+exp(B-C*Tair)), data = Ratio_GPP_Reco,
                start = list(A= 1.2, B=1.1315, C=0.119), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_GPP_Reco<- function (x) {1.2626138/(1+exp(-1.2851635-0.19*x))}

#Precipitation
Ratio_GPP_Reco$Annual_Preci[Ratio_GPP_Reco$Annual_Preci == 0] <- NA
Fun_Preci<-nlsLM(values~A*(1-exp(B*Annual_Preci)), data = Ratio_GPP_Reco,
                 start = list(A= 1.2, B=-0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_GPP_Reco<- function (x) {1.148760390*(1-exp(-0.008826998*x))}

# 1.5.2 Compute random forest with climate variables and stand age
Ratio_GPP_Recoyoung<- subset(Ratio_GPP_Reco, Stand_Age < 50)
Ratio_GPP_Recoold<- subset(Ratio_GPP_Reco, Stand_Age > 50)
Ratio_GPP_Reco$f_P<- f_P_Ratio_GPP_Reco(Ratio_GPP_Reco$Annual_Preci)
Ratio_GPP_Reco$f_Tair<- f_Tair_Ratio_GPP_Reco(Ratio_GPP_Reco$Tair)
Ratio_GPP_Reco$f_Age<- f_Age_Ratio_GPP_Reco(Ratio_GPP_Reco$Stand_Age)
Ratio_GPP_Reco<- Ratio_GPP_Reco[c("Site_ID",  "Annual_Preci", "Tair","Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_GPP_Reco<- na.omit(Ratio_GPP_Reco)
Ratio_GPP_Reco$prediction <- NA
for(id in unique(Ratio_GPP_Reco$Site_ID)){
  train.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID != id,]
  test.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "Rg", "f_P", "f_Tair", "f_Age")]
  Ratio_GPP_Reco.ruf <- randomUniformForest(values~Annual_Preci+ Tair+ Stand_Age+ Rg + f_P + f_Tair + f_Age, 
                                 data = train.df,
                                 importance = T,
                                 ntree = 2000, na.action="omit")
  Ratio_GPP_Reco.ruf.pred = predict(object = Ratio_GPP_Reco.ruf, X = test.df)
  Ratio_GPP_Reco$prediction[Ratio_GPP_Reco$Site_ID == id] <- Ratio_GPP_Reco.ruf.pred
}

summary(Ratio_GPP_Reco.ruf)
statsModel = model.stats(Ratio_GPP_Reco$prediction, Ratio_GPP_Reco$values, regression = TRUE)
RMSE_Ratio_GPP_Reco <- (sum((Ratio_GPP_Reco$prediction-Ratio_GPP_Reco$values)^2)/length(Ratio_GPP_Reco$values))^(1/2)

# 1.6 Analysis for Ratio NEP-GPP

# 1.6.1 Compute transform function

#Age
Fun_Ratio_NEP_GPP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPP,
                         start = list(A=0.165450, B= -0.003772, C=-1.319022, D=-0.148503), control = list(maxiter = 500))
coef(Fun_Ratio_NEP_GPP)
f_Age_Ratio_NEP_GPP<- function (x) {0.165451439*(exp(-0.003771699*x)) -1.319022091*(exp(-0.148502307*x))}

#Tair
Fun_Tair<-nlsLM(values~A/(1+exp(1.1315-0.119*Tair)), data = Ratio_NEP_GPP,
                start = list(A= 0.2), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPP<- function (x) {0.2069852/(1+exp(1.1315-0.119*x))}

#Precipitation
Ratio_NEP_GPP$Annual_Preci[Ratio_NEP_GPP$Annual_Preci == 0] <- NA
Fun_Preci<-nlsLM(values~A*(1-exp(B*Annual_Preci)), data = Ratio_NEP_GPP,
                 start = list(A= 0.2, B=-0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPP<- function (x) { 0.095920571*(1-exp(-0.001707187*x))}

# 1.6.2 Compute random forest with climate variables and stand age
Ratio_NEP_GPPyoung<- subset(Ratio_NEP_GPP, Stand_Age < 50)
Ratio_NEP_GPPold<- subset(Ratio_NEP_GPP, Stand_Age > 50)
Ratio_NEP_GPP$f_P<- f_P_Ratio_NEP_GPP(Ratio_NEP_GPP$Annual_Preci)
Ratio_NEP_GPP$f_Tair<- f_Tair_Ratio_NEP_GPP(Ratio_NEP_GPP$Tair)
Ratio_NEP_GPP$f_Age<- f_Age_Ratio_NEP_GPP(Ratio_NEP_GPP$Stand_Age)
Ratio_NEP_GPP<- Ratio_NEP_GPP[c("Site_ID",  "Annual_Preci", "Tair","Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPP<- na.omit(Ratio_NEP_GPP)
Ratio_NEP_GPP$prediction <- NA
for(id in unique(Ratio_NEP_GPP$Site_ID)){
  train.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != id,]
  test.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "Rg", "f_P", "f_Tair", "f_Age")]
  Ratio_NEP_GPP.ruf <- randomUniformForest(values~Annual_Preci+ Tair+ Stand_Age+ Rg + f_P + f_Tair + f_Age, 
                                 data = train.df,
                                 importance = T,
                                 ntree = 2000, na.action="omit")
  Ratio_NEP_GPP.ruf.pred = predict(object = Ratio_NEP_GPP.ruf, X = test.df)
  Ratio_NEP_GPP$prediction[Ratio_NEP_GPP$Site_ID == id] <- Ratio_NEP_GPP.ruf.pred
}

summary(Ratio_NEP_GPP.ruf)
statsModel = model.stats(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values, regression = TRUE)
RMSE_Ratio_NEP_GPP <- (sum((Ratio_NEP_GPP$prediction-Ratio_NEP_GPP$values)^2)/length(Ratio_NEP_GPP$values))^(1/2)

# 1.7 Analysis for Ratio NEP-GPPmax

# 1.7.1 Compute transform function

#Age
Fun_Ratio_NEP_GPPmax<-nlsLM(values~A+B*Stand_Age^(C)*exp(D*Stand_Age)+E/(1+exp(-Stand_Age*H)), data =  Ratio_NEP_GPPmax, 
                           start = list(A = -2.198508, B = -0.000103, C = 13.611050, D=-2.748321, E=2.312422, H=0.276304), control = list(maxiter = 500))
coef(Fun_Ratio_NEP_GPPmax)
f_Age_Ratio_NEP_GPPmax<- function (x) {-2.143196e+00-2.588261e-06*x^(2.402659e+01)*exp(-5.176585e+00*x)+2.280457e+00/(1+exp(-x*9.254190e-01))}

#Tair
Fun_Tair<-nlsLM(values~A/(1+exp(B-C*Tair)), data = Ratio_NEP_GPPmax,
                start = list(A= 0.2, B=1.1315, C=0.119), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPPmax<- function (x) {0.1907508/(1+exp(1.1315-0.119*x))}

#Precipitation
Ratio_NEP_GPPmax$Annual_Preci[Ratio_NEP_GPPmax$Annual_Preci == 0] <- NA
Fun_Preci<-nlsLM(values~A*(1-exp(B*Annual_Preci)), data = Ratio_NEP_GPPmax,
                 start = list(A= 0.2, B=-0.000664), control = list(maxiter = 500))
coef(Fun_Preci)
f_P_Ratio_NEP_GPPmax<- function (x) {0.156992052*(1-exp(-0.002907185*x))}

# 1.7.2 Compute random forest with climate variables and stand age
Ratio_NEP_GPPmaxyoung<- subset(Ratio_NEP_GPPmax, Stand_Age < 50)
Ratio_NEP_GPPmaxold<- subset(Ratio_NEP_GPPmax, Stand_Age > 50)
Ratio_NEP_GPPmax$f_P<- f_P_Ratio_NEP_GPPmax(Ratio_NEP_GPPmax$Annual_Preci)
Ratio_NEP_GPPmax$f_Tair<- f_Tair_Ratio_NEP_GPPmax(Ratio_NEP_GPPmax$Tair)
Ratio_NEP_GPPmax$f_Age<- f_Age_Ratio_NEP_GPPmax(Ratio_NEP_GPPmax$Stand_Age)
Ratio_NEP_GPPmax<- Ratio_NEP_GPPmax[c("Site_ID",  "Annual_Preci", "Tair","Rg", "Stand_Age", "f_P", "f_Tair", "f_Age", "values")]
Ratio_NEP_GPPmax<- na.omit(Ratio_NEP_GPPmax)
Ratio_NEP_GPPmax$prediction <- NA
for(id in unique(Ratio_NEP_GPPmax$Site_ID)){
  train.df <- Ratio_NEP_GPPmax[Ratio_NEP_GPPmax$Site_ID != id,]
  test.df <- Ratio_NEP_GPPmax[Ratio_NEP_GPPmax$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "Rg", "f_P", "f_Tair", "f_Age")]
  Ratio_NEP_GPPmax.ruf <- randomUniformForest(values~Annual_Preci+ Tair+ Stand_Age+ Rg + f_P + f_Tair + f_Age, 
                                 data = train.df,
                                 importance = T,
                                 ntree = 2000, na.action="omit")
  Ratio_NEP_GPPmax.ruf.pred = predict(object = Ratio_NEP_GPPmax.ruf, X = test.df)
  Ratio_NEP_GPPmax$prediction[Ratio_NEP_GPPmax$Site_ID == id] <- Ratio_NEP_GPPmax.ruf.pred
}

summary(Ratio_NEP_GPPmax.ruf)
statsModel = model.stats(Ratio_NEP_GPPmax$prediction, Ratio_NEP_GPPmax$values, regression = TRUE)
RMSE_Ratio_NEP_GPPmax <- (sum((Ratio_NEP_GPPmax$prediction-Ratio_NEP_GPPmax$values)^2)/length(Ratio_NEP_GPPmax$values))^(1/2)

# 1.8. Plot relative contribution output

#Subset data
df<-read.csv("Output/VarImp_Flux.csv", header = TRUE)
df_NEP<-df[df$Flux %in% c("NEP"),]
df_GPP<-df[df$Flux %in% c("GPP"),]
df_Reco<-df[df$Flux %in% c("Respiration"),]
df_NEP_GPP<-df[df$Flux %in% c("Ratio NEP-GPP"),]
df_GPP_Reco<-df[df$Flux %in% c("Ratio GPP-ER"),]
df_NEP_GPPmax<-df[df$Flux %in% c("Ratio NEP-GPPclimax"),]

# Create plot per flux
gg1<- ggplot(df_NEP, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")+
  ylab("Relative contribution to NEP (%)")

gg2<- ggplot(df_GPP, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")+
  ylab("Relative contribution to GPP (%)")

gg3<- ggplot(df_Reco, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")+
  ylab("Relative contribution to Reco (%)")

gg4<- ggplot(df_GPP_Reco, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")+
  ylab("Relative contribution to GPP/Reco (%)")

gg5<- ggplot(df_NEP_GPP, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")+
  ylab("Relative contribution to NEP/GPP (%)")

gg6<- ggplot(df_NEP_GPPmax, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")+
  ylab("Relative contribution to NEP/GPPclimax (%)")

#Plot all plots together
pdf("Latex/Figures/VarImp_Flux.eps", width = 15, height = 12) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file

# 2. Compute residual of fitted model vs. climate variables

# 2.1 NEP

# Compute residual
NEP$Res<- residuals(Fun_NEP)
NEP$Annual_Preci[NEP$Annual_Preci == 0] <- NA

#Plot residuals
R2_Tair<-(cor(residuals(Fun_NEP), NEP$Tair, use="pairwise.complete.obs"))^2
R2_Preci<-(cor(residuals(Fun_NEP), NEP$Annual_Preci, use="pairwise.complete.obs"))^2

gg1 <- ggplot(data = NEP, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.19", x = 23, y = -400) +
  facet_wrap(~Type_Flux)

gg2 <- ggplot(data = NEP, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.07", x = 3200, y = -400) +
  facet_wrap(~Type_Flux)

# 2.2 GPP

# Compute residual
GPP$Res<- residuals(Fun_GPP)
GPP$Annual_Preci[GPP$Annual_Preci == 0] <- NA

#Plot residuals
R2_Tair<-(cor(residuals(Fun_GPP), GPP$Tair, use="pairwise.complete.obs"))^2
R2_Preci<-(cor(residuals(Fun_GPP), GPP$Annual_Preci, use="pairwise.complete.obs"))^2

gg3 <- ggplot(data = GPP, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.52", x = 23, y = -800) +
  facet_wrap(~Type_Flux)

gg4 <- ggplot(data = GPP, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.51", x = 3200, y = -800) +
  facet_wrap(~Type_Flux)

# 2.3 Reco

# Compute residual
Reco$Res<- residuals(Fun_Reco)
Reco$Annual_Preci[Reco$Annual_Preci == 0] <- NA

#Plot residuals
R2_Tair<-(cor(residuals(Fun_Reco), Reco$Tair, use="pairwise.complete.obs"))^2
R2_Preci<-(cor(residuals(Fun_Reco), Reco$Annual_Preci, use="pairwise.complete.obs"))^2

gg5 <- ggplot(data = Reco, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.47", x = 23, y = -700) +
  facet_wrap(~Type_Flux)

gg6 <- ggplot(data = Reco, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.59", x = 3200, y = -700) +
  facet_wrap(~Type_Flux)

# 2.4 Ratio GPP-Reco
levels(Ratio_GPP_Reco$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-Reco", "Ratio NEP-GPP", "Ratio NEP-GPPmax")

# Compute residual
Ratio_GPP_Reco$Res<- residuals(Fun_Ratio_GPP_Reco)
Ratio_GPP_Reco$Annual_Preci[Ratio_GPP_Reco$Annual_Preci == 0] <- NA

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
  annotate("text", label = "R-squared = 0.26", x = 23, y = -0.4) +
  facet_wrap(~Type_Flux)

gg8 <- ggplot(data = Ratio_GPP_Reco, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.11", x = 3200, y = -0.4) +
  facet_wrap(~Type_Flux)

# 2.5 Ratio NEP-GPP
levels(Ratio_NEP_GPP$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-Reco", "Ratio NEP-GPP", "Ratio NEP-GPPmax")

# Compute residual
Ratio_NEP_GPP$Res<- residuals(Fun_Ratio_NEP_GPP)
Ratio_NEP_GPP$Annual_Preci[Ratio_NEP_GPP$Annual_Preci == 0] <- NA

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
  annotate("text", label = "R-squared = 0.25", x = 23, y = -0.55) +
  facet_wrap(~Type_Flux)

gg10 <- ggplot(data = Ratio_NEP_GPP, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.12", x = 3200, y = -0.55) +
  facet_wrap(~Type_Flux)

# 2.6 Ratio NEP-GPPclimax

# Compute residual
Ratio_NEP_GPPmax$Res<- residuals(Fun_Ratio_NEP_GPPmax)
Ratio_NEP_GPPmax$Annual_Preci[Ratio_NEP_GPPmax$Annual_Preci == 0] <- NA

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
  annotate("text", label = "R-squared = 0.20", x = 23, y = -2) +
  facet_wrap(~Type_Flux)

gg12 <- ggplot(data = Ratio_NEP_GPPmax, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.16", x = 3200, y = -2) +
  facet_wrap(~Type_Flux)

# 2.6 Plot all graphs together
pdf("Latex/Figures/Residual_Flux.eps", width = 17, height = 14) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, gg10, gg11, gg12, nrow=6) # Write the grid.arrange in the file
dev.off() # Close the file