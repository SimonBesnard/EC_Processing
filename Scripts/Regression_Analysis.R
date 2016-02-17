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
library (MASS)
library (reshape2)

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

#1.2 Add GPP data to dataframe
NEP$GPP<- GPP$values
NEP_Mean_Site$GPP<- GPP_Mean_Site$values
Ratio_GPP_Reco$GPP<- GPP$values
Ratio_GPP_Reco_Mean_Site$GPP<- GPP_Mean_Site$values
Ratio_NEP_GPP$GPP<- GPP$values
Ratio_NEP_GPP_Mean_Site$GPP<- GPP_Mean_Site$values

# 1.3 Analysis for NEP

## All years per site
#-----------------------------------------------------------------
#Compute transform function
id<-unique(NEP$Site_ID)
Age<- c()
GPP<- c()
Tair<- c()
for (i in id){
  Fun_Age<-   try(nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = NEP[NEP$Site_ID != i,],
                        start = list(A= 2.867e+02, B=-4.941e-03, C=-1.014e+03, D=-1.802e-01), control = list(maxiter = 500)), silent=TRUE); 
  Age[[i]]<-  if (inherits(Fun_Age, "nls")) sim = predict(Fun_Age, newdata=NEP[NEP$Site_ID == i,]) else NA;
  Fun_Tair<- try(nlsLM(values~A*Tair^2+B*Tair+C, data = NEP[NEP$Site_ID != i,], 
                       start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500)), silent=TRUE);
  Tair[[i]]<-  if (inherits(Fun_Tair, "nls")) sim = predict(Fun_Tair, newdata=NEP[NEP$Site_ID == i,]) else NA;
  Fun_GPP<-  try(nlsLM(values~A*(1-exp(k*GPP)), data = NEP[NEP$Site_ID != i,], 
                       start = list(A=8.000e+02, k=-2.403e-04), control = list(maxiter = 500)), silent=TRUE);
  GPP[[i]]<- if (inherits(Fun_GPP, "nls")) sim = predict(Fun_GPP, newdata=NEP[NEP$Site_ID == i,]) else NA; 
}

# Append transform climate variables and stand age
NEP$f_Age<- melt(Age)$value
NEP$f_Tair<- melt(Tair)$value
NEP$f_GPP<- melt(GPP)$value

# Perform PCA to check correlation between predictors
PCA_NEP<- NEP[c("f_Age", "GPP", "f_Tair", "Annual_Preci", "SPI_CRU", "CEC_Total_1km", "Clay_1km", "Stand_Age")]
data.pca<-dudi.pca(PCA_NEP,center=T,scale=T,scannf=T,nf=5) #calling the PCA
data.pca$eig
cumsum(data.pca$eig)/sum(data.pca$eig)
data.pca$co #correlation table between PCA axes and original variables
s.corcircle(data.pca$co,xax=1,yax=3)
s.corcircle(data.pca$co,xax=1,yax=4)
s.label(data.pca$li,xax=1,yax=2) #plotting PCA results
s.label(data.pca$li,xax=1,yax=3)
s.label(data.pca$li,xax=1,yax=4)
scatter(data.pca) #biplot graph
s.corcircle(data.pca$co,xax=1,yax=2) #plotting correlation circles

# Stepwise regression
lm.NEP<-lm(values ~ (f_Age + f_GPP + f_Tair + SPI_CRU)^2, data=NEP)
step.NEP<- stepAIC(lm.NEP, direction = "backward")
summary(step.NEP)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:SPI_CRU + 
                           f_GPP:f_Tair,
                         data= NEP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
NEP$prediction <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "f_Tair", "f_Age", "f_GPP", "SPI_CRU", "CEC_Total_1km", "Clay_1km", "MAT_An")]
  lm.NEP<- lm(values~ f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:SPI_CRU + 
                f_GPP:f_Tair, data=train.df)
  NEP.pred = predict(object = lm.NEP, newdata = test.df)
  NEP$prediction[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$prediction, NEP$values, use="complete")^2
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
f_Age_NEP_Mean_Site<- function (x) {2.548650e+02*(exp(-2.543848e-03*x)) -1.006867e+03*(exp(-2.173415e-01*x))}

#Tair
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = NEP_Mean_Site, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_NEP_Mean_Site<- function (x) {-0.5140868*x^2+24.5816118*x+38.0432864}

# GPP
Fun_NEP_Photo<-nlsLM(values~A*GPP^2+B*GPP+C, data = NEP_Mean_Site, 
                     start = list(A=-1.077e-04, B= 5.378e-01, C=-2.785e+02), control = list(maxiter = 500))
coef(Fun_NEP_Photo)
f_Photo_NEP<- function (x) {-1.117996e-04*x^2+5.623442e-01*x -2.832764e+02}

# Append transform climate variables and stand age
NEP_Mean_Site$f_Tair<- f_Tair_NEP_Mean_Site(NEP_Mean_Site$Tair)
NEP_Mean_Site$f_Age<- f_Age_NEP_Mean_Site(NEP_Mean_Site$Stand_Age)
NEP_Mean_Site$f_GPP<- f_Photo_NEP(NEP_Mean_Site$GPP)

# Stepwise regression
lm.NEP_Mean_Site<-lm(values ~ (f_Age + f_GPP + f_Tair + SPI_CRU)^2, data=NEP_Mean_Site)
step.NEP_Mean_Site<- step(lm.NEP_Mean_Site, direction = "backward")
summary(step.NEP_Mean_Site)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_Age + f_GPP + SPI_CRU + f_Age:SPI_CRU, 
                         data= NEP_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
NEP_Mean_Site$prediction <- NA
for(id in unique(NEP_Mean_Site$Site_ID)){
  train.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID != id,]
  test.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "f_Age", "f_GPP", "SPI_CRU", "MAT_An",
                                                          "CEC_Total_1km")]
  lm.NEP_Mean_Site<- lm(values ~  f_Age + f_GPP + SPI_CRU + f_Age:SPI_CRU, data=train.df)
  NEP_Mean_Site.pred = predict(object = lm.NEP_Mean_Site, newdata = test.df)
  NEP_Mean_Site$prediction[NEP_Mean_Site$Site_ID == id] <- NEP_Mean_Site.pred
}

R2_NEP_Mean_Site<- cor(NEP_Mean_Site$prediction, NEP_Mean_Site$values, use="complete")^2
RMSE_NEP_Mean_Site <- rmse(NEP_Mean_Site$prediction, NEP_Mean_Site$values)
NSE_NEP_Mean_Site<-NSE(NEP_Mean_Site$prediction, NEP_Mean_Site$values, na.rm=TRUE)
Bias_NEP_Mean_Site<-pbias(NEP_Mean_Site$prediction, NEP_Mean_Site$values)

# 1.2 Analysis for Ratio GPP-Reco

## All years per site
#-----------------------------------------------------------------
# Compute transform function

#Age
Fun_Ratio_GPP_Reco<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = Ratio_GPP_Reco, 
                         start = list(A= 1.1582, k= -0.2312), control = list(maxiter = 500))
coef(Fun_Ratio_GPP_Reco)
f_Age_Ratio_GPP_Reco<- function (x) {1.2134719*(1-exp(-0.2257333*x))}

#Tair
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_GPP_Reco, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_GPP_Reco<- function (x) {-0.00125467*x^2+0.03968507*x+ 0.99669771}

# GPP
Fun_GPP_ER_Photo<-nlsLM(values~A*(GPP^B)*(exp(k*GPP)), data = Ratio_GPP_Reco, 
                        start = list(A = 0.0662924, B =0.4576941, k =-0.0002671), control = list(maxiter = 500))
coef(Fun_GPP_ER_Photo)
f_Photo_Ratio_GPP_Reco<- function (x) {0.0591990517*(x^0.4768037616)*(exp(-0.0002775358*x))}

# Append transform climate variables and stand age
Ratio_GPP_Reco$f_Tair<- f_Tair_Ratio_GPP_Reco(Ratio_GPP_Reco$Tair)
Ratio_GPP_Reco$f_Age<- f_Age_Ratio_GPP_Reco(Ratio_GPP_Reco$Stand_Age)
Ratio_GPP_Reco$f_GPP<- f_Photo_Ratio_GPP_Reco(Ratio_GPP_Reco$GPP)

# Stepwise regression
lm.Ratio_GPP_Reco<-lm(values ~ (f_Age + f_GPP + f_Tair + SPI_CRU)^2, data=Ratio_GPP_Reco)
step.Ratio_GPP_Reco<- stepAIC(lm.Ratio_GPP_Reco, direction = "backward")
summary(step.Ratio_GPP_Reco)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:SPI_CRU + 
                           f_Tair:SPI_CRU, 
                         data= Ratio_GPP_Reco, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_GPP_Reco$prediction <- NA
for(id in unique(Ratio_GPP_Reco$Site_ID)){
  train.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID != id,]
  test.df <- Ratio_GPP_Reco[Ratio_GPP_Reco$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "f_Age", "f_GPP", "SPI_CRU",
                                                            "MAT_An", "CEC_Total_1km")]
  lm.Ratio_GPP_Reco<- lm(values ~  f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:SPI_CRU + 
                           f_Tair:SPI_CRU, data=train.df)
  Ratio_GPP_Reco.pred = predict(object = lm.Ratio_GPP_Reco, newdata = test.df)
  Ratio_GPP_Reco$prediction[Ratio_GPP_Reco$Site_ID == id] <- Ratio_GPP_Reco.pred
}

R2_Ratio_GPP_Reco<- cor(Ratio_GPP_Reco$prediction, Ratio_GPP_Reco$values, use="complete")^2
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
f_Age_Ratio_GPP_Reco_Mean_Site<- function (x) {1.205008*(1-exp(-0.243221*x))}

#Tair
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_GPP_Reco_Mean_Site, 
                start = list(A= -0.001214, B=0.036857, C= 0.979401), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_GPP_Reco_Mean_Site<- function (x) {-0.0007484326*x^2+ 0.0259862124*x + 1.0399136528}

# GPP
Fun_GPP_ER_Photo<-nlsLM(values~A*(GPP^B)*(exp(k*GPP)), data = Ratio_GPP_Reco_Mean_Site, 
                        start = list(A = 0.0662924, B =0.4576941, k =-0.0002671), control = list(maxiter = 500))
coef(Fun_GPP_ER_Photo)
f_Photo_Ratio_GPP_Reco<- function (x) {0.0451066043*(x^0.5238655238)*(exp(-0.0003126905*x))}

# Append transformed variables
Ratio_GPP_Reco_Mean_Site$f_Tair<- f_Tair_Ratio_GPP_Reco_Mean_Site(Ratio_GPP_Reco_Mean_Site$Tair)
Ratio_GPP_Reco_Mean_Site$f_Age<- f_Age_Ratio_GPP_Reco_Mean_Site(Ratio_GPP_Reco_Mean_Site$Stand_Age)
Ratio_GPP_Reco_Mean_Site$f_GPP<- f_Photo_Ratio_GPP_Reco(Ratio_GPP_Reco_Mean_Site$GPP)

# Stepwise regression
lm.Ratio_GPP_Reco_Mean_Site<-lm(values ~ (f_Age + f_GPP + f_Tair + SPI_CRU)^2, data=Ratio_GPP_Reco_Mean_Site)
step.Ratio_GPP_Reco_Mean_Site<- stepAIC(lm.Ratio_GPP_Reco_Mean_Site, direction = "backward")
summary(step.Ratio_GPP_Reco_Mean_Site)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_Age + f_GPP + f_Tair + SPI_CRU + f_GPP:f_Tair,  
                         data= Ratio_GPP_Reco_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_GPP_Reco_Mean_Site$prediction <- NA
for(id in unique(Ratio_GPP_Reco_Mean_Site$Site_ID)){
  train.df <- Ratio_GPP_Reco_Mean_Site[Ratio_GPP_Reco_Mean_Site$Site_ID != id,]
  test.df <- Ratio_GPP_Reco_Mean_Site[Ratio_GPP_Reco_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "f_Age", "f_GPP",
                                                                                "SPI_CRU", "MAT_An", "CEC_Total_1km")]
  lm.Ratio_GPP_Reco_Mean_Site<- lm(values ~f_Age + f_GPP + f_Tair + SPI_CRU + f_GPP:f_Tair,  data=train.df)
  Ratio_GPP_Reco_Mean_Site.pred = predict(object = lm.Ratio_GPP_Reco_Mean_Site, newdata = test.df)
  Ratio_GPP_Reco_Mean_Site$prediction[Ratio_GPP_Reco_Mean_Site$Site_ID == id] <- Ratio_GPP_Reco_Mean_Site.pred
}

R2_Ratio_GPP_Reco_Mean_Site<- cor(Ratio_GPP_Reco_Mean_Site$prediction, Ratio_GPP_Reco_Mean_Site$values, use= "complete")^2
RMSE_Ratio_GPP_Reco_Mean_Site <- (sum((Ratio_GPP_Reco_Mean_Site$prediction-Ratio_GPP_Reco_Mean_Site$values)^2)/length(Ratio_GPP_Reco_Mean_Site$values))^(1/2)
NSE_Ratio_GPP_Reco_Mean_Site<-NSE(Ratio_GPP_Reco_Mean_Site$prediction, Ratio_GPP_Reco_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_GPP_Reco_Mean_Site<-pbias(Ratio_GPP_Reco_Mean_Site$prediction, Ratio_GPP_Reco_Mean_Site$values)  

# 1.3 Analysis for Ratio NEP-GPP

## All years per site
#-----------------------------------------------------------------
  
# Compute transform function

#Age
Fun_Ratio_NEP_GPP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPP,
                         start = list(A=0.165450, B= -0.003772, C=-1.319022, D=-0.148503), control = list(maxiter = 500))
coef(Fun_Ratio_NEP_GPP)
f_Age_Ratio_NEP_GPP<- function (x) {0.225119297*(exp(-0.005007755*x)) -1.530424388*(exp(-0.170444337*x))}

#Tair
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPP, 
                start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPP<- function (x) {-0.001135343*x^2+0.036432176*x-0.068136762}

# GPP
Fun_NEP_GPP_Photo<- nlsLM(values~A*(1+((B*((GPP/C)^D)-1)/(exp(GPP/C)))), data = Ratio_NEP_GPP, 
                          start = list(A = 0.188, B =  -4.871, C =  365.122, D=  -0.393), control = list(maxiter = 500))
coef(Fun_NEP_GPP_Photo)
f_Photo_Ratio_NEP_GPP<- function (x) {0.1982498*(1+((-3.4772396*((x/426.2889796)^ -0.5526806)-1)/(exp(x/426.2889796))))}

# Append transform climate variables and stand age
Ratio_NEP_GPP$f_Tair<- f_Tair_Ratio_NEP_GPP(Ratio_NEP_GPP$Tair)
Ratio_NEP_GPP$f_Age<- f_Age_Ratio_NEP_GPP(Ratio_NEP_GPP$Stand_Age)
Ratio_NEP_GPP$f_GPP<- f_Photo_Ratio_NEP_GPP(Ratio_NEP_GPP$GPP)

# Perform PCA analysis to check variable correlation
PCA_Ratio_NEP_GPP<- Ratio_NEP_GPP[c("MAT_An", "f_Age", "f_GPP", "f_Tair", "Annual_Preci", "SPI_CRU", "CEC_Total_1km", "Clay_1km")]
data.pca<-dudi.pca(PCA_Ratio_NEP_GPP,center=T,scale=T,scannf=T,nf=5) #calling the PCA
data.pca$eig
cumsum(data.pca$eig)/sum(data.pca$eig)
data.pca$co #correlation table between PCA axes and original variables
s.corcircle(data.pca$co,xax=1,yax=3)
s.corcircle(data.pca$co,xax=1,yax=4)
s.label(data.pca$li,xax=1,yax=2) #plotting PCA results
s.label(data.pca$li,xax=1,yax=3)
s.label(data.pca$li,xax=1,yax=4)
scatter(data.pca) #biplot graph
s.corcircle(data.pca$co,xax=1,yax=2) #plotting correlation circles

# Stepwise regression
lm.Ratio_NEP_GPP<-lm(values ~ (f_Age + f_GPP + f_Tair + SPI_CRU)^2, data=Ratio_NEP_GPP)
step.Ratio_NEP_GPP<- stepAIC(lm.Ratio_NEP_GPP, direction = "backward")
summary(step.Ratio_NEP_GPP)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:SPI_CRU + 
                           f_GPP:f_Tair + f_GPP:SPI_CRU, 
                         data= Ratio_NEP_GPP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPP$prediction <- NA
for(id in unique(Ratio_NEP_GPP$Site_ID)){
  train.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != id,]
  test.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "f_Age", "f_GPP", "SPI_CRU", "MAT_An",
                                                          "CEC_Total_1km", "Clay_1km")]
  lm.Ratio_NEP_GPP<- lm(values ~ f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:SPI_CRU + 
                          f_GPP:f_Tair + f_GPP:SPI_CRU, data=train.df)
  Ratio_NEP_GPP.pred = predict(object = lm.Ratio_NEP_GPP, newdata = test.df)
  Ratio_NEP_GPP$prediction[Ratio_NEP_GPP$Site_ID == id] <- Ratio_NEP_GPP.pred
}

R2_Ratio_NEP_GPP<- cor(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values, use="complete")^2
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
f_Age_Ratio_NEP_GPP_Mean_Site<- function (x) {0.1497866064*(exp(-0.0004478665*x)) -1.4211923834*(exp(-0.1848797631*x))}

#Tair
Fun_Tair<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPP_Mean_Site, 
                start = list(A= -0.001129, B=0.035511, C= -0.080841), control = list(maxiter = 500))
coef(Fun_Tair)
f_Tair_Ratio_NEP_GPP_Mean_Site<- function (x) {-0.0006538801*x^2+0.0253992259*x-0.0582247503}

# GPP
f_Photo_Ratio_NEP_GPP<-nlsLM(values~A*GPP^2+B*GPP +C, data = Ratio_NEP_GPP_Mean_Site, 
                                   start = list(A=-1.337e-07, B= 5.752e-04, C=-3.631e-01), control = list(maxiter = 500))
coef(f_Photo_Ratio_NEP_GPP)
f_Photo_Ratio_NEP_GPP<- function(x) {-1.553290e-07*x^2+6.838167e-04*x -4.443732e-01}

# Append transform climate variables and stand age
Ratio_NEP_GPP_Mean_Site$f_Tair<- f_Tair_Ratio_NEP_GPP_Mean_Site(Ratio_NEP_GPP_Mean_Site$Tair)
Ratio_NEP_GPP_Mean_Site$f_Age<- f_Age_Ratio_NEP_GPP_Mean_Site(Ratio_NEP_GPP_Mean_Site$Stand_Age)
Ratio_NEP_GPP_Mean_Site$f_GPP<- f_Photo_Ratio_NEP_GPP(Ratio_NEP_GPP_Mean_Site$GPP)

# Stepwise regression
lm.Ratio_NEP_GPP_Mean_Site<-lm(values ~ (f_Age + f_GPP + f_Tair + SPI_CRU)^2, data=Ratio_NEP_GPP_Mean_Site)
step.Ratio_NEP_GPP_Mean_Site<- stepAIC(lm.Ratio_NEP_GPP_Mean_Site, direction = "backward")
summary(step.Ratio_NEP_GPP_Mean_Site)

# Compute Importance variable
bootswiss <- boot.relimp(values~ f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:f_GPP + 
                           f_Age:SPI_CRU + f_GPP:f_Tair,
                         data= Ratio_NEP_GPP_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(bootswiss))

# Assess model performance
Ratio_NEP_GPP_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPP_Mean_Site$Site_ID)){
  train.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != id,]
  test.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "f_Age", "f_GPP",
                                                                              "SPI_CRU", "MAT_An", "CEC_Total_1km")]
  lm.Ratio_NEP_GPP_Mean_Site<- lm(values ~ f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:f_GPP + 
                                    f_Age:SPI_CRU + f_GPP:f_Tair, data=train.df)
  Ratio_NEP_GPP_Mean_Site.pred = predict(object = lm.Ratio_NEP_GPP_Mean_Site, newdata = test.df)
  Ratio_NEP_GPP_Mean_Site$prediction[Ratio_NEP_GPP_Mean_Site$Site_ID == id] <- Ratio_NEP_GPP_Mean_Site.pred
}

R2_Ratio_NEP_GPP_Mean_Site<- cor(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values, use="complete")^2
RMSE_Ratio_NEP_GPP_Mean_Site <- (sum((Ratio_NEP_GPP_Mean_Site$prediction-Ratio_NEP_GPP_Mean_Site$values)^2)/length(Ratio_NEP_GPP_Mean_Site$values))^(1/2)
NSE_Ratio_NEP_GPP_Mean_Site<-NSE(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_NEP_GPP_Mean_Site<-pbias(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values)  

# 2. Plot observed vs. actual for the different flux

# 2.1 All years per site

# Prepare dataset for plotting
pred_NEP<- NEP[c("Type_Flux", "values", "prediction", "Stand_Age")]
pred_Ratio_GPP_Reco<- Ratio_GPP_Reco[c("Type_Flux", "values", "prediction", "Stand_Age")]
pred_Ratio_NEP_GPP<- Ratio_NEP_GPP[c("Type_Flux", "values", "prediction", "Stand_Age")]
pred_All<- rbind(pred_NEP, pred_Ratio_GPP_Reco, pred_Ratio_NEP_GPP)
levels(pred_All$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-ER", "Ratio NEP-GPP", 
                                "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4")
df_NEP<-pred_All[pred_All$Type_Flux %in% c("NEP"),]
df_NEP_GPP<-pred_All[pred_All$Type_Flux %in% c("Ratio NEP-GPP"),]
df_GPP_Reco<-pred_All[pred_All$Type_Flux %in% c("Ratio GPP-ER"),]

#Plot prediction vs. observation

# NEP
gg1<- ggplot(df_NEP, aes(x=prediction, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-700, 800)+
  ylim(-700, 800)+
  annotate("text", label = "R-squared = 0.49", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 178.00 gC.m-2.y-1", x = -370, y = 550, size =4)

#Ratio GPP-Reco
gg2<- ggplot(df_GPP_Reco, aes(x=prediction, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(0, 2)+
  ylim(0, 2)+
  annotate("text", label = "R-squared = 0.43", x = 0.30, y = 1.9, size =4) +
  annotate("text", label = "RMSE = 0.19", x = 0.22, y = 1.70, size =4)

# Ratio NEP-GPP
gg3<- ggplot(df_NEP_GPP, aes(x=prediction, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-1.5, 1)+
  ylim(-1.5, 1)+
  annotate("text", label = "R-squared = 0.72", x = -1.1, y = 0.8, size =4) +
  annotate("text", label = "RMSE = 0.14", x = -1.20, y = 0.55, size =4)

# Create an arrange plot object
pdf("Latex/Figures/Pred_Flux_All_Site.eps", width = 5, height = 12) # Open a new pdf file
source("Function/Legend_Grid_Arrange.R")
grid_arrange_shared_legend(gg1, gg2, gg3, nrow=4) # Write the grid.arrange in the file

# 2.2 Average site

# Prepare dataset for plotting
pred_NEP_Mean_Site<- NEP_Mean_Site[c("Type_Flux", "values", "prediction", "Stand_Age")]
pred_Ratio_GPP_Reco_Mean_Site<- Ratio_GPP_Reco_Mean_Site[c("Type_Flux", "values", "prediction", "Stand_Age")]
pred_Ratio_NEP_GPP_Mean_Site<- Ratio_NEP_GPP_Mean_Site[c("Type_Flux", "values", "prediction", "Stand_Age")]
pred_All<- rbind(pred_NEP_Mean_Site, pred_Ratio_GPP_Reco_Mean_Site, 
                 pred_Ratio_NEP_GPP_Mean_Site)
levels(pred_All$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-ER", "Ratio NEP-GPP", 
                                "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4")

df_NEP_Mean_Site<-pred_All[pred_All$Type_Flux %in% c("NEP"),]
df_NEP_GPP_Mean_Site<-pred_All[pred_All$Type_Flux %in% c("Ratio NEP-GPP"),]
df_GPP_Reco_Mean_site<-pred_All[pred_All$Type_Flux %in% c("Ratio GPP-ER"),]

#Plot prediction vs. observation

# NEP
gg1<- ggplot(df_NEP_Mean_Site, aes(x=prediction, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-700, 800)+
  ylim(-700, 800)+
  annotate("text", label = "R-squared = 0.60", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 160.8 gC.m-2.y-1", x = -370, y = 550, size =4)

#Ratio GPP-Reco
gg2<- ggplot(df_GPP_Reco_Mean_site, aes(x=prediction, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(0, 2)+
  ylim(0, 2)+
  annotate("text", label = "R-squared = 0.51", x = 0.30, y = 1.9, size =4) +
  annotate("text", label = "RMSE = 0.20", x = 0.22, y = 1.70, size =4)

# Ratio NEP-GPP
gg3<- ggplot(df_NEP_GPP_Mean_Site, aes(x=prediction, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-1.5, 1)+
  ylim(-1.5, 1)+
  annotate("text", label = "R-squared = 0.73", x = -1.1, y = 0.8, size =4) +
  annotate("text", label = "RMSE = 0.16", x = -1.20, y = 0.55, size =4)

#Plot all plots together
pdf("Latex/Figures/Pred_Flux_Mean_Site.eps", width = 5, height = 12) # Open a new pdf file
source("Function/Legend_Grid_Arrange.R")
grid_arrange_shared_legend(gg1, gg2, gg3, nrow=4) # Write the grid.arrange in the file
