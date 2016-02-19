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
library(caret)

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
  Fun_Age<-   try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = NEP[NEP$Site_ID != i,], 
                        start = list(A=192.93829, k=-0.08976, offset=-700), control = list(maxiter = 500)), silent=TRUE); 
  Age[[i]]<-  if (inherits(Fun_Age, "nls")) sim = predict(Fun_Age, newdata=NEP[NEP$Site_ID == i,]) else NA;
  Fun_Tair<- try(nlsLM(values~A*Tair^2+B*Tair+C, data = NEP[NEP$Site_ID != i,], 
                       start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500)), silent=TRUE);
  Tair[[i]]<-  if (inherits(Fun_Tair, "nls")) sim = predict(Fun_Tair, newdata=NEP[NEP$Site_ID == i,]) else NA;
  Fun_GPP<-  try(nlsLM(values~offset + A*(1-exp(k*GPP)), data = NEP[NEP$Site_ID != i,], 
                       start = list(A=8.000e+02, k=-2.403e-04, offset=-700), control = list(maxiter = 500)), silent=TRUE);
  GPP[[i]]<- if (inherits(Fun_GPP, "nls")) sim = predict(Fun_GPP, newdata=NEP[NEP$Site_ID == i,]) else NA; 
}

# Append transform climate variables and stand age
NEP$f_Age<- melt(Age)$value
NEP$f_Tair<- melt(Tair)$value
NEP$f_GPP<- melt(GPP)$value

# Check for variable correlation
correlationMatrix <- cor(NEP[c("Stand_Age", "f_Age", "GPP", "Tair", "SPI_CRU", 
                               "Clay_1km", "CEC_Total_1km")])
print(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
print(highlyCorrelated)

# Stepwise regression
lm.NEP<-lm(values ~ (f_Age + GPP + Stand_Age + Tair + SPI_CRU + Clay_1km)^2, data=NEP)
step.NEP<- stepAIC(lm.NEP, direction = "backward",  k = 2)
summary(step.NEP)

# Check for residuals normality and homoskedasticity
res1=residuals(lm.NEP)
shapiro.test(res1)
plot(step.NEP) # ncreasing trend in the scale location plot shows heteroskedasticity

# Assess model performance
NEP$prediction <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "f_Tair", "Stand_Age", "f_Age", "GPP", "f_GPP", "SPI_CRU", "CEC_Total_1km",
                                      "Clay_1km", "MAT_An")]
  lm.NEP<- lm(values~   f_Age + GPP + Stand_Age + Tair + SPI_CRU + 
                Clay_1km + f_Age:GPP + f_Age:Stand_Age + f_Age:Tair + f_Age:SPI_CRU + 
                f_Age:Clay_1km + GPP:Stand_Age + GPP:Tair,
              data=train.df)
  NEP.pred = predict(object = lm.NEP, newdata = test.df)
  NEP$prediction[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$prediction, NEP$values, use="complete")^2
RMSE_NEP <- rmse(NEP$prediction, NEP$values)
NSE_NEP<-NSE(NEP$prediction, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$prediction, NEP$values)

# Compute Importance variable
VarImp_NEP <- boot.relimp(values~ Stand_Age + f_Age + GPP + Tair + SPI_CRU + 
                           Clay_1km + Stand_Age:f_Age + Stand_Age:GPP + Stand_Age:SPI_CRU + 
                           f_Age:GPP + f_Age:Tair + f_Age:SPI_CRU +
                           GPP:Tair + GPP:Clay_1km + Tair:SPI_CRU,
                         data= NEP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(VarImp_NEP))


## Mean site
#-----------------------------------------------------------------
#Compute transform function
id<-unique(NEP_Mean_Site$Site_ID)
Age<- c()
GPP<- c()
Tair<- c()
for (i in id){
  Fun_Age<-    try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = NEP_Mean_Site[NEP_Mean_Site$Site_ID != i,], 
                         start = list(A=192.93829, k=-0.08976, offset=-700), control = list(maxiter = 500)), silent=TRUE); 
  Age[[i]]<-  if (inherits(Fun_Age, "nls")) sim = predict(Fun_Age, newdata=NEP_Mean_Site[NEP_Mean_Site$Site_ID == i,]) else NA;
  Fun_Tair<- try(nlsLM(values~A*Tair^2+B*Tair+C, data = NEP_Mean_Site[NEP_Mean_Site$Site_ID != i,], 
                       start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500)), silent=TRUE);
  Tair[[i]]<-  if (inherits(Fun_Tair, "nls")) sim = predict(Fun_Tair, newdata=NEP_Mean_Site[NEP_Mean_Site$Site_ID == i,]) else NA;
  Fun_GPP<-  try(nlsLM(values~A*GPP^2+B*GPP +C, data = NEP_Mean_Site[NEP_Mean_Site$Site_ID != i,], 
                       start = list(A=-1.077e-04, B= 5.378e-01, C=-2.785e+02), control = list(maxiter = 500)), silent=TRUE);
  GPP[[i]]<- if (inherits(Fun_GPP, "nls")) sim = predict(Fun_GPP, newdata=NEP_Mean_Site[NEP_Mean_Site$Site_ID == i,]) else NA; 
}

# Append transform climate variables and stand age
NEP_Mean_Site$f_Age<- melt(Age)$value
NEP_Mean_Site$f_Tair<- melt(Tair)$value
NEP_Mean_Site$f_GPP<- melt(GPP)$value

# Check for variable correlation
correlationMatrix <- cor(NEP_Mean_Site[c("f_Tair", "Stand_Age", "f_Age", "f_GPP", "SPI_CRU", 
                               "Clay_1km", "CEC_Total_1km")])
print(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
print(highlyCorrelated)

# Stepwise regression
lm.NEP_Mean_Site<-lm(values ~ (f_Age + f_GPP + f_Tair + Stand_Age + SPI_CRU + CEC_Total_1km)^2, data=NEP_Mean_Site)
step.NEP_Mean_Site<- step(lm.NEP_Mean_Site, direction = "backward")
summary(step.NEP_Mean_Site)

# Check for residuals normality and homoskedasticity
res1=residuals(lm.NEP_Mean_Site)
shapiro.test(res1)
plot(step.NEP_Mean_Site) # ncreasing trend in the scale location plot shows heteroskedasticity

# Assess model performance
NEP_Mean_Site$prediction <- NA
for(id in unique(NEP_Mean_Site$Site_ID)){
  train.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID != id,]
  test.df <- NEP_Mean_Site[NEP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "Stand_Age", "f_Age",
                                                          "GPP", "f_GPP", "SPI_CRU", "MAT_An",
                                                          "CEC_Total_1km", "Clay_1km")]
  lm.NEP_Mean_Site<- lm(values ~ f_Age + f_GPP + f_Tair + Stand_Age + SPI_CRU + 
                          CEC_Total_1km + f_Age:Stand_Age + f_Age:SPI_CRU + f_GPP:Stand_Age + 
                          f_GPP:CEC_Total_1km + f_Tair:SPI_CRU + SPI_CRU:CEC_Total_1km, data=train.df)
  NEP_Mean_Site.pred = predict(object = lm.NEP_Mean_Site, newdata = test.df)
  NEP_Mean_Site$prediction[NEP_Mean_Site$Site_ID == id] <- NEP_Mean_Site.pred
}

R2_NEP_Mean_Site<- cor(NEP_Mean_Site$prediction, NEP_Mean_Site$values, use="complete")^2
RMSE_NEP_Mean_Site <- rmse(NEP_Mean_Site$prediction, NEP_Mean_Site$values)
NSE_NEP_Mean_Site<-NSE(NEP_Mean_Site$prediction, NEP_Mean_Site$values, na.rm=TRUE)
Bias_NEP_Mean_Site<-pbias(NEP_Mean_Site$prediction, NEP_Mean_Site$values)

# Compute Importance variable
VarImp_NEP_Mean_Site <- boot.relimp(values~ f_Age + f_GPP + SPI_CRU + f_Age:SPI_CRU, 
                         data= NEP_Mean_Site, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(VarImp_NEP_Mean_Site))

# 1.2 Analysis for Ratio NEP-GPP

## All years per site
#-----------------------------------------------------------------
  
#Compute transform function
id<-unique(Ratio_NEP_GPP$Site_ID)
Age<- c()
GPP<- c()
Tair<- c()
for (i in id){
  Fun_Age<-   try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != i,], 
                        start = list(A= 0.11795, k= -0.03746, offset= -1.5), control = list(maxiter = 500)), silent=TRUE);  
  Age[[i]]<-  if (inherits(Fun_Age, "nls")) sim = predict(Fun_Age, newdata=Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == i,]) else NA;
  Fun_Tair<- try(nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != i,], 
                       start = list(A= -0.001129, B=0.035511, C= -0.080841), control = list(maxiter = 500)), silent=TRUE);
  Tair[[i]]<-  if (inherits(Fun_Tair, "nls")) sim = predict(Fun_Tair, newdata=Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == i,]) else NA;
  Fun_GPP<-   try(nlsLM(values~offset + A*(1-exp(k*GPP)), data = Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != i,], 
                        start = list(A=0.3152035, k=-0.0003991, offset= -1.5), control = list(maxiter = 500)), silent=TRUE);
  GPP[[i]]<- if (inherits(Fun_GPP, "nls")) sim = predict(Fun_GPP, newdata=Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == i,]) else NA; 
}

# Append transform climate variables and stand age
Ratio_NEP_GPP$f_Age<- melt(Age)$value
Ratio_NEP_GPP$f_Tair<- melt(Tair)$value
Ratio_NEP_GPP$f_GPP<- melt(GPP)$value

# Check for variable correlation
correlationMatrix <- cor(Ratio_NEP_GPP[c("Stand_Age", "f_Age", "f_GPP", "f_Tair", "SPI_CRU", 
                               "Clay_1km", "CEC_Total_1km")])
print(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
print(highlyCorrelated)

# Stepwise regression
lm.Ratio_NEP_GPP<-lm(values ~ (f_Age + f_GPP + f_Tair + Stand_Age + SPI_CRU + Clay_1km)^2, data=Ratio_NEP_GPP)
step.Ratio_NEP_GPP<- stepAIC(lm.Ratio_NEP_GPP, direction = "backward")
summary(step.Ratio_NEP_GPP)

# Check for residuals normality and homoskedasticity
res1=residuals(lm.Ratio_NEP_GPP)
shapiro.test(res1)
plot(step.Ratio_NEP_GPP) # ncreasing trend in the scale location plot shows heteroskedasticity

# Assess model performance
Ratio_NEP_GPP$prediction <- NA
for(id in unique(Ratio_NEP_GPP$Site_ID)){
  train.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != id,]
  test.df <- Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "f_Age", "GPP", "f_GPP", "SPI_CRU", "MAT_An",
                                                          "CEC_Total_1km", "Clay_1km")]
  lm.Ratio_NEP_GPP<- lm(values ~ f_Age + f_GPP + f_Tair + Stand_Age + SPI_CRU + 
                          Clay_1km + f_Age:Stand_Age + f_Age:SPI_CRU +
                          f_GPP:f_Tair + f_GPP:Stand_Age + f_Tair:SPI_CRU, data=train.df)
  Ratio_NEP_GPP.pred = predict(object = lm.Ratio_NEP_GPP, newdata = test.df)
  Ratio_NEP_GPP$prediction[Ratio_NEP_GPP$Site_ID == id] <- Ratio_NEP_GPP.pred
}

R2_Ratio_NEP_GPP<- cor(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values, use="complete")^2
RMSE_Ratio_NEP_GPP <- (sum((Ratio_NEP_GPP$prediction-Ratio_NEP_GPP$values)^2)/length(Ratio_NEP_GPP$values))^(1/2)
NSE_Ratio_NEP_GPP<-NSE(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values, na.rm=TRUE)
Bias_Ratio_NEP_GPP<-pbias(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values) 

# Compute Importance variable
VarImp_NEP_GPP <- boot.relimp(values~ f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:SPI_CRU + 
                           f_GPP:f_Tair + f_GPP:SPI_CRU, 
                         data= Ratio_NEP_GPP, 
                         b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(VarImp_NEP_GPP))

## Mean site
-----------------------------------------------------------------

#Compute transform function
id<-unique(Ratio_NEP_GPP_Mean_Site$Site_ID)
Age<- c()
GPP<- c()
Tair<- c()
for (i in id){
  Fun_Age<-   try(nlsLM(values~offset + A*(1+((B*((Stand_Age/C)^D)-1)/(exp(Stand_Age/C)))), data = Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != i,], 
                        start = list(A = 2.949e-04, B = 7.755e+02, C = 3.011e+01, D= 2.273e+00, offset= 1.5), control = list(maxiter = 500)), silent=TRUE);  
  Age[[i]]<-  if (inherits(Fun_Age, "nls")) sim = predict(Fun_Age, newdata=Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == i,]) else NA;
  Fun_Tair<- try(nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != i,], 
                       start = list(A= -0.001129, B=0.035511, C= -0.080841), control = list(maxiter = 500)), silent=TRUE);
  Tair[[i]]<-  if (inherits(Fun_Tair, "nls")) sim = predict(Fun_Tair, newdata=Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == i,]) else NA;
  Fun_GPP<-   try(nlsLM(values~offset + A*(1-exp(k*GPP)), data = Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != i,], 
                        start = list(A=0.3152035, k=-0.0003991, offset= -1.5), control = list(maxiter = 500)), silent=TRUE);
  GPP[[i]]<- if (inherits(Fun_GPP, "nls")) sim = predict(Fun_GPP, newdata=Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == i,]) else NA; 
}

# Append transform climate variables and stand age
Ratio_NEP_GPP_Mean_Site$f_Age<- melt(Age)$value
Ratio_NEP_GPP_Mean_Site$f_Tair<- melt(Tair)$value
Ratio_NEP_GPP_Mean_Site$f_GPP<- melt(GPP)$value

# Check for variable correlation
correlationMatrix <- cor(Ratio_NEP_GPP_Mean_Site[c("Stand_Age", "f_Age", "f_GPP", "f_Tair", "SPI_CRU", 
                                         "Clay_1km", "CEC_Total_1km")])
print(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
print(highlyCorrelated)

# Stepwise regression
lm.Ratio_NEP_GPP_Mean_Site<-lm(values ~ (f_Age + f_GPP + f_Tair + Stand_Age + SPI_CRU + CEC_Total_1km)^2, data=Ratio_NEP_GPP_Mean_Site)
step.Ratio_NEP_GPP_Mean_Site<- stepAIC(lm.Ratio_NEP_GPP_Mean_Site, direction = "backward")
summary(step.Ratio_NEP_GPP_Mean_Site)

# Check for residuals normality and homoskedasticity
res1=residuals(lm.Ratio_NEP_GPP_Mean_Site)
shapiro.test(res1)
plot(step.Ratio_NEP_GPP_Mean_Site) # ncreasing trend in the scale location plot shows heteroskedasticity

# Assess model performance
Ratio_NEP_GPP_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPP_Mean_Site$Site_ID)){
  train.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != id,]
  test.df <- Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "f_Age", "GPP", "f_GPP", "SPI_CRU", "MAT_An",
                                                          "CEC_Total_1km", "Clay_1km")]
  lm.Ratio_NEP_GPP_Mean_Site<- lm(values ~ f_Age + f_GPP + f_Tair + Stand_Age + SPI_CRU + 
                          Clay_1km + f_Age:Stand_Age + f_Age:SPI_CRU +
                          f_GPP:f_Tair + f_GPP:Stand_Age + f_Tair:SPI_CRU, data=train.df)
  Ratio_NEP_GPP_Mean_Site.pred = predict(object = lm.Ratio_NEP_GPP_Mean_Site, newdata = test.df)
  Ratio_NEP_GPP_Mean_Site$prediction[Ratio_NEP_GPP_Mean_Site$Site_ID == id] <- Ratio_NEP_GPP_Mean_Site.pred
}

R2_Ratio_NEP_GPP_Mean_Site<- cor(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values, use="complete")^2
RMSE_Ratio_NEP_GPP_Mean_Site <- (sum((Ratio_NEP_GPP_Mean_Site$prediction-Ratio_NEP_GPP_Mean_Site$values)^2)/length(Ratio_NEP_GPP_Mean_Site$values))^(1/2)
NSE_Ratio_NEP_GPP_Mean_Site<-NSE(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_NEP_GPP_Mean_Site<-pbias(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values) 

# Compute Importance variable
VarImp_NEP_GPP_Mean_Site <- boot.relimp(values~ f_Age + f_GPP + f_Tair + SPI_CRU + f_Age:SPI_CRU + 
                                f_GPP:f_Tair + f_GPP:SPI_CRU, 
                              data= Ratio_NEP_GPP_Mean_Site, 
                              b = 100,  
                              type = "lmg",
                              rank = TRUE, diff = TRUE, rela = T)
print(booteval.relimp(VarImp_NEP_GPP_Mean_Site)) 

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
