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
library (rbokeh)

#1 Explain variabilty of the fluxes using linear regression analysis

# 1.1 Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
NEP<-readRDS("Output/NEP.rds")
NEP_Mean_Site<-readRDS("Output/NEP_Mean_Site.rds")
GPP<-readRDS("Output/GPP.rds")
GPP_Mean_Site<-readRDS("Output/GPP_Mean_Site.rds")
Ratio_NEP_GPP<-readRDS("Output/Ratio_NEP_GPP.rds")
Ratio_NEP_GPP_Mean_Site<- readRDS("Output/Ratio_NEP_GPP_Mean_Site.rds")

#1.2 Add GPP data to dataframe
NEP$GPP<- GPP$values
NEP_Mean_Site$GPP<- GPP_Mean_Site$values
Ratio_NEP_GPP$GPP<- GPP$values
Ratio_NEP_GPP_Mean_Site$GPP<- GPP_Mean_Site$values

# 1.3 Analysis for NEP

## All years per site
#-----------------------------------------------------------------

#Compute stand age transform function

# weighted
st.w <- list(A = 800, k = -0.2, offset = -700)
obj.w <- function(p, NEP) with(NEP, sum((values - rhs(p, NEP))^2/Uncert^2))
fit.optim.w <- optim(st.w, obj.w, NEP = NEP)
fit.optim.w$par; fit.optim.w$value # 3

fit.nlsLM.w <- nlsLM(values ~ rhs(c(A, k, offset), NEP), data = NEP, weights = 1/NEP$Uncert^2, 
                     start = st.w)
coef(fit.nlsLM.w); deviance(fit.nlsLM.w) # 4 = same as 3


NEP$f_Age<-NA
for(id in unique(NEP$Site_ID)){
  lm.Age<- try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = NEP[NEP$Site_ID != id,], 
                     start = list(A=192.93829, k=-0.08976, offset=-700), control = list(maxiter = 500), weights = 1/NEP[NEP$Site_ID != id,]$Uncert^2), silent=TRUE); 
  Age.pred = predict(object = lm.Age, newdata = NEP[NEP$Site_ID == id,])
  NEP$f_Age[NEP$Site_ID == id] <- Age.pred
}

# Check for variable correlation
correlationMatrix <- cor(NEP[c("Stand_Age", "f_Age", "GPP", "MAT_CRU", "SPI_CRU", 
                               "CEC_Total_1km", "MAT_An", "Trend_An", "NHx", "Soil_C_1km")], use="complete")
print(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.60)
print(highlyCorrelated)

# Stepwise regression
lm.NEP<-lm(values ~ (f_Age + I(GPP)^2 + Stand_Age + I(MAT_CRU)^2 + SPI_CRU + NHx)^2, data=NEP)
step.NEP<- stepAIC(lm.NEP, direction = "backward",  k = 2)
summary(step.NEP)

# Check for residuals normality and homoskedasticity
res1=residuals(lm.NEP)
shapiro.test(res1)
plot(step.NEP) # ncreasing trend in the scale location plot shows heteroskedasticity

# Assess model performance
NEP$prediction <- NA
for(id in unique(NEP$Site_ID)){
  lm.NEP<- lm(values~ f_Age + I(GPP) + Stand_Age + I(MAT_CRU) + 
                SPI_CRU + NHx + f_Age:Stand_Age + f_Age:SPI_CRU + 
                f_Age:NHx + I(GPP):Stand_Age + I(GPP):I(MAT_CRU) + 
                I(GPP):NHx + SPI_CRU:NHx,
              data=NEP[NEP$Site_ID != id,])
  NEP.pred = predict(object = lm.NEP, newdata = NEP[NEP$Site_ID == id,])
  NEP$prediction[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$prediction, NEP$values, use="complete")^2
RMSE_NEP <- rmse(NEP$prediction, NEP$values)
NSE_NEP<-NSE(NEP$prediction, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$prediction, NEP$values)

# Compute Importance variable
#  t-statistic 
VarImp_NEP<- varImp(lm(values~  f_Age + I(GPP) + Stand_Age + I(MAT_CRU) + 
                         SPI_CRU + NHx + f_Age:I(GPP) + f_Age:Stand_Age + f_Age:SPI_CRU + 
                         f_Age:NHx + I(GPP):Stand_Age + I(GPP):I(MAT_CRU) + I(GPP):SPI_CRU + 
                         I(GPP):NHx + Stand_Age:SPI_CRU + SPI_CRU:NHx,
                        data=NEP), scale=T)

#LMG method
VarImp_NEP<- boot.relimp(values ~ f_Age + I(GPP) + Stand_Age + MAT_CRU + 
                           SPI_CRU + NHx + f_Age:Stand_Age + f_Age:SPI_CRU + 
                           f_Age:NHx +  I(GPP):MAT_CRU + 
                           I(GPP):NHx + SPI_CRU:NHx, 
                         data=NEP,
                         b = 100, type = "lmg", rank = T, diff = F, rela = TRUE)
print(booteval.relimp(VarImp_NEP))

## Mean site
#-----------------------------------------------------------------
#Compute transform function
NEP_Mean_Site$f_Age<-NA
for(id in unique(NEP_Mean_Site$Site_ID)){
  lm.Age<-try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = NEP_Mean_Site[NEP_Mean_Site$Site_ID != i,], 
                         start = list(A=192.93829, k=-0.08976, offset=-700), control = list(maxiter = 500)), silent=TRUE); 
  Age.pred = predict(object = lm.Age, newdata = NEP_Mean_Site[NEP_Mean_Site$Site_ID == id,])
  NEP_Mean_Site$f_Age[NEP_Mean_Site$Site_ID == id] <- Age.pred
}


# Check for variable correlation
correlationMatrix <- cor(NEP_Mean_Site[c("Stand_Age", "f_Age", "f_GPP", "f_Tair", "MAT_CRU", "SPI_CRU", 
                                         "Clay_Silt", "CEC_Total_1km")], use="complete")
print(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
print(highlyCorrelated)

# Stepwise regression
lm.NEP_Mean_Site<-lm(values ~ (f_Age + GPP + Stand_Age + MAT_CRU + SPI_CRU + Clay_Silt + Trend_An)^2, data=NEP_Mean_Site)
step.NEP_Mean_Site<- step(lm.NEP_Mean_Site, direction = "backward")
summary(step.NEP_Mean_Site)

# Check for residuals normality and homoskedasticity
res1=residuals(lm.NEP_Mean_Site)
shapiro.test(res1)
plot(step.NEP_Mean_Site) # ncreasing trend in the scale location plot shows heteroskedasticity

# Assess model performance
NEP_Mean_Site$prediction <- NA
for(id in unique(NEP_Mean_Site$Site_ID)){
  lm.NEP_Mean_Site<- lm(values~ GPP + Stand_Age + SPI_CRU + 
                          Trend_An + f_Age:SPI_CRU + 
                          GPP:SPI_CRU + GPP:Clay_Silt + GPP:Trend_An + 
                          Stand_Age:Trend_An + 
                          Clay_Silt:Trend_An,
              data=NEP_Mean_Site[NEP_Mean_Site$Site_ID != id,])
  NEP.Mean_Site.pred = predict(object = lm.NEP_Mean_Site, newdata = NEP_Mean_Site[NEP_Mean_Site$Site_ID == id,])
  NEP_Mean_Site$prediction[NEP_Mean_Site$Site_ID == id] <- NEP.Mean_Site.pred
}

R2_NEP_Mean_Site<- cor(NEP_Mean_Site$prediction, NEP_Mean_Site$values, use="complete")^2
RMSE_NEP_Mean_Site <- rmse(NEP_Mean_Site$prediction, NEP_Mean_Site$values)
NSE_NEP_Mean_Site<-NSE(NEP_Mean_Site$prediction, NEP_Mean_Site$values, na.rm=TRUE)
Bias_NEP_Mean_Site<-pbias(NEP_Mean_Site$prediction, NEP_Mean_Site$values)

# Compute Importance variable
#t-statistic
VarImp_NEP_Mean_Site <-varImp(lm(values~f_Age + GPP + Stand_Age + MAT_CRU + SPI_CRU + 
                                   Clay_Silt + Trend_An + f_Age:GPP + f_Age:MAT_CRU + f_Age:SPI_CRU + 
                                   f_Age:Trend_An + GPP:SPI_CRU + GPP:Clay_Silt + GPP:Trend_An + 
                                   Stand_Age:MAT_CRU + Stand_Age:SPI_CRU + Stand_Age:Trend_An + 
                                   SPI_CRU:Trend_An + Clay_Silt:Trend_An,
                                 data=NEP_Mean_Site))

#LMG method
VarImp_NEP_Mean_Site<- boot.relimp(values~  GPP + Stand_Age + f_Tair + SPI_CRU + 
                                     f_Age:Stand_Age + f_Age:SPI_CRU + 
                                     f_Age:Trend_An + GPP:Stand_Age + GPP:f_Tair + 
                                     GPP:Trend_An+SPI_CRU:Trend_An, 
                         data=NEP_Mean_Site,
                         b = 100, type = "lmg", rank = T, diff = F, rela = TRUE)
print(VarImp_NEP_Mean_Site)

# 1.2 Analysis for Ratio NEP-GPP

## All years per site
#-----------------------------------------------------------------
  
#Compute transform function
Ratio_NEP_GPP$prediction <- NA
for(id in unique(Ratio_NEP_GPP$Site_ID)){
  lm.Age<-   try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != i,], 
                        start = list(A= 0.11795, k= -0.03746, offset= -1.5), control = list(maxiter = 500)), silent=TRUE);  
  Age.pred = predict(object = lm.Age, newdata = Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == id,])
  Ratio_NEP_GPP$f_Age[Ratio_NEP_GPP$Site_ID == id] <- Age.pred
}

# Check for variable correlation
correlationMatrix <- cor(Ratio_NEP_GPP[c("Stand_Age", "f_Age", "f_GPP", "f_Tair", "MAT_CRU", "SPI_CRU", 
                                         "Clay_Silt", "CEC_Total_1km", "Soil_Quality")], use="complete")
print(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
print(highlyCorrelated)

# Stepwise regression
lm.Ratio_NEP_GPP<-lm(values ~ (f_Age + I(GPP)^2 + Stand_Age + I(MAT_CRU)^2 + SPI_CRU + NHx)^2, data=Ratio_NEP_GPP)
step.Ratio_NEP_GPP<- stepAIC(lm.Ratio_NEP_GPP, direction = "backward")
summary(step.Ratio_NEP_GPP)

# Check for residuals normality and homoskedasticity
res1=residuals(lm.Ratio_NEP_GPP)
shapiro.test(res1)
plot(step.Ratio_NEP_GPP) # ncreasing trend in the scale location plot shows heteroskedasticity

# Assess model performance
Ratio_NEP_GPP$prediction <- NA
for(id in unique(Ratio_NEP_GPP$Site_ID)){
  lm.Ratio_NEP_GPP<- lm(values ~ f_Age + I(GPP) + I(MAT_CRU) + 
                          SPI_CRU + NHx + f_Age:I(GPP) + f_Age:SPI_CRU + f_Age:NHx + 
                          I(GPP):I(MAT_CRU) + Stand_Age:I(MAT_CRU), data=Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID != id,])
  Ratio_NEP_GPP.pred = predict(object = lm.Ratio_NEP_GPP, newdata = Ratio_NEP_GPP[Ratio_NEP_GPP$Site_ID == id,])
  Ratio_NEP_GPP$prediction[Ratio_NEP_GPP$Site_ID == id] <- Ratio_NEP_GPP.pred
}

R2_Ratio_NEP_GPP<- cor(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values, use="complete")^2
RMSE_Ratio_NEP_GPP <- (sum((Ratio_NEP_GPP$prediction-Ratio_NEP_GPP$values)^2)/length(Ratio_NEP_GPP$values))^(1/2)
NSE_Ratio_NEP_GPP<-NSE(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values, na.rm=TRUE)
Bias_Ratio_NEP_GPP<-pbias(Ratio_NEP_GPP$prediction, Ratio_NEP_GPP$values) 

# Compute Importance variable
# t-statistic
VarImp_NEP_GPP <- varImp(lm(values~f_Age + I(GPP) + I(MAT_CRU) + 
                              SPI_CRU + NHx + f_Age:I(GPP) + f_Age:SPI_CRU + f_Age:NHx + 
                              I(GPP):I(MAT_CRU) + Stand_Age:I(MAT_CRU),
                            data=Ratio_NEP_GPP))

#LMG method
VarImp_NEP_GPP<- boot.relimp(values~  f_Age + I(GPP) + I(MAT_CRU) + 
                               SPI_CRU + NHx + f_Age:I(GPP) + f_Age:SPI_CRU + f_Age:NHx + 
                               I(GPP):I(MAT_CRU) + Stand_Age:I(MAT_CRU), 
                                   data=Ratio_NEP_GPP,
                                   b = 100, type = "lmg", rank = T, diff = F, rela = TRUE)
print(VarImp_NEP_GPP)

## Mean site
-----------------------------------------------------------------

#Compute transform function
Ratio_NEP_GPP_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPP_Mean_Site$Site_ID)){
  lm.Age<-   try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != i,], 
                       start = list(A= 0.11795, k= -0.03746, offset= -1.5), control = list(maxiter = 500)), silent=TRUE);  
  Age.pred = predict(object = lm.Age, newdata = Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == id,])
  Ratio_NEP_GPP_Mean_Site$f_Age[Ratio_NEP_GPP_Mean_Site$Site_ID == id] <- Age.pred
}

# Check for variable correlation
correlationMatrix <- cor(Ratio_NEP_GPP_Mean_Site[c("Stand_Age", "f_Age", "f_GPP", "f_Tair", "MAT_CRU", "SPI_CRU", 
                                                   "Clay_Silt", "CEC_Total_1km", "Soil_Quality")], use="complete")
print(correlationMatrix)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
print(highlyCorrelated)

# Stepwise regression
lm.Ratio_NEP_GPP_Mean_Site<-lm(values ~ (f_Age + I(GPP)^2 + Stand_Age + I(MAT_CRU) + SPI_CRU + NHx)^2, data=Ratio_NEP_GPP_Mean_Site)
step.Ratio_NEP_GPP_Mean_Site<- stepAIC(lm.Ratio_NEP_GPP_Mean_Site, direction = "backward")
summary(step.Ratio_NEP_GPP_Mean_Site)

# Check for residuals normality and homoskedasticity
res1=residuals(lm.Ratio_NEP_GPP_Mean_Site)
shapiro.test(res1)
plot(step.Ratio_NEP_GPP_Mean_Site) # ncreasing trend in the scale location plot shows heteroskedasticity

# Assess model performance
Ratio_NEP_GPP_Mean_Site$prediction <- NA
for(id in unique(Ratio_NEP_GPP_Mean_Site$Site_ID)){
  lm.Ratio_NEP_GPP_Mean_Site<- lm(values ~   f_Age + I(GPP) + Stand_Age + I(MAT_CRU) + 
                                    NHx + f_Age:I(GPP) + f_Age:SPI_CRU + f_Age:NHx + 
                                    I(GPP):I(MAT_CRU) + Stand_Age:I(MAT_CRU), 
                                  data= Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID != id,])
  Ratio_NEP_GPP_Mean_Site.pred = predict(object = lm.Ratio_NEP_GPP_Mean_Site, newdata =  Ratio_NEP_GPP_Mean_Site[Ratio_NEP_GPP_Mean_Site$Site_ID == id,])
  Ratio_NEP_GPP_Mean_Site$prediction[Ratio_NEP_GPP_Mean_Site$Site_ID == id] <- Ratio_NEP_GPP_Mean_Site.pred
}

R2_Ratio_NEP_GPP_Mean_Site<- cor(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values, use="complete")^2
RMSE_Ratio_NEP_GPP_Mean_Site <- (sum((Ratio_NEP_GPP_Mean_Site$prediction-Ratio_NEP_GPP_Mean_Site$values)^2)/length(Ratio_NEP_GPP_Mean_Site$values))^(1/2)
NSE_Ratio_NEP_GPP_Mean_Site<-NSE(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values, na.rm=TRUE)
Bias_Ratio_NEP_GPP_Mean_Site<-pbias(Ratio_NEP_GPP_Mean_Site$prediction, Ratio_NEP_GPP_Mean_Site$values) 

# Compute Importance variable
# t-statistic
VarImp_NEP_GPP_Mean_Site <-varImp(lm(values~  f_Age + I(GPP) + Stand_Age + I(MAT_CRU) + 
                                       NHx + f_Age:I(GPP) + f_Age:SPI_CRU + f_Age:NHx + 
                                       I(GPP):I(MAT_CRU) + I(GPP):NHx + Stand_Age:I(MAT_CRU),
                                     data=Ratio_NEP_GPP_Mean_Site), scale=T)

#LMG method
VarImp_NEP_GPP_Mean_Site<- boot.relimp(values~  f_Age + I(GPP) + Stand_Age + I(MAT_CRU) + 
                                         NHx + f_Age:I(GPP) + f_Age:SPI_CRU + f_Age:NHx + 
                                         I(GPP):I(MAT_CRU) + I(GPP):NHx + Stand_Age:I(MAT_CRU), 
                              data=Ratio_NEP_GPP_Mean_Site,
                              b = 100, type = "lmg", rank = T, diff = F, rela = TRUE)
print(VarImp_NEP_GPP_Mean_Site)

# 2. Plot observed vs. actual for the different flux

# 2.1 All years per site

# Prepare dataset for plotting
pred_NEP<- NEP[c("Type_Flux", "values", "prediction", "Stand_Age")]
pred_Ratio_NEP_GPP<- Ratio_NEP_GPP[c("Type_Flux", "values", "prediction", "Stand_Age")]
pred_All<- rbind(pred_NEP, pred_Ratio_NEP_GPP)
levels(pred_All$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "GPP_ER", "NEP_GPP")
df_NEP<-pred_All[pred_All$Type_Flux %in% c("NEP"),]
df_NEP_GPP<-pred_All[pred_All$Type_Flux %in% c("Ratio NEP-GPP"),]

#Plot prediction vs. observation

# NEP
gg1<- ggplot(NEP, aes(x=prediction, y=values, colour=Stand_Age))+
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
  annotate("text", label = "R-squared = 0.68", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 142.6 gC.m-2.y-1", x = -370, y = 550, size =4)

# Ratio NEP-GPP
gg2<- ggplot(Ratio_NEP_GPP, aes(x=prediction, y=values, colour=Stand_Age))+
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
  annotate("text", label = "R-squared = 0.61", x = -1.1, y = 0.8, size =4) +
  annotate("text", label = "RMSE = 0.16", x = -1.20, y = 0.55, size =4)

# Create an arrange plot object
pdf("Latex/Figures/Pred_Flux_All_Site.eps", width = 5, height = 12) # Open a new pdf file
source("Function/Legend_Grid_Arrange.R")
grid_arrange_shared_legend(gg1, gg2, nrow=3) # Write the grid.arrange in the file

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

#3. Plot prediction vs. observation

# NEP
gg1<- ggplot(NEP_Mean_Site, aes(x=prediction, y=values, colour=Stand_Age))+
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
  annotate("text", label = "R-squared = 0.71", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 139.4 gC.m-2.y-1", x = -370, y = 550, size =4)

# Ratio NEP-GPP
gg2<- ggplot(Ratio_NEP_GPP_Mean_Site, aes(x=prediction, y=values, colour=Stand_Age))+
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
  annotate("text", label = "RMSE = 0.16", x = -1.20, y = 0.55, size =4)

#Plot all plots together
pdf("Latex/Figures/Pred_Flux_Mean_Site.eps", width = 5, height = 12) # Open a new pdf file
source("Function/Legend_Grid_Arrange.R")
grid_arrange_shared_legend(gg1, gg2, nrow=3) # Write the grid.arrange in the file

# 4. Perform boked plots

# NEP
NEP_Bokeh <- figure(xlab = "Prediction",
            ylab = "Observation", legend_location = "top_left", xlim=c(-700, 800),
              ylim=c(-700, 850)) %>%
  ly_points(prediction, values, data = NEP,
            color = Stand_Age,
            hover = list(Site_ID, GPP, Stand_Age))%>%
  ly_abline(0,1)

#NEP - Average site
NEP_Av_Bokeh <- figure(xlab = "Prediction",
                    ylab = "Observation", legend_location = "top_left", xlim=c(-700, 800),
                    ylim=c(-700, 850)) %>%
  ly_points(prediction, values, data = NEP_Mean_Site,
            color = Stand_Age,
            hover = list(Site_ID, GPP, Stand_Age))%>%
  ly_abline(0, 1)

#CUEe
CUE_Bokeh <- figure(xlab = "Prediction",
                    ylab = "Observation", legend_location = "top_left", xlim=c(-1.5, 1),
                      ylim=c(-1.5, 1)) %>%
  ly_points(prediction, values, data = Ratio_NEP_GPP,
            color = Stand_Age, 
            hover = list(Site_ID, GPP, Stand_Age))%>%
  ly_abline(0, 1)

#CUEe - Average sute
CUE_Av_Bokeh <- figure(xlab = "Prediction",
                    ylab = "Observation", legend_location = "top_left", xlim=c(-1.5, 1),
                      ylim=c(-1.5, 1)) %>%
  ly_points(prediction, values, data = Ratio_NEP_GPP_Mean_Site,
            color = Stand_Age,
            hover = list(Site_ID, GPP, Stand_Age))%>%
  ly_abline(0, 1)

# Link plot together
Pred_Obs_NEP<- grid_plot(list(NEP_Bokeh, NEP_Av_Bokeh), same_axes = F, link_data = F, nrow=1)
rbokeh2html(Pred_Obs_NEP, file="/media/simonbesnard/External_SB/SimonBesnard/PhD_MPI/Presentation/Figures/Pred_Obs_NEP.html")
Pred_Obs_CUE<- grid_plot(list(CUE_Bokeh, CUE_Av_Bokeh), same_axes = F, link_data = F, nrow=1)
rbokeh2html(Pred_Obs_CUE, file="/media/simonbesnard/External_SB/SimonBesnard/PhD_MPI/Presentation/Figures/Pred_Obs_CUE.html")
