## Script to compute statistical analysis annual carbon flux
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library (ggplot2)
library(scales)
library(lubridate)
library (dplyr)
library (plyr)
library(tidyr)
library(manipulate)
library(gridExtra)
library (bootstrap)
library(randomForest)
library (reshape)
library(GGally)
library(e1071) 
library(ROCR)


#1 Explain variabilty of the fluxes

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")

#1.1. Subset dataset
Flux_High<-dfAll_Sites[dfAll_Sites$Int_Replacement %in% c("High"),]
GPP_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco_High<-Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
Ratio_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]
NEP_High<-Flux_High[Flux_High$Type_Flux %in% c("NEP"),]

# Harvest
GPP_High_Harvest<-GPP_High[GPP_High$Disturbance %in% c("Harvest"),]
Reco_High_Harvest<-Reco_High[Reco_High$Disturbance %in% c("Harvest"),]
NEP_High_Harvest<-NEP_High[NEP_High$Disturbance %in% c("Harvest"),]
Ratio_High_Harvest<-Ratio_High[Ratio_High$Disturbance %in% c("Harvest"),]

#Fire
GPP_High_Fire<-GPP_High[GPP_High$Disturbance %in% c("Wildfire"),]
Reco_High_Fire<-Reco_High[Reco_High$Disturbance %in% c("Wildfire"),]
Ratio_High_Fire<-Ratio_High[Ratio_High$Disturbance %in% c("Wildfire"),]
NEP_High_Fire<-NEP_High[NEP_High$Disturbance %in% c("Wildfire"),]

# 1. 2. Compute Random Forest on the annual flux

# 1.2.1 Compute R2-squared for linear correlaltion

#Subset the data
GPP_High<-data.frame(cbind(GPP_High$Site_ID, GPP_High$values, GPP_High$Annual_Preci, GPP_High$Stand_Age,
                           GPP_High$Tair, GPP_High$Tsoil,
                           GPP_High$Rg, GPP_High$Rn, GPP_High$LE, GPP_High$ET))

colnames(GPP_High)<-c("Site_ID", "GPP", "Precipitation", "Stand_Age", "Tair", "Tsoil", "Rg", "Rn", "LE", "ET")

# plot panels for each covariate colored by the logical chas variable.
gg1<-ggpairs(with(GPP_High, data.frame(GPP, Precipitation, Tair, Tsoil, Rg, Rn, LE, Stand_Age)))+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust=-1))+
  theme_bw(base_size = 12, base_family = "Helvetica")
print(gg1)

# 1.2.2. Compute multivariate analysis with Random Forest

# Create Training dataset
GPP_High<-na.omit(GPP_High)
tvec<-unique(GPP_High$Site_ID)
nruns <- length(tvec)
crossclass<-sample(nruns,length(tvec),TRUE)
nobs<-nrow(GPP_High)
crossPredict<-rep(NA,nobs)

#Run a RandomForest with leave one out ID CV
for (i in 1:nruns) {
  indtrain<-which(GPP_High$Site_ID %in% tvec[crossclass!=i])
  indvalidate<-setdiff(1:nobs,indtrain)
  rf<-randomForest(formula = GPP ~ Precipitation + Tair + Tsoil + LE + Rn + Stand_Age, data=GPP_High, subset=indtrain, ntree=5000, importance=T)
  crossPredict[indvalidate]<-predict(rf,GPP_High[indvalidate,])
}

#Compute correlation
cor_rf_GPP_CV<-(cor(crossPredict, GPP_High$GPP, use='pairwise'))^2

# Plot importance of variable from random forest
plot(rf, log="y")
varImpPlot(rf)

## Model prediction

# Estimate MSE and RMSE
pred<-crossPredict
actual<-GPP_High$GPP
result<-data.frame(actual=actual,predicted=pred)
paste('Function Call: ', rf$call)
paste('Mean Squared error: ',mean(rf$mse))
paste('Root Mean Squared error: ',mean(sqrt(rf$mse)))

#Plot observation vs. prediction
ggplot(result)+
  geom_point(aes(x=actual,y=predicted,color=predicted-actual),alpha=0.7)+
  ggtitle('Plotting Error')+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw(base_size = 12, base_family = "Helvetica")

