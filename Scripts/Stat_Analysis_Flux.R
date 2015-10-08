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
library(caret)
library (verification)

#1 Explain variabilty of the fluxes

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")

#1.1. Subset dataset
Flux_High<-dfAll_Sites[dfAll_Sites$Int_Replacement %in% c("High"),]

#Subset type of disturbance
GPP_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
NEP_High<-Flux_High[Flux_High$Type_Flux %in% c("NEP"),]
Reco_High<-Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
Ratio_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]

# 1. 2. Compute Random Forest on the annual flux

# 1.2.1 Compute R2-squared for linear correlaltion

# Transform inf values into NA values
is.na(NEP_High) <- sapply(NEP_High, is.infinite)
var<-c("Annual_Preci", "Tair", "Tsoil", "Rg", "Rn", "LE", "ET")
NEP_High[,var] <- sapply(NEP_High[,var],function(x) ifelse(x==0,NA,x))

# plot panels for each covariate colored by the logical chas variable.
gg1<-ggpairs(with(NEP_High, data.frame(NEP, Precipitation, Tair, Tsoil, Rn, LE, ET, Stand_Age)))+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust=-1))+
  theme_bw(base_size = 12, base_family = "Helvetica")
print(gg1)

# 1.2.2. Compute multivariate analysis with Random Forest

# Create training and test dataset
NEP_High<-na.omit(NEP_High)
subs <- unique(NEP_High$Site_ID)
train<- vector(mode = "list", length = length(subs))
test<- vector(mode = "list", length = length(subs))

# Run a RandomForest with leave one out ID CV and Predict the outcome of the testing data
predicted <- NULL
for(i in seq_along(subs)){
  train[[i]] <- subset(NEP_High[NEP_High$Site_ID != subs[i],])
  rf <- randomForest(values~ Annual_Preci + Tair + Tsoil + Rg + Stand_Age + Disturbance, data=train[[i]],
                     ntree=5000, keep.forest=T, importance=T)
  test[[i]] <- subset(NEP_High[NEP_High$Site_ID == subs[i],])
  predicted[[i]] <- predict(rf, newdata=test[[i]])
  }

# Compute R-squared estimates
predi<-melt(predicted)
actual <- NEP_High$values
rsq <- 1-sum((actual-predi$value)^2)/sum((actual-mean(actual))^2)
print(rsq)

# Plot importance of variable from random forest
imp<-varImp(rf, useModel = T)

## Model prediction

# Estimate MSE and RMSE
pred<-predi
actual<-NEP_High$values
result<-data.frame(Observation=actual,Prediction=pred)

#Plot observation vs. prediction
gg2<-ggplot(result)+
  geom_point(aes(x=Observation,y=Prediction.value))+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw(base_size = 12, base_family = "Helvetica")
print(gg2)
ggsave("Latex/Figures/Model_pred_NEP.eps", height = 12, width = 15)
