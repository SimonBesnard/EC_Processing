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
library(hexbin)
library(RColorBrewer)
library(plot3D)
library(ggRandomForests)
library(randomForestSRC)
library (reshape)

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

#Subset the data
GPP_High<-data.frame(cbind(GPP_High$Site_ID, GPP_High$values, GPP_High$Annual_Preci, GPP_High$Stand_Age,
                           GPP_High$Tair, GPP_High$Tsoil,
                           GPP_High$Rg))

colnames(GPP_High)<-c("GPP", "Precipitation", "Stand_Age", "Tair", "Tsoil", "Rg")

# plot panels for each covariate colored by the logical chas variable.
dta <- melt(GPP_High, id.vars=c("GPP"))

ggplot(dta, aes(x=GPP, y=value))+
  geom_point(alpha=.4)+
  geom_rug(data=dta %>% filter(is.na(value)))+
  # labs(y="", x=st.labs["GPP"]) +
  scale_color_brewer(palette="Set2")+
  facet_wrap(~variable, scales="free_y", ncol=2)+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw(base_size = 12, base_family = "Helvetica") 

## Create a basic Random Forest

# Create Training & Test Sets
GPP_High<-na.omit(GPP_High)
rownames(GPP_High)<-1:nrow(GPP_High)
rows <- sample(x=1:nrow(GPP_High),size=0.7 *  nrow(GPP_High))
train <- GPP_High[rows,]
test <-GPP_High[!rownames(GPP_High) %in% rows,]

## Create a Random forest with leave one out ID CV
nobs <- nrow(na.omit(train))
crossclass <- sample(5,nobs,TRUE)
crossPredict <- rep(NA,nobs)
for (i in 1:5) {
  indtrain<-which(crossclass!=i)
  indvalidate<-setdiff(1:nobs,indtrain)
  rf_GPP_CV<-randomForest(GPP~., data=train[indtrain,], ntree = 1000, importance = T, na.action = "na.omit")
  crossPredict[indvalidate]<-predict(rf_GPP_CV, train[indvalidate,])
}

#Compute correlation
cor_rf_GPP_CV<-(cor(crossPredict, train$GPP, use='pairwise'))^2

# Plot variable importance importance .
var_importance <- data_frame(variable=setdiff(colnames(train), "GPP"),
                             importance=as.vector(importance(rf_GPP)))

var_importance <- arrange(var_importance, desc(importance))
var_importance$variable <- factor(var_importance$variable, levels=var_importance$variable)

p <- ggplot(var_importance, aes(x=variable, weight=importance, fill=variable))
p <- p + geom_bar()
p <- p + xlab("Demographic Attribute") + ylab("Variable Importance (Mean Decrease in Gini Index)")
p <- p + scale_fill_discrete(name="Variable Name")
p + theme(axis.text.x=element_blank(),
          axis.text.y=element_text(size=12),
          axis.title=element_text(size=16),
          plot.title=element_text(size=18),
          legend.title=element_text(size=16),
          legend.text=element_text(size=12))+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw(base_size = 12, base_family = "Helvetica") 

## Model prediction

# Estimate MSE and RMSE
pred<-predict(object=rf_GPP_CV, newdata=test)
actual<-test$GPP
result<-data.frame(actual=actual,predicted=pred)
paste('Function Call: ', rf_GPP_CV$call)
paste('Mean Squared error: ',mean(rf_GPP_CV$mse))
paste('Root Mean Squared error: ',mean(sqrt(rf_GPP_CV$mse)))

#Plot observation vs. prediction
ggplot(result)+
  geom_point(aes(x=actual,y=predicted,color=predicted-actual),alpha=0.7)+
  ggtitle('Plotting Error')+
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_bw(base_size = 12, base_family = "Helvetica") 