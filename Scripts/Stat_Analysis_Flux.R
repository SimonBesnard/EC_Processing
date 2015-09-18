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

#1.1. Anova analysis

#Subset dataset
Flux_High<-dfAll_Sites[dfAll_Sites$Int_Replacement %in% c("High"),]
GPP_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco_High<-Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
NEE_High<-Flux_High[Flux_High$Type_Flux %in% c("NEE"),]
GPP_High_Harvest<-GPP_High[GPP_High$Disturbance %in% c("Harvest"),]
Reco_High_Harvest<-Reco_High[Reco_High$Disturbance %in% c("Harvest"),]
GPP_High_Fire<-GPP_High[GPP_High$Disturbance %in% c("Wildfire"),]
Reco_High_Fire<-Reco_High[Reco_High$Disturbance %in% c("Wildfire"),]
Ratio_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]
Ratio_High_Harvest<-Ratio_High[Ratio_High$Disturbance %in% c("Harvest"),]
Ratio_High_Fire<-Ratio_High[Ratio_High$Disturbance %in% c("Wildfire"),]

#Compute linear model
fit = lm(values ~ Annual_Preci + Stand_Age +Tair + Tsoil + Rg, data=GPP_High)

m5 = fit
m4 = update(m5, ~ . - Tair)
m3 = update(m4, ~ . - Tsoil)
m2=  update(m3, ~ . - Rg)
m1=  update(m2, ~ . - Stand_Age)
m0=  update(m1, ~ . - Annual_Preci)

af<-anova(m0,m1,m2,m3,m4,m5)
summary(af)
afss <- af$"Sum of Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))

#Compute performance of the model
af<-step(fit, direction = c("forward"), steps = 2000)
summary(af)

# 1.2 Principal Components Analysis

#Subset the data
Reco_High_PCA<-data.frame(cbind(Reco_High$values, Reco_High$Annual_Preci, Reco_High$Stand_Age,
                                Reco_High$Tair, Reco_High$Tsoil,
                                Reco_High$Rg))

colnames(Reco_High_PCA)<-c("Reco", "Precipitation", "Stand_Age", "Tair", "Tsoil", "Rg")

#Plot the pairwise scatterplots
gg7<-ggpairs(Reco_High_PCA,
             columns = c(1,2,3,4,5,6),
             legends=T,
             lower = list(continuous = "points"),
             diag = list(continuous = "density"),
             axisLabels = "none",
             title = "Correlation matrix - Reco")+ 
  theme_bw(base_size = 14, base_family = "Helvetica")
print(gg7)

# Perform the principal component analysis
arc.pca1<-prcomp(na.omit(Reco_High_PCA), scores=T, cor=T)
summary(arc.pca1, loadings=T)

# create data frame with scores
scores = as.data.frame(arc.pca1$x)

# plot of observations
gg8<-ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "tomato", alpha = 0.8, size = 4)
print(gg8)

# function to create a circle of correlations
circle <- function(center = c(0, 0), npoints = 100) {
  r = 1
  tt = seq(0, 2 * pi, length = npoints)
  xx = center[1] + r * cos(tt)
  yy = center[1] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}
corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(na.omit(Reco_High_PCA), arc.pca1$x))

# data frame with arrows coordinates
arrows = data.frame(x1 = c(0, 0, 0, 0, 0, 0), y1 = c(0, 0, 0, 0, 0, 0), x2 = correlations$PC1, 
                    y2 = correlations$PC2)

# geom_path will do open circles
gg9<-ggplot() + geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
  geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "red") + 
  geom_text(data = correlations, aes(x = PC1, y = PC2, label = rownames(correlations))) + 
  geom_hline(yintercept = 0, colour = "gray65") + geom_vline(xintercept = 0, colour = "gray65")+
  xlim(-1.1, 1.1) + ylim(-1.1, 1.1) + labs(x = "pc1 aixs", y = "pc2 axis")+
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggtitle("PCA-Reco")
print(gg9)

# 2. Compute Random Forest on the annual flux

#Subset the data
GPP_High<-data.frame(cbind(GPP_High$values, GPP_High$Annual_Preci, GPP_High$Stand_Age,
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
set.seed(42)
GPP_High<-na.omit(GPP_High)
rownames(GPP_High)<-1:nrow(GPP_High)
rows <- sample(x=1:nrow(GPP_High),size=0.7 *  nrow(GPP_High))
train <- GPP_High[rows,]
test <-GPP_High[!rownames(GPP_High) %in% rows,]

#Train a Random Foret 
set.seed(42)
rf_GPP<-randomForest(GPP~., data=train, ntree = 1000)

# Compute correlation
cor_rf_GPP<-(cor(predict(rf_GPP,train), train$GPP, use="pairwise"))^2

## Create a final Random forest
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