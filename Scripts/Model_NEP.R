## Script to model NEP vs. stand age
## Author: Simon Besnard
## 12.10.2015
###################################
## Load the necessary packages
library(gridExtra)
library (ggplot2)
library(scales)
library (dplyr)
library (plyr)
library(tidyr)
library (hydroGOF)
library(grid)

# 1. Model NEP with semi-empirical model and GPP product

# 1.1 Prepare dataset

# Import GPP product
dir <- file.path(path, 'GPP_Extraction/data')
list <- list.files(dir, pattern=glob2rx('*.txt'), full.names=TRUE)
GPP_Site <- list()
for(i in seq_along(list)) {
  GPP_Site[[i]] = read.table(list[i], header=T,  check.names=FALSE, sep=",")[-c(1:2), ]
}

# Compute monthly GPP sum per site
for (i in seq_along(GPP_Site)){
  names(GPP_Site[[i]])[1]<- "ID"
  GPP_Site[[i]]<- gather(GPP_Site[[i]], variable, value, -ID )
  GPP_Site[[i]]<- GPP_Site[[i]][c("ID", "variable", "value")]
  colnames(GPP_Site[[i]])<- c("ID", "year", "GPP_RS")
  GPP_Site[[i]]<- ddply(GPP_Site[[i]], .(ID, year),
                        summarize,
                        GPP_RS= sum(GPP_RS, na.rm=T)*8)
}

# Merge list of dataframe
GPP_Site<- do.call("rbind", GPP_Site)

# Group dataframe per year and per site
GPP_Site$ID <- as.character(GPP_Site$ID)
GPP_Site$ID <- sapply(GPP_Site$ID, function(id) as.numeric(strsplit(id,"[.]")[[1]][2]))
GPP_Site <- GPP_Site[with(GPP_Site, order(ID)),]

# Include Site ID into the dataframe
dir <- file.path(path, 'GPP_Extraction')
list <- list.files(dir, pattern=glob2rx('*.csv'), full.names=TRUE)
Site_ID<- read.csv(list, header=T,  check.names=FALSE, sep=",")
GPP_Site<- merge(GPP_Site, Site_ID, by.x="ID", by.y="ID", all.x=TRUE)
GPP_Site<- GPP_Site[c("Site_ID", "year", "GPP_RS")]

# Join flux data and GPP product data 
NEP<-readRDS("Output/NEP.rds") # Import NEP dataframe
GPP<-readRDS("Output/GPP.rds") # Import NEP dataframe
NEP$GPP<- GPP$values
params <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$Tair))))
GPP$GPPp<-with(as.list(params), GPP1000*((1-exp(-k*GPP$Annual_Preci))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))
NEP$GPPmax<- GPP$GPPmax
NEP<- merge(NEP, GPP_Site, by=c("Site_ID", "year"), all.x=TRUE) #Merge GPP product

# 1.2 Compute NEP models from GPP eddy covariance

# 1.2.1 Model 1 - From NEP

# Compute transformed variables 
f_Age_NEP<- function (x) {3.418801e+02*(exp(-5.334230e-03*x)) -1.042789e+03*(exp(-1.589315e-01*x))} # Age model
f_Tair_NEP<- function (x) {-1.073267*x^2 +37.232210*x + 15.090029} # Tair model
f_Photo_NEP<- function (x) {-1.142427e-04*x^2+5.591564e-01*x -2.833088e+02} # GPP model

NEP$f_Tair<- f_Tair_NEP(NEP$Tair)
NEP$f_Age<- f_Age_NEP(NEP$Stand_Age)
NEP$f_GPP<- f_Photo_NEP(NEP$GPP)

# Perform stepwise regression
lm.NEP<-lm(values ~ (Annual_Preci + f_Tair + f_Age + f_GPP + SPI_CRU + CEC_Total_1km)^2, data=NEP)
step.NEP<- step(lm.NEP, direction = "backward")
print(step.NEP)
summary(step.NEP)

# Assess model performance
NEP$pred_Model1 <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "f_Age", "f_GPP", "SPI_CRU", "CEC_Total_1km")]
  lm.NEP<- lm(values ~Annual_Preci + f_Tair + f_Age + f_GPP + 
                SPI_CRU + CEC_Total_1km + Annual_Preci:f_Age + Annual_Preci:SPI_CRU + 
                f_Tair:f_GPP + f_Tair:SPI_CRU + f_Age:SPI_CRU + f_GPP:CEC_Total_1km + 
                SPI_CRU:CEC_Total_1km, data=train.df)
  NEP.pred = predict(object = lm.NEP, newdata = test.df)
  NEP$pred_Model1[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$pred_Model1, NEP$values, use = "complete.obs")^2
RMSE_NEP <- rmse(NEP$pred_Model1, NEP$values)
NSE_NEP<-NSE(NEP$pred_Model1, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$pred_Model1, NEP$values)

# 1.2.2 Model 2 - From NEP to GPP ratio

# Compute transformed variables 
f_Age_Ratio_NEP_GPP<- function (x) {0.225119297*(exp(-0.005007755*x)) -1.530424388*(exp(-0.170444337*x))} # Age model
NEP$f_Age_GPP<- f_Age_Ratio_NEP_GPP(NEP$Stand_Age)*NEP$GPPmax

# Perform stepwise regression
lm.NEP<-lm(values ~ (Annual_Preci + f_Tair + f_Age_GPP + f_GPP + SPI_CRU + CEC_Total_1km)^2, data=NEP)
step.NEP<- step(lm.NEP, direction = "backward")
print(step.NEP)
summary(step.NEP)

# Assess model performance
NEP$pred_Model2 <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_Tair", "f_Age_GPP", "f_GPP", "SPI_CRU", "CEC_Total_1km")]
  lm.NEP<- lm(values ~Annual_Preci + f_Tair + f_Age_GPP + f_GPP + 
                SPI_CRU + CEC_Total_1km + Annual_Preci:f_Age_GPP + f_Tair:f_GPP + 
                f_Tair:SPI_CRU + f_Age_GPP:SPI_CRU + f_GPP:CEC_Total_1km + 
                SPI_CRU:CEC_Total_1km, data=train.df)
  NEP.pred = predict(object = lm.NEP, newdata = test.df)
  NEP$pred_Model2[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$pred_Model2, NEP$values, use = "complete.obs")^2
RMSE_NEP <- rmse(NEP$pred_Model2, NEP$values)
NSE_NEP<-NSE(NEP$pred_Model2, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$pred_Model2, NEP$values)

# 1.2.3 Model 2 - From GPP to ER ratio

# Compute transformed variables 
f_Age_Ratio_GPP_Reco<- function (x) {1.2134719*(1-exp(-0.2257333*x))} # Age model
NEP$f_Age_GPP<- NEP$GPP- (NEP$GPPmax/f_Age_Ratio_GPP_Reco(NEP$Stand_Age))

# Perform stepwise regression
lm.NEP<-lm(values ~ (Annual_Preci + f_Tair + f_Age_GPP + f_GPP + SPI_CRU + CEC_Total_1km)^2, data=NEP)
step.NEP<- step(lm.NEP, direction = "backward")
print(step.NEP)
summary(step.NEP)

# Assess model performance
NEP$pred_Model3 <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "f_Tair", "f_Age_GPP", "f_GPP", "SPI_CRU", "CEC_Total_1km")]
  lm.NEP<- lm(values ~ Annual_Preci + f_Tair + f_Age_GPP + f_GPP + 
                SPI_CRU + CEC_Total_1km + Annual_Preci:SPI_CRU + f_Tair:f_Age_GPP + 
                f_Tair:f_GPP + f_Age_GPP:f_GPP + f_Age_GPP:SPI_CRU + f_GPP:CEC_Total_1km + 
                SPI_CRU:CEC_Total_1km, data=train.df)
  NEP.pred = predict(object = lm.NEP, newdata = test.df)
  NEP$pred_Model3[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$pred_Model3, NEP$values, use = "complete.obs")^2
RMSE_NEP <- rmse(NEP$pred_Model3, NEP$values)
NSE_NEP<-NSE(NEP$pred_Model3, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$pred_Model3, NEP$values)

# 1.2.4 Plot cross validated prediction vs. observations

# Model 1
gg1<- ggplot(NEP, aes(x=pred_Model1, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-700, 800)+
  ylim(-700, 800)+
  annotate("text", label = "R-squared = 0.43", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 189.9 gC.m-2.y-1", x = -390, y = 550, size =4)
  

# Model 2
gg2<- ggplot(NEP, aes(x=pred_Model2, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-700, 800)+
  ylim(-700, 800)+ 
  annotate("text", label = "R-squared = 0.44", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 187.5 gC.m-2.y-1", x = -390, y = 550, size =4)


# Model 3 
gg3<- ggplot(NEP, aes(x=pred_Model3, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-700, 800)+
  ylim(-700, 800)+
  annotate("text", label = "R-squared = 0.41", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 191.4 gC.m-2.y-1", x = -390, y = 550, size =4)


#Plot all plots together
pdf("Latex/Figures/Pred_NEP_EC.eps", width = 5, height = 12) # Open a new pdf file
source("Function/Legend_Grid_Arrange.R")
grid_arrange_shared_legend(gg1, gg2, gg3, nrow=4) # Write the grid.arrange in the file

# 1.3 Compute NEP models from GPP RS product

# 1.2.1 Model 1 - From NEP

# Compute transformed variables 
f_Age_NEP<- function (x) {3.437699e+02*(exp(-5.336560e-03*x)) -1.042647e+03*(exp(-1.578798e-01*x))} # Age model
f_Tair_NEP<- function (x) {-1.08456*x^2 +37.42899*x + 15.74958} # Tair model
f_P_NEP<- function (x) {2.762248e+02*(1-exp(-1.954383e-03*x))} # Precipitation model
f_Photo_NEP<- function (x) {-1.161296e-04*x^2+5.706779e-01*x -2.928120e+02} # GPP model

NEP$f_P<- f_P_NEP(NEP$Annual_Preci)
NEP$f_Tair<- f_Tair_NEP(NEP$Tair)
NEP$f_Age<- f_Age_NEP(NEP$Stand_Age)
NEP$f_GPP<- f_Photo_NEP(NEP$GPP_RS)

# Perform stepwise regression
lm.NEP<-lm(values ~ (f_P + f_Tair + f_Age + f_GPP)^2, data=NEP)
step.NEP<- step(lm.NEP, direction = "both")
print(step.NEP)
summary(step.NEP)

# Assess model performance
NEP$pred_Model1 <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age", "f_GPP")]
  lm.NEP<- lm(values ~ f_P + f_Tair + f_Age + f_GPP + f_P:f_Age + 
                f_P:f_GPP + f_Tair:f_Age, data=train.df)
  step.NEP<- step(lm.NEP, direction = "backward")
  NEP.pred = predict(object = step.NEP, newdata = test.df)
  NEP$pred_Model1[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$pred_Model1, NEP$values, use = "complete.obs")^2
RMSE_NEP <- rmse(NEP$pred_Model1, NEP$values)
NSE_NEP<-NSE(NEP$pred_Model1, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$pred_Model1, NEP$values)

# 1.2.2 Model 2 - From NEP to GPP ratio

# Compute transformed variables 
f_Age_Ratio_NEP_GPP<- function (x) {0.227464785*(exp(-0.005010447*x)) -1.529965842*(exp(-0.169444307*x))} # Age model
NEP$f_Age_GPP<- f_Age_Ratio_NEP_GPP(NEP$Stand_Age)*NEP$GPPmax

# Perform stepwise regression
lm.NEP<-lm(values ~ (f_P + f_Tair + f_Age_GPP + f_GPP)^2, data=NEP)
step.NEP<- step(lm.NEP, direction = "both")
print(step.NEP)
summary(step.NEP)

# Assess model performance
NEP$pred_Model2 <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age_GPP", "f_GPP")]
  lm.NEP<- lm(values ~f_P + f_Tair + f_Age_GPP + f_GPP + f_P:f_GPP + 
                f_Tair:f_Age_GPP + f_Age_GPP:f_GPP, data=train.df)
  step.NEP<- step(lm.NEP, direction = "backward")
  NEP.pred = predict(object = step.NEP, newdata = test.df)
  NEP$pred_Model2[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$pred_Model2, NEP$values, use = "complete.obs")^2
RMSE_NEP <- rmse(NEP$pred_Model2, NEP$values)
NSE_NEP<-NSE(NEP$pred_Model2, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$pred_Model2, NEP$values)

# 1.2.3 Model 2 - From GPP to ER ratio

# Compute transformed variables 
f_Age_Ratio_GPP_Reco<- function (x) {1.2159729*(1-exp(-0.2245381*x))} # Age model
NEP$f_Age_GPP<- NEP$GPP_RS- (NEP$GPPmax/f_Age_Ratio_GPP_Reco(NEP$Stand_Age))

# Perform stepwise regression
lm.NEP<-lm(values ~ (f_P + f_Tair + f_Age_GPP + f_GPP)^2, data=NEP)
step.NEP<- step(lm.NEP, direction = "both")
print(step.NEP)
summary(step.NEP)

# Assess model performance
NEP$pred_Model3 <- NA
for(id in unique(NEP$Site_ID)){
  train.df <- NEP[NEP$Site_ID != id,]
  test.df <- NEP[NEP$Site_ID == id, c("Annual_Preci", "Tair", "Stand_Age", "f_P", "f_Tair", "f_Age_GPP", "f_GPP")]
  lm.NEP<- lm(values ~f_P + f_Tair + f_Age_GPP + f_GPP + f_P:f_Tair + 
                f_P:f_Age_GPP + f_P:f_GPP + f_Age_GPP:f_GPP, data=train.df)
  step.NEP<- step(lm.NEP, direction = "backward")
  NEP.pred = predict(object = step.NEP, newdata = test.df)
  NEP$pred_Model3[NEP$Site_ID == id] <- NEP.pred
}

R2_NEP<- cor(NEP$pred_Model3, NEP$values, use = "complete.obs")^2
RMSE_NEP <- rmse(NEP$pred_Model3, NEP$values)
NSE_NEP<-NSE(NEP$pred_Model3, NEP$values, na.rm=TRUE)
Bias_NEP<-pbias(NEP$pred_Model3, NEP$values)

# 1.2.4 Plot cross validated prediction vs. observations

# Model 1
gg1<- ggplot(NEP, aes(x=pred_Model1, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-700, 800)+
  ylim(-700, 800)+
  annotate("text", label = "R-squared = 0.40", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 200.3 gC.m-2.y-1",  x = -390, y = 550, size =4)


# Model 2
gg2<- ggplot(NEP, aes(x=pred_Model2, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-700, 800)+
  ylim(-700, 800)+ 
  annotate("text", label = "R-squared = 0.31", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 210.7 gC.m-2.y-1",  x = -390, y = 550, size =4)


# Model 3 
gg3<- ggplot(NEP, aes(x=pred_Model3, y=values, colour=Stand_Age))+
  geom_point(shape=3)+
  geom_abline(slope=1)+
  facet_wrap(~Type_Flux) +
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(hjust = 1),
        legend.position="none", 
        legend.box="horizontal",
        legend.key = element_blank(),
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  xlab("Predicted")+
  ylab("Observed")+
  scale_colour_gradient(name= "Stand age", low= "#67a9cf", high ="#ef8a62", space="Lab",
                        breaks=c(1,250),labels=c("Young", "Old"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  xlim(-700, 800)+
  ylim(-700, 800)+
  annotate("text", label = "R-squared = 0.15", x = -500, y = 700, size =4) +
  annotate("text", label = "RMSE = 245.9 gC.m-2.y-1",  x = -390, y = 550, size =4)


#Plot all plots together
pdf("Latex/Figures/Pred_NEP_RS.eps", width = 5, height = 12) # Open a new pdf file
source("Function/Legend_Grid_Arrange.R")
grid_arrange_shared_legend(gg1, gg2, gg3, nrow=4) # Write the grid.arrange in the file

