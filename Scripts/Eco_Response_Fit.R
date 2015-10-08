## Script to fit ecosystem carbon flux response
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library (ggplot2)
library(scales)
library (dplyr)
library (plyr)
library(gridExtra)
library (bootstrap)
library (nlstools)
library(nls2)
library (boot)
library(tidyr)
library (reshape)
library (hydroGOF)
library (minpack.lm)
library (cvTools)

#.1.Function fit choice for ecosystem response

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")

# 1.1. Subset dataframe for fitting process 
Flux_High<-dfAll_Sites[dfAll_Sites$Int_Replacement %in% c("High"),]

#Subset type of disturbance
# Harvest
Harvest_High<-Flux_High[Flux_High$Disturbance %in% c("Harvest"),]

#Fire
Fire_High<-Flux_High[Flux_High$Disturbance %in% c("Wildfire"),]

# Add mature forest to the sub-datasets
Flux_Mature<-dfAll_Sites[dfAll_Sites$Int_Replacement %in% c("None"),]
Flux_Mature<-Flux_Mature[Flux_Mature$Stand_Age >100,]
Harvest_High<-rbind(Harvest_High,Flux_Mature)
Fire_High<-rbind(Fire_High,Flux_Mature)

#Subset type of flux for each disturbance

#NEP
NEP_High_Harvest<-Harvest_High[Harvest_High$Type_Flux %in% c("NEP"),]
NEP_High_Fire<-Fire_High[Fire_High$Type_Flux %in% c("NEP"),]
NEP_High_Mean_Harvest<-ddply(NEP_High_Harvest, .(Site_ID, Type_Flux),
                          summarise,
                          Stand_Age= median(Stand_Age, na.rm=T),
                          values=mean(values, na.rm=T),
                          Precipitation=mean(Annual_Preci, na.rm=T),
                          Tair=mean(Tair, na.rm=T))

NEP_High_Mean_Fire<-ddply(NEP_High_Fire, .(Site_ID, Type_Flux),
                          summarise,
                          Stand_Age= mean(Stand_Age, na.rm=T),
                          values=median(values, na.rm=T),
                          Precipitation=mean(Annual_Preci, na.rm=T),
                          Tair=mean(Tair, na.rm=T))

#GPP
GPP_High_Harvest<-Harvest_High[Harvest_High$Type_Flux %in% c("GPP"),]
GPP_High_Fire<-Fire_High[Fire_High$Type_Flux %in% c("GPP"),]
GPP_High_Mean_Harvest<-ddply(GPP_High_Harvest, .(Site_ID, Type_Flux),
                             summarise,
                             Stand_Age= median(Stand_Age, na.rm=T),
                             values=mean(values, na.rm=T),
                             Precipitation=mean(Annual_Preci, na.rm=T),
                             Tair=mean(Tair, na.rm=T))

GPP_High_Mean_Fire<-ddply(GPP_High_Fire, .(Site_ID, Type_Flux),
                          summarise,
                          Stand_Age= mean(Stand_Age, na.rm=T),
                          values=median(values, na.rm=T),
                          Precipitation=mean(Annual_Preci, na.rm=T),
                          Tair=mean(Tair, na.rm=T))

#Reco
Reco_High_Harvest<-Harvest_High[Harvest_High$Type_Flux %in% c("Respiration"),]
Reco_High_Fire<-Fire_High[Fire_High$Type_Flux %in% c("Respiration"),]
Reco_High_Mean_Harvest<-ddply(Reco_High_Harvest, .(Site_ID, Type_Flux),
                             summarise,
                             Stand_Age= median(Stand_Age, na.rm=T),
                             values=mean(values, na.rm=T),
                             Precipitation=mean(Annual_Preci, na.rm=T),
                             Tair=mean(Tair, na.rm=T))

Reco_High_Mean_Fire<-ddply(Reco_High_Fire, .(Site_ID, Type_Flux),
                          summarise,
                          Stand_Age= mean(Stand_Age, na.rm=T),
                          values=median(values, na.rm=T),
                          Precipitation=mean(Annual_Preci, na.rm=T),
                          Tair=mean(Tair, na.rm=T))


#Ratio GPP-Reco
Ratio_High_Harvest<-Harvest_High[Harvest_High$Type_Flux %in% c("GPP_ER"),]
Ratio_High_Fire<-Fire_High[Fire_High$Type_Flux %in% c("GPP_ER"),]
Ratio_High_Mean_Harvest<-ddply(Ratio_High_Harvest, .(Site_ID, Type_Flux),
                             summarise,
                             Stand_Age= median(Stand_Age, na.rm=T),
                             values=mean(values, na.rm=T),
                             Precipitation=mean(Annual_Preci, na.rm=T),
                             Tair=mean(Tair, na.rm=T))

Ratio_High_Mean_Fire<-ddply(Ratio_High_Fire, .(Site_ID, Type_Flux),
                          summarise,
                          Stand_Age= median(Stand_Age, na.rm=T),
                          values=mean(values, na.rm=T),
                          Precipitation=mean(Annual_Preci, na.rm=T),
                          Tair=mean(Tair, na.rm=T))

# 1.2 Test the best fitting function to ecosystem response - Modelling efficiency (Nash-Sutcliffe Efficiency test)

# Load function 
source("Function/NSE1.R")
source("Function/NSE2.R")

# Create a list of dataframe
df.list1<-list(NEP_High_Harvest, NEP_High_Mean_Harvest, GPP_High_Harvest, GPP_High_Mean_Harvest, 
               Reco_High_Harvest, Reco_High_Mean_Harvest, Ratio_High_Harvest, Ratio_High_Mean_Harvest)
df.list2<-list(GPP_High_Mean_Fire, GPP_High_Fire, 
               Reco_High_Fire, Reco_High_Mean_Fire, Ratio_High_Fire, Ratio_High_Mean_Fire)

# Create a dataframe with the output results
Out_NSE_Harvest<- lapply(df.list1, stat1)
Out_NSE_Fire<-lapply(df.list2, stat2)
names(Out_NSE_Harvest)<-c("NEP-Harvest", "Mean NEP-Harvest", "GPP-Harvest", "Mean GPP-Harvest", "Reco-Harvest", "Mean Reco-Harvest",
                  "Ratio-Harvest", "Mean Ratio-Harvest")
names(Out_NSE_Fire)<-c("GPP-Fire", "Mean GPP-Fire", "Reco-Fire", "Mean Reco-Fire","Ratio-Fire", "Mean Ratio-Fire")
df1<-melt(Out_NSE_Harvest)
colnames(df1)[2]<-"Flux"
df1$Function<-c("Gamma", "Second Polymonial", "Third Polymonial", "Amiro", "Asymptote")
df2<-melt(Out_NSE_Fire)
colnames(df2)[2]<-"Flux"
df2$Function<-c("Gamma", "Second Polymonial", "Third Polymonial", "Amiro")
df.stat<-rbind(df1, df2)

# Plot this output of the test
df.stat$Flux <- factor(df.stat$Flux, levels = c("NEP-Harvest", "Mean NEP-Harvest", "GPP-Harvest", "Mean GPP-Harvest", "Reco-Harvest", "Mean Reco-Harvest",
                                           "Ratio-Harvest", "Mean Ratio-Harvest","GPP-Fire", "Mean GPP-Fire", "Reco-Fire", "Mean Reco-Fire", 
                                           "Ratio-Fire", "Mean Ratio-Fire"))

gg2<-ggplot(df.stat, aes(x =Function, y = value)) +
  geom_point() +
  facet_wrap(~ Flux, scales = "free", nrow=2)+
  xlab("") + ylab("Nash-Sutcliffe Efficiency")+
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(gg2)
ggsave("Latex/Figures/MEF_results.eps", height = 12, width = 15)

# 2. Plot ecosystem response with the best fit function
source("Function/CI_Est.R")

# 2.1 NEP Harvest

# Compute the best fit function
Fun_NEP_Harvest<-nls(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data=NEP_High_Harvest,
                     start = list(A=0.02, B=-0.6, C= 50, D=200))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP_Harvest), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP_High_Harvest$Stand_Age),max(NEP_High_Harvest$Stand_Age),length=50)
pred1 <- approx(NEP_High_Harvest$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP_High_Harvest$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP_High_Harvest$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg3<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = NEP_High_Harvest, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  # ylim(0,3000)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 2.2 GPP Harvest

# Compute the best fit function
Fun_GPP_Harvest<-nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data=GPP_High_Harvest,
                     start = list(A = 1000, B = 0.170, k = -0.00295))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_GPP_Harvest), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(GPP_High_Harvest$Stand_Age),max(GPP_High_Harvest$Stand_Age),length=50)
pred1 <- approx(GPP_High_Harvest$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(GPP_High_Harvest$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(GPP_High_Harvest$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg4<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = GPP_High_Harvest, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(-1,3000)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 2.3 Reco Harvest

# Compute the best fit function
Fun_Reco_Harvest<-nls(values~A*Stand_Age^2+B*Stand_Age+C, data=Reco_High_Harvest,
        start = list(A=-0.4, B=50, C= 300))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Reco_Harvest), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Reco_High_Harvest$Stand_Age),max(Reco_High_Harvest$Stand_Age),length=50)
pred1 <- approx(Reco_High_Harvest$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Reco_High_Harvest$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Reco_High_Harvest$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg5<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = Reco_High_Harvest, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(0,3000)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 2.4 Ratio GPP/Reco Harvest

# Compute the best fit function
Fun_Ratio_Harvest<-nls(values~A*(1-exp(k*Stand_Age)), data=Ratio_High_Harvest,
                       start = list(A=1500, k= -0.224))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_Harvest), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_High_Harvest$Stand_Age),max(Ratio_High_Harvest$Stand_Age),length=50)
pred1 <- approx(Ratio_High_Harvest$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_High_Harvest$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_High_Harvest$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg6<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = Ratio_High_Harvest, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Ratio GPP-Reco")+ 
  ylim(-0.0005,2)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  geom_hline(yintercept=1, linetype="dashed", colour="grey", size=0.8)

# 2.5 NEP Fire

# Compute the best fit function
Fun_NEP_Fire<-nls(values~A*Stand_Age^2+B*Stand_Age+C, data=NEP_High_Fire, 
                  start = list(A=-0.4, B=50, C= 300))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP_Fire), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP_High_Fire$Stand_Age),max(NEP_High_Fire$Stand_Age),length=50)
pred1 <- approx(NEP_High_Fire$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP_High_Fire$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP_High_Fire$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg7<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = NEP_High_Fire, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  # ylim(0,3000)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 2.6 GPP Fire

# Compute the best fit function
Fun_GPP_Fire<-nls(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data=GPP_High_Fire, 
                  start = list(A=0.02, B=-0.6, C= 50, D=200))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_GPP_Fire), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(GPP_High_Fire$Stand_Age),max(GPP_High_Fire$Stand_Age),length=50)
pred1 <- approx(GPP_High_Fire$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(GPP_High_Fire$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(GPP_High_Fire$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg8<-ggplot(predVals, aes(x, lower, upper)) +  
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = GPP_High_Fire, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(0,3000)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 2.7 Reco Fire

# Compute the best fit function
Fun_Reco_Fire<-nls(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data=Reco_High_Fire, 
                   start = list(A=0.02, B=-0.6, C= 50, D=200))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Reco_Fire), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Reco_High_Fire$Stand_Age),max(Reco_High_Fire$Stand_Age),length=50)
pred1 <- approx(Reco_High_Fire$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Reco_High_Fire$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Reco_High_Fire$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg9<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = Reco_High_Fire, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(0,3000)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 2.8 Ratio GPP/Reco Fire

# Compute the best fit function
Fun_Ratio_Fire<-nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data=Ratio_High_Fire,
                    start = list(A = 1000, B = 0.170, k = -0.00295))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_Fire), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_High_Fire$Stand_Age),max(Ratio_High_Fire$Stand_Age),length=50)
pred1 <- approx(Ratio_High_Fire$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_High_Fire$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_High_Fire$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg10<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = Ratio_High_Fire, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Ratio GPP-Reco")+ 
  ylim(0,3)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  geom_hline(yintercept=1, linetype="dashed", colour="grey", size=0.8)

# 2.9.Plot all carbon fluxes plots together

#2.9.1 Harvest

# Create an arrange plot object
pdf("Latex/Figures/Flux_Harvest.eps", width = 15, height = 12) # Open a new pdf file
grid.arrange(gg3, gg4, gg5, gg6, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file


#2.9.2 Fire

# Create an arrange plot object
pdf("Latex/Figures/Flux_Fire.eps", width = 15, height = 12) # Open a new pdf file
grid.arrange(gg7, gg8, gg9, gg10, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file

# Resisdual
plot(residuals(Fun_Ratio_Harvest), Ratio_High_Harvest$Tair)
(cor(residuals(Fun_Ratio_Harvest), Ratio_High_Harvest$Tair))
