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

#Subset data set
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]

NEP<- Flux_High[Flux_High$Type_Flux %in% c("NEP"),]
GPP<- Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco<- Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
Ratio_NEP_GPP<- Flux_High[Flux_High$Type_Flux %in% c("NEP_GPP"),]
Ratio_GPP_Reco<- Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]

#Compute GPPmax based on the lieth model
params <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$Tair))))
GPP$GPPp<-with(as.list(params), GPP1000*((1-exp(-k*GPP$Annual_Preci))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))
GPP$GPPmax[GPP$GPPmax == 0] <- NA
GPP$Scalar<- GPP$values/GPP$GPPmax

# Create dataset mean per site
NEP_Mean_Site<-ddply(NEP, .(Site_ID, Type_Flux),
                          summarise,
                          Stand_Age= mean(Stand_Age, na.rm=T),
                          values=mean(values, na.rm=T),
                          Precipitation=mean(Annual_Preci, na.rm=T),
                          Tair=mean(Tair, na.rm=T))

GPP_Mean_Site<-ddply(GPP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Precipitation=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))

Reco_Mean_Site<-ddply(Reco, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Precipitation=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))

Ratio_NEP_GPP_Mean_Site<-ddply(Ratio_NEP_GPP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Precipitation=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))

Ratio_GPP_Reco_Mean_Site<-ddply(Ratio_GPP_Reco, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Precipitation=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))



# 1.2 Test the best fitting function to ecosystem response - Modelling efficiency (Nash-Sutcliffe Efficiency test)

# Load function 
source("Function/NSE.R")

# Create a list of dataframe
df.list<-list(NEP, NEP_Mean_Site, GPP, GPP_Mean_Site, Reco, Reco_Mean_Site, Ratio_GPP_Reco, Ratio_GPP_Reco_Mean_Site, Ratio_NEP_GPP, Ratio_NEP_GPP_Mean_Site)
stat(Ratio_NEP_GPP_Mean_Site)

# Create a dataframe with the output results
Out_NSE<- lapply(df.list, stat)
names(Out_NSE_Harvest)<-c("NEP-Harvest", "Mean NEP-Harvest", "GPP-Harvest", "Mean GPP-Harvest", "Reco-Harvest", "Mean Reco-Harvest",
                  "Ratio-Harvest", "Mean Ratio-Harvest")
names(Out_NSE)<-c("GPP-Fire", "Mean GPP-Fire", "Reco-Fire", "Mean Reco-Fire","Ratio-Fire", "Mean Ratio-Fire")
df1<-melt(Out_NSE)
colnames(df1)[2]<-"Flux"
df1$Function<-c("Gamma", "Second Polymonial", "Third Polymonial", "Amiro", "Asymptote")

# Plot this output of the test
df1$Flux <- factor(df1$Flux, levels = c("NEP-Harvest", "Mean NEP-Harvest", "GPP-Harvest", "Mean GPP-Harvest", "Reco-Harvest", "Mean Reco-Harvest",
                                           "Ratio-Harvest", "Mean Ratio-Harvest","GPP-Fire", "Mean GPP-Fire", "Reco-Fire", "Mean Reco-Fire", 
                                           "Ratio-Fire", "Mean Ratio-Fire"))

gg1<-ggplot(df1, aes(x =Function, y = value)) +
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

# 2.1 NEP

# Compute the best fit function
Fun_NEP<-nls(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data = NEP, 
             start = list(A=0.02, B=-0.6, C= 50, D=200))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP$Stand_Age),max(NEP$Stand_Age),length=50)
pred1 <- approx(NEP$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg2<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = NEP, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
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

# 2.2 GPP

# Add artificial values to GPP
newrow = list(NA, NA, NA, NA, "GPP", numeric(1), NA, NA, NA, NA, NA, NA, NA, NA, numeric(1), NA, NA, NA, NA, NA, NA)

# Compute the best fit function
GPP<-rbind(newrow, GPP)
Fun_GPP<-nls(values~A*Stand_Age^2+B*Stand_Age, data = GPP, 
             start = list(A=-0.4, B=50))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_GPP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(GPP$Stand_Age),max(GPP$Stand_Age),length=50)
pred1 <- approx(GPP$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(GPP$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(GPP$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg3<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = GPP, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
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

# 2.3 Reco Harvest

# Compute the best fit function
Fun_Reco<-nls(values~A*Stand_Age^2+B*Stand_Age+C, data = Reco, 
             start = list(A=-0.4, B=50, C= 300))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Reco), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Reco$Stand_Age),max(Reco_Mean_Site$Stand_Age),length=50)
pred1 <- approx(Reco$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Reco$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Reco$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg4<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = Reco, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
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

# 2.4 Ratio GPP/Reco Harvest

# Compute the best fit function
Fun_Ratio_GPP_Reco_Mean_Site<-nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = Ratio_GPP_Reco_Mean_Site, 
                                 start = list(A = 1000, B = 0.170, k = -0.00295))
# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_GPP_Reco_Mean_Site), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_GPP_Reco_Mean_Site$Stand_Age),max(Ratio_GPP_Reco_Mean_Site$Stand_Age),length=50)
pred1 <- approx(Ratio_GPP_Reco_Mean_Site$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_GPP_Reco_Mean_Site$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_GPP_Reco_Mean_Site$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg5<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = Ratio_GPP_Reco_Mean_Site, aes(x = Stand_Age, y = values, size=Precipitation, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Ratio GPP-Reco")+ 
  ylim(0,2)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  geom_hline(yintercept=1, colour='grey', lty="dashed", size=0.8)

# 2.5 Ratio NEP/GPP Harvest

# Compute the best fit function
Fun_Ratio_NEP_GPP_Mean_Site<-nls(values ~ SSweibull(Stand_Age, Asym, Drop, lrc, pwr), data=Ratio_NEP_GPP)

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPP_Mean_Site), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPP$Stand_Age),max(Ratio_NEP_GPP_Mean_Site$Stand_Age),length=50)
pred1 <- approx(Ratio_NEP_GPP$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPP$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPP$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg6<-ggplot(predVals, aes(x, lower, upper)) +
  geom_line(aes(y = fit), colour="black", size=0.8)+
  geom_line(mapping = aes(y = upper), lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), lty = "dashed", size=0.8)+
  geom_point(data = Ratio_NEP_GPP, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Ratio NEP-GPP")+ 
  # ylim(0,3000)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  geom_hline(yintercept=0, colour='grey', lty="dashed", size=0.8)
  
# 2.9.Plot all carbon fluxes plots together

# Create an arrange plot object
pdf("Latex/Figures/Annual_Flux.eps", width = 15, height = 12) # Open a new pdf file
grid.arrange(gg2, gg3, gg4, gg5, gg6, nrow=3) # Write the grid.arrange in the file
dev.off() # Close the file
