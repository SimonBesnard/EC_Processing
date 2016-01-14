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
library(grid)

#1.Function fit choice for ecosystem response

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
#Min model
params <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$Tair))))
GPP$GPPp<-with(as.list(params), GPP1000*((1-exp(-k*GPP$Annual_Preci))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))
#Multiplicative model
# params <- c(a=-0.119, b=1.315, c=-0.000664)
# GPP$GPPmat<-with(as.list(params), 1/(1+exp(a*GPP$Tair+b))) 
# GPP$GPPp<-with(as.list(params),(1-exp(c*GPP$Annual_Preci)))
# GPP<-transform(GPP, GPPmax = 3000*GPPmat*GPPp)
GPP$NEP_GPPmax<- NEP$values/GPP$GPPmax
Ratio_NEP_GPPmax<- GPP
Ratio_NEP_GPPmax<-Ratio_NEP_GPPmax[, !(colnames(Ratio_NEP_GPPmax) %in% c("values", "Type_Flux", "GPPmat", "GPPmax", "GPPp"))]
Ratio_NEP_GPPmax<-gather(Ratio_NEP_GPPmax, variable, values, -Annual_Preci, -year, 
       -Ecosystem, -Climate, -Disturbance,
       -Stand_Age, -Site_ID, -Stand_Replacement, -Int_Replacement,
       -Tair, -Rg, -Study, -Lat, -Long)
Ratio_NEP_GPPmax<- Ratio_NEP_GPPmax[c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "variable", "values", "Annual_Preci", 
                             "Tair","Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")]
colnames(Ratio_NEP_GPPmax)<-c("Site_ID", "year", "Stand_Replacement", "Int_Replacement", "Type_Flux", "values", "Annual_Preci", 
                         "Tair", "Rg", "Stand_Age", "Disturbance", "Climate", "Ecosystem", "Study", "Lat", "Long")

# Create dataset mean per site
NEP_Mean_Site<-ddply(NEP, .(Site_ID, Type_Flux),
                          summarise,
                          Stand_Age= mean(Stand_Age, na.rm=T),
                          values=mean(values, na.rm=T),
                          Annual_Preci=mean(Annual_Preci, na.rm=T),
                          Tair=mean(Tair, na.rm=T))

GPP_Mean_Site<-ddply(GPP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Annual_Preci=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))

Reco_Mean_Site<-ddply(Reco, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Annual_Preci=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))

Ratio_NEP_GPP_Mean_Site<-ddply(Ratio_NEP_GPP, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Annual_Preci=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))

Ratio_GPP_Reco_Mean_Site<-ddply(Ratio_GPP_Reco, .(Site_ID, Type_Flux),
                     summarise,
                     Stand_Age= mean(Stand_Age, na.rm=T),
                     values=mean(values, na.rm=T),
                     Annual_Preci=mean(Annual_Preci, na.rm=T),
                     Tair=mean(Tair, na.rm=T))

Ratio_NEP_GPPmax_Mean_Site<-ddply(Ratio_NEP_GPPmax, .(Site_ID, Type_Flux),
                                summarise,
                                Stand_Age= median(Stand_Age, na.rm=T),
                                values=mean(values, na.rm=T),
                                Annual_Preci=mean(Annual_Preci, na.rm=T),
                                Tair=mean(Tair, na.rm=T))

# 1.2 Test the best fitting function to ecosystem response - Modelling efficiency (Nash-Sutcliffe Efficiency test)

# Load function 
source("Function/NSE_NEP.R")
source("Function/NSE_GPP.R")
source("Function/NSE_Reco.R")
source("Function/NSE_Ratio_GPP_ER.R")
source("Function/NSE_Ratio_NEP_GPP.R")
source("Function/NSE_Ratio_NEP_GPPmax.R")

# Create a list of dataframe
df.list1<-list(NEP, NEP_Mean_Site)
df.list2<-list(GPP, GPP_Mean_Site)
df.list3<-list(Reco, Reco_Mean_Site)
df.list4<-list(Ratio_GPP_Reco, Ratio_GPP_Reco_Mean_Site)
df.list5<-list(Ratio_NEP_GPP, Ratio_NEP_GPP_Mean_Site)
df.list6<-list(Ratio_NEP_GPPmax, Ratio_NEP_GPPmax_Mean_Site)

Out1<- lapply(df.list1, stat_NEP)
Out2<- lapply(df.list2, stat_GPP)
Out3<- lapply(df.list3, stat_Reco)
Out4<- lapply(df.list4, stat_Ratio_GPP_ER)
Out5<- lapply(df.list5, stat_Ratio_NEP_GPP)
Out6<- lapply(df.list6, stat_Ratio_NEP_GPPmax)

# 2. Plot ecosystem response with the best fit function with all years per site
source("Function/CI_Est.R")

# 2.1 NEP

# Compute the best fit function
Fun_NEP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = NEP,
               start = list(A= 2.734e+02, B=-4.834e-03, C=-9.307e+02, D=-1.582e-01), control = list(maxiter = 500))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP$Stand_Age),max(NEP$Stand_Age),length=1000)
pred1 <- approx(NEP$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg1<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(NEP, !is.na(Annual_Preci)), aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(NEP, is.na(Annual_Preci)), aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 2.2 GPP

# Compute the best fit function
Fun_GPP<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = GPP, 
               start = list(A=1327.462, k=-0.134), control = list(maxiter = 500))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_GPP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(GPP$Stand_Age),max(GPP$Stand_Age),length=1000)
pred1 <- approx(GPP$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(GPP$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(GPP$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg2<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(GPP, !is.na(Annual_Preci)), aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(GPP, is.na(Annual_Preci)), aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(0,4000)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 2.3 Reco

# Compute the best fit function
Fun_Reco<-nlsLM(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = Reco, 
                start = list(A = 407.499524, B =0.379616, k =-0.005368), control = list(maxiter = 500))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Reco), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Reco$Stand_Age),max(Reco$Stand_Age),length=1000)
pred1 <- approx(Reco$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Reco$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Reco$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg3<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(Reco, !is.na(Annual_Preci)), aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Reco, is.na(Annual_Preci)), aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(0,4000)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+ 
  facet_grid(Type_Flux~., scales="free_x")


# 2.4 Ratio GPP/Reco

# Compute the best fit function
Fun_Ratio_GPP_Reco<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = Ratio_GPP_Reco, 
                          start = list(A= 1.1635, k= -0.2284), control = list(maxiter = 500))
# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_GPP_Reco), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_GPP_Reco$Stand_Age),max(Ratio_GPP_Reco$Stand_Age),length=1000)
pred1 <- approx(Ratio_GPP_Reco$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_GPP_Reco$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_GPP_Reco$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg4<-ggplot(predVals, aes(x, lower, upper)) +
  geom_hline(yintercept=1, colour='grey', lty="dashed", size=0.8)+
  geom_point(data = subset(Ratio_GPP_Reco, !is.na(Annual_Preci)), aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_GPP_Reco, is.na(Annual_Preci)), aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Ratio GPP-Reco")+ 
  ylim(0,2)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

# 2.5 Ratio NEP/GPP

# Compute the best fit function
Fun_Ratio_NEP_GPP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPP,
                         start = list(A=0.204818, B= -0.006599, C=-1.411068, D=-0.145227), control = list(maxiter = 500))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPP$Stand_Age),max(Ratio_NEP_GPP$Stand_Age),length=1000)
pred1 <- approx(Ratio_NEP_GPP$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPP$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPP$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

max(predVals$upper)

# Plot using ggplot
gg5<-ggplot(predVals, aes(x, lower, upper)) +
  geom_hline(yintercept=0, colour='grey', lty="dashed", size=0.8)+
  geom_point(data = subset(Ratio_NEP_GPP, !is.na(Annual_Preci)), aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_NEP_GPP, is.na(Annual_Preci)), aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Ratio NEP-GPP")+ 
  ylim(-1.5,1)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

# 2.6 Ratio NEP/GPPclimax 

# Compute the best fit function
Fun_Ratio_NEP_GPPmax<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPPmax,
                            start = list(A=-0.741553, B= -0.116557, C=0.236692, D=-0.006249), control = list(maxiter = 500))
# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPPmax), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPPmax$Stand_Age),max(Ratio_NEP_GPPmax$Stand_Age),length=1000)
pred1 <- approx(Ratio_NEP_GPPmax$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPPmax$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPPmax$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg6<-ggplot(predVals, aes(x, lower, upper)) +
  geom_hline(yintercept=0, colour='grey', lty="dashed", size=0.8)+
  geom_point(data = subset(Ratio_NEP_GPPmax, !is.na(Annual_Preci)),  aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_NEP_GPPmax, is.na(Annual_Preci)),  aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Ratio NEP-GPPclimax")+ 
  # ylim(0,3000)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

# 2.8. Plot all carbon fluxes plots together

# Create an arrange plot object
source("Function/Legend_Grid_Arrange.R")
pdf("Latex/Figures/Annual_Flux_All_Years.eps", width = 15, height = 12) # Open a new pdf file
grid_arrange_shared_legend(gg1, gg2, gg3, gg4, gg5, gg6, nrow=3) # Write the grid.arrange in the file
dev.off() # Close the file

# 3. Plot ecosystem response with the best fit function with the mean per site
source("Function/CI_Est.R")

# 3.1 NEP

# Compute the best fit function
Fun_NEP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = NEP_Mean_Site,
               start = list(A= 2.379e+02, B= -3.295e-03, C=-7.879e+02, D=-1.519e-01), control = list(maxiter = 500))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP_Mean_Site$Stand_Age),max(NEP_Mean_Site$Stand_Age),length=1000)
pred1 <- approx(NEP_Mean_Site$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP_Mean_Site$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP_Mean_Site$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg1<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(NEP_Mean_Site, !is.na(Annual_Preci)),  aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(NEP_Mean_Site, is.na(Annual_Preci)),  aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  # ylim(0,3000)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")


# 3.2 GPP

# Compute the best fit function
Fun_GPP<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = GPP_Mean_Site, 
              start = list(A=1287.1816, k= -0.1344), control = list(maxiter = 500))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_GPP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(GPP_Mean_Site$Stand_Age),max(GPP_Mean_Site$Stand_Age),length=1000)
pred1 <- approx(GPP_Mean_Site$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(GPP_Mean_Site$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(GPP_Mean_Site$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg2<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(GPP_Mean_Site, !is.na(Annual_Preci)),  aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(GPP_Mean_Site, is.na(Annual_Preci)),  aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(0,4000)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 3.3 Reco Harvest

# Compute the best fit function
Fun_Reco<-nlsLM(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = Reco_Mean_Site, 
                start = list(A = 631.614933, B = 0.154252, k = -0.001269))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Reco), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Reco_Mean_Site$Stand_Age),max(Reco_Mean_Site$Stand_Age),length=1000)
pred1 <- approx(Reco_Mean_Site$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Reco_Mean_Site$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Reco_Mean_Site$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg3<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(Reco_Mean_Site, !is.na(Annual_Preci)),  aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Reco_Mean_Site, is.na(Annual_Preci)),  aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(0,4000)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~., scales="free_x")

# 3.4 Ratio GPP/Reco Harvest

# Compute the best fit function
Fun_Ratio_GPP_Reco<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = Ratio_GPP_Reco_Mean_Site, 
                          start = list(A= 1.1582, k= -0.2312), control = list(maxiter = 500))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_GPP_Reco), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_GPP_Reco_Mean_Site$Stand_Age),max(Ratio_GPP_Reco_Mean_Site$Stand_Age),length=1000)
pred1 <- approx(Ratio_GPP_Reco_Mean_Site$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_GPP_Reco_Mean_Site$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_GPP_Reco_Mean_Site$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg4<-ggplot(predVals, aes(x, lower, upper)) +
  geom_hline(yintercept=1, colour='grey', lty="dashed", size=0.8)+
  geom_point(data = subset(Ratio_GPP_Reco_Mean_Site, !is.na(Annual_Preci)),  aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_GPP_Reco_Mean_Site, is.na(Annual_Preci)),  aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Ratio GPP-Reco")+ 
  ylim(0,2)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

# 3.5 Ratio NEP/GPP

# Compute the best fit function
Fun_Ratio_NEP_GPP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPP_Mean_Site,
                         start = list(A=0.165450, B= -0.003772, C=-1.319022, D=-0.148503), control = list(maxiter = 500))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPP_Mean_Site$Stand_Age),max(Ratio_NEP_GPP$Stand_Age),length=1000)
pred1 <- approx(Ratio_NEP_GPP_Mean_Site$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPP_Mean_Site$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPP_Mean_Site$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg5<-ggplot(predVals, aes(x, lower, upper)) +
  geom_hline(yintercept=0, colour='grey', lty="dashed", size=0.8)+
  geom_point(data = subset(Ratio_NEP_GPP_Mean_Site, !is.na(Annual_Preci)),  aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_NEP_GPP_Mean_Site, is.na(Annual_Preci)),  aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Ratio NEP-GPP")+ 
  ylim(-1.2,1)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

# 3.6 Ratio NEP/GPPclimax 

# Compute the best fit function
Fun_Ratio_NEP_GPPmax<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPPmax_Mean_Site,
                            start = list(A=-0.776705, B= -0.161076, C=0.189838, D=-0.002193), control = list(maxiter = 500))

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPPmax), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPPmax_Mean_Site$Stand_Age),max(Ratio_NEP_GPPmax_Mean_Site$Stand_Age),length=1000)
pred1 <- approx(Ratio_NEP_GPPmax_Mean_Site$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPPmax_Mean_Site$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPPmax_Mean_Site$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg6<-ggplot(predVals, aes(x, lower, upper)) +
  geom_hline(yintercept=0, colour='grey', lty="dashed", size=0.8)+
  geom_point(data = subset(Ratio_NEP_GPPmax_Mean_Site, !is.na(Annual_Preci)),  aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_NEP_GPPmax_Mean_Site, is.na(Annual_Preci)),  aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Stand age") + ylab("Ratio NEP-GPPclimax")+ 
  # ylim(0,3000)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

# 3.8.Plot all carbon fluxes plots together

# Create an arrange plot object
source("Function/Legend_Grid_Arrange.R")
pdf("Latex/Figures/Annual_Flux_Mean_Site.eps", width = 15, height = 12) # Open a new pdf file
grid_arrange_shared_legend(gg1, gg2, gg3, gg4, gg5, gg6, nrow=3) # Write the grid.arrange in the file
dev.off() # Close the file

