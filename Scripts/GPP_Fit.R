## Script to fit ecosystem carbon flux response with GPP
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library (dplyr)
library (plyr)
library (nlstools)
library(nls2)
library (boot)
library(tidyr)
library (reshape)
library (hydroGOF)
library (minpack.lm)
library (cvTools)
library(grid)
library(ggplot2)
library(gridExtra)

#1. Function fit choice for ecosystem response vs. GPP

# 1.1. Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
NEP<-readRDS("Output/NEP.rds")
NEP_Mean_Site<-readRDS("Output/NEP_Mean_Site.rds")
GPP<-readRDS("Output/GPP.rds")
GPP_Mean_Site<-readRDS("Output/GPP_Mean_Site.rds")
Ratio_GPP_Reco<-readRDS("Output/Ratio_GPP_Reco.rds")
Ratio_GPP_Reco_Mean_Site<-readRDS("Output/Ratio_GPP_Reco_Mean_Site.rds")
Ratio_NEP_GPP<-readRDS("Output/Ratio_NEP_GPP.rds")
Ratio_NEP_GPP_Mean_Site<- readRDS("Output/Ratio_NEP_GPP_Mean_Site.rds")

#1.2 Add GPP data to dataframe
NEP$GPP<- GPP$values
NEP_Mean_Site$GPP<- GPP_Mean_Site$values
Ratio_GPP_Reco$GPP<- GPP$values
Ratio_GPP_Reco_Mean_Site$GPP<- GPP_Mean_Site$values
Ratio_NEP_GPP$GPP<- GPP$values
Ratio_NEP_GPP_Mean_Site$GPP<- GPP_Mean_Site$values

# 1.3 Test the best fitting function to ecosystem response - Modelling efficiency (Nash-Sutcliffe Efficiency test)

# Load function 
source("Function/NSE_NEP_Photo.R")
source("Function/NSE_GPP_ER_Photo.R")
source("Function/NSE_NEP_GPP_Photo.R")

# Create a list of dataframe
df.list1<-list(NEP, NEP_Mean_Site)
df.list2<-list(Ratio_GPP_Reco, Ratio_GPP_Reco_Mean_Site)
df.list3<-list(Ratio_NEP_GPP, Ratio_NEP_GPP_Mean_Site)

# Compute modelling efficiency
Out1<- lapply(df.list1, stat_NEP_Photo)
Out2<- lapply(df.list2, stat_GPP_ER_Photo)
Out3<- lapply(df.list3, stat_NEP_GPP_Photo)

# 2. Plot ecosystem response vs GPP all years per site
source("Function/CI_Est.R")

# 2.1 NEP

# Compute the best fit function
Fun_NEP<-nlsLM(values~A*GPP^2+B*GPP +C, data = NEP, 
               start = list(A=-1.077e-04, B= 5.378e-01, C=-2.785e+02), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP<- nlsBoot(Fun_NEP, niter=100)
Boot_NEP$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP$GPP),max(NEP$GPP),length=100000)
pred1 <- approx(NEP$GPP, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP$GPP, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP$GPP, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg1<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(NEP, !is.na(Annual_Preci)), aes(x = GPP, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(NEP, is.na(Annual_Preci)), aes(x = GPP, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("GPP (g.m-2.y-1)") + ylab("NEP (g.m-2.y-1)")+ 
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

# 2.2 Ratio GPP-ER

# Compute the best fit function
Fun_Ratio_GPP_ER<-nlsLM(values~A*(GPP^B)*(exp(k*GPP)), data = Ratio_GPP_Reco, 
                        start = list(A = 0.0662924, B =0.4576941, k =-0.0002671), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_GPP_ER<- nlsBoot(Fun_Ratio_GPP_ER, niter=100)
Boot_GPP_ER$bootCI
  
# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_GPP_ER), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_GPP_Reco$GPP),max(Ratio_GPP_Reco$GPP),length=100000)
pred1 <- approx(Ratio_GPP_Reco$GPP, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_GPP_Reco$GPP, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_GPP_Reco$GPP, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg2<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(Ratio_GPP_Reco, !is.na(Annual_Preci)), aes(x = GPP, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_GPP_Reco, is.na(Annual_Preci)), aes(x = GPP, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("GPP (g.m-2.y-1)") + ylab("Ratio GPP-ER")+ 
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

# 2.3 Ratio NEP-GPP

# Compute the best fit function
Fun_Ratio_NEP_GPP<-nlsLM(values~A*(1+((B*((GPP/C)^D)-1)/(exp(GPP/C)))), data = Ratio_NEP_GPP, 
                         start = list(A = 0.188, B =  -4.871, C =  365.122, D=  -0.393), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP_GPP<- nlsBoot(Fun_Ratio_NEP_GPP, niter=100)
Boot_NEP_GPP$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPP$GPP),max(Ratio_NEP_GPP$GPP),length=100000)
pred1 <- approx(Ratio_NEP_GPP$GPP, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPP$GPP, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPP$GPP, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg3<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(Ratio_NEP_GPP, !is.na(Annual_Preci)), aes(x = GPP, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_NEP_GPP, is.na(Annual_Preci)), aes(x = GPP, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("GPP (g.m-2.y-1)") + ylab("Ratio NEP-GPP")+ 
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

# 2.4. Plot

# Create an arrange plot object
source("Function/Legend_Grid_Arrange.R")
grid_arrange_shared_legend(gg1, gg2, gg3, nrow=4) # Write the grid.arrange in the file

# 2. Plot ecosystem response vs GPP mean per site
source("Function/CI_Est.R")

# 3.1 NEP

# Compute the best fit function
Fun_NEP_Mean_Site<-nlsLM(values~A*GPP^2+B*GPP +C, data = NEP_Mean_Site, 
               start = list(A=-1.077e-04, B= 5.378e-01, C=-2.785e+02), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP<- nlsBoot(Fun_NEP_Mean_Site, niter=100)
Boot_NEP$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP_Mean_Site), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP_Mean_Site$GPP),max(NEP_Mean_Site$GPP),length=100000)
pred1 <- approx(NEP_Mean_Site$GPP, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP_Mean_Site$GPP, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP_Mean_Site$GPP, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg1<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(NEP_Mean_Site, !is.na(Annual_Preci)), aes(x = GPP, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(NEP_Mean_Site, is.na(Annual_Preci)), aes(x = GPP, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("GPP (g.m-2.y-1)") + ylab("NEP (g.m-2.y-1)")+ 
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

# 3.2 Ratio GPP-ER

# Compute the best fit function
Fun_Ratio_GPP_ER<-nlsLM(values~A*(GPP^B)*(exp(k*GPP)), data = Ratio_GPP_Reco_Mean_Site, 
                        start = list(A = 0.0662924, B =0.4576941, k =-0.0002671), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_GPP_ER<- nlsBoot(Fun_Ratio_GPP_ER, niter=100)
Boot_GPP_ER$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_GPP_ER), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_GPP_Reco_Mean_Site$GPP),max(Ratio_GPP_Reco_Mean_Site$GPP),length=100000)
pred1 <- approx(Ratio_GPP_Reco_Mean_Site$GPP, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_GPP_Reco_Mean_Site$GPP, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_GPP_Reco_Mean_Site$GPP, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg2<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(Ratio_GPP_Reco_Mean_Site, !is.na(Annual_Preci)), aes(x = GPP, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_GPP_Reco_Mean_Site, is.na(Annual_Preci)), aes(x = GPP, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("GPP (g.m-2.y-1)") + ylab("Ratio GPP-ER")+ 
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

# 3.3 Ratio NEP-GPP

# Compute the best fit function
Fun_Ratio_NEP_GPP_Mean_Site<-nlsLM(values~A*GPP^2+B*GPP +C, data = Ratio_NEP_GPP_Mean_Site, 
                                   start = list(A=-1.337e-07, B= 5.752e-04, C=-3.631e-01), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP_GPP<- nlsBoot(Fun_Ratio_NEP_GPP_Mean_Site, niter=100)
Boot_NEP_GPP$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPP_Mean_Site), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPP_Mean_Site$GPP),max(Ratio_NEP_GPP_Mean_Site$GPP),length=100000)
pred1 <- approx(Ratio_NEP_GPP_Mean_Site$GPP, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPP_Mean_Site$GPP, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPP_Mean_Site$GPP, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg3<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = subset(Ratio_NEP_GPP_Mean_Site, !is.na(Annual_Preci)), aes(x = GPP, y = values, size=Annual_Preci, colour=Tair), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_NEP_GPP_Mean_Site, is.na(Annual_Preci)), aes(x = GPP, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#ffffcc", high ="#b10026", space="Lab")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("GPP (g.m-2.y-1)") + ylab("Ratio NEP-GPP")+ 
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

# 3.4. Plot

# Create an arrange plot object
source("Function/Legend_Grid_Arrange.R")
grid_arrange_shared_legend(gg1, gg2, gg3, nrow=4) # Write the grid.arrange in the file

