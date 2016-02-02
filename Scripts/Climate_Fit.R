## Script to fit ecosystem carbon flux response with climate condition
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
Ratio_GPP_Reco<-readRDS("Output/Ratio_GPP_Reco.rds")
Ratio_GPP_Reco_Mean_Site<-readRDS("Output/Ratio_GPP_Reco_Mean_Site.rds")
Ratio_NEP_GPP<-readRDS("Output/Ratio_NEP_GPP.rds")
Ratio_NEP_GPP_Mean_Site<- readRDS("Output/Ratio_NEP_GPP_Mean_Site.rds")

# 1.2 Test the best fitting function to ecosystem response - Modelling efficiency (Nash-Sutcliffe Efficiency test)

# Load function 
source("Function/NEP_MAT.R")
source("Function/NEP_GPP_MAT.R")
source("Function/GPP_ER_MAT.R")

# Create a list of dataframe
df.list1<-list(NEP, NEP_Mean_Site)
df.list2<-list(Ratio_GPP_Reco, Ratio_GPP_Reco_Mean_Site)
df.list3<-list(Ratio_NEP_GPP, Ratio_NEP_GPP_Mean_Site)

# Compute modelling efficiency - Temperature
Out1<- lapply(df.list1, stat_Temp_NEP)
Out2<- lapply(df.list2, stat_Temp_GPP_ER)
Out3<- lapply(df.list3, stat_Temp_NEP_GPP)

# 2. Plot ecosystem response vs climate - all years per site
source("Function/CI_Est.R")

# 2.1 NEP vs. temperature

# Compute the best fit function
Fun_NEP<-nlsLM(values~A*Tair^2+B*Tair+C, data = NEP, 
               start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP<- nlsBoot(Fun_NEP, niter=100)
Boot_NEP$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP$Tair),max(NEP$Tair),length=100000)
pred1 <- approx(NEP$Tair, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP$Tair, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP$Tair, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg1<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = NEP, aes(x = Tair, y = values), inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  xlab("Annual air temperature (°C)") + ylab("NEP (g.m-2.y-1)")+ 
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

# 2.2 Ratio GPP-ER

# Compute the best fit function
Fun_Ratio_GPP_ER<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_GPP_Reco, 
                        start = list(A= -0.001214, B=0.036857, C= 0.979401), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_GPP_ER<- nlsBoot(Fun_Ratio_GPP_ER, niter=100)
Boot_GPP_ER$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_GPP_ER), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_GPP_Reco$Tair),max(Ratio_GPP_Reco$Tair),length=100000)
pred1 <- approx(Ratio_GPP_Reco$Tair, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_GPP_Reco$Tair, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_GPP_Reco$Tair, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg2<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = Ratio_GPP_Reco, aes(x = Tair, y = values), inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  xlab("Annual air temperature (°C)") + ylab("Ratio GPP-ER")+ 
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

# 2.3 Ratio NEP-GPP

# Compute the best fit function
Fun_Ratio_NEP_GPP<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPP, 
                         start = list(A= -0.001129, B=0.035511, C= -0.080841), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP_GPP<- nlsBoot(Fun_Ratio_NEP_GPP, niter=100)
Boot_NEP_GPP$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPP$Tair),max(Ratio_NEP_GPP$Tair),length=100000)
pred1 <- approx(Ratio_NEP_GPP$Tair, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPP$Tair, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPP$Tair, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg3<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data =Ratio_NEP_GPP, aes(x = Tair, y = values), inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  xlab("Annual air temperature (°C)") + ylab("Ratio NEP-GPP")+ 
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

# 2.4. Plot

# Create an arrange plot object
pdf("Latex/Figures/Flux_Clim_All_Years.eps", width = 10, height = 10) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, nrow=3) # Write the grid.arrange in the file
dev.off() # Close the file

# 3. Plot ecosystem response vs Temperature - mean per site
source("Function/CI_Est.R")

# 3.1 NEP vs. temperature

# Compute the best fit function
Fun_NEP_Mean_Site<-nlsLM(values~A*Tair^2+B*Tair+C, data = NEP_Mean_Site, 
                 start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP_Mean_Site<- nlsBoot(Fun_NEP_Mean_Site, niter=100)
Boot_NEP_Mean_Site$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP_Mean_Site), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP_Mean_Site$Tair),max(NEP_Mean_Site$Tair),length=100000)
pred1 <- approx(NEP_Mean_Site$Tair, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP_Mean_Site$Tair, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP_Mean_Site$Tair, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg1<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = NEP_Mean_Site, aes(x = Tair, y = values), inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  xlab("Annual air temperature (°C)") + ylab("NEP (g.m-2.y-1)")+ 
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

# 3.2 Ratio GPP-ER

# Compute the best fit function
Fun_Ratio_GPP_ER_Mean_Site<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_GPP_Reco_Mean_Site, 
                                  start = list(A= -0.001214, B=0.036857, C= 0.979401), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_GPP_ER_Mean_Site<- nlsBoot(Fun_Ratio_GPP_ER_Mean_Site, niter=100)
Boot_GPP_ER_Mean_Site$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_GPP_ER_Mean_Site), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_GPP_Reco_Mean_Site$Tair),max(Ratio_GPP_Reco_Mean_Site$Tair),length=100000)
pred1 <- approx(Ratio_GPP_Reco_Mean_Site$Tair, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_GPP_Reco_Mean_Site$Tair, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_GPP_Reco_Mean_Site$Tair, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg2<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = Ratio_GPP_Reco_Mean_Site, aes(x = Tair, y = values), inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  xlab("Annual air temperature (°C)") + ylab("Ratio GPP-ER")+ 
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

# 3.3 Ratio NEP-GPP

# Compute the best fit function
Fun_Ratio_NEP_GPP_Mean_Site<-nlsLM(values~A*Tair^2+B*Tair+C, data = Ratio_NEP_GPP_Mean_Site, 
                         start = list(A= -0.001129, B=0.035511, C= -0.080841), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP_GPP<- nlsBoot(Fun_Ratio_NEP_GPP_Mean_Site, niter=100)
Boot_NEP_GPP$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPP_Mean_Site), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPP_Mean_Site$Tair),max(Ratio_NEP_GPP_Mean_Site$Tair),length=100000)
pred1 <- approx(Ratio_NEP_GPP_Mean_Site$Tair, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPP_Mean_Site$Tair, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPP_Mean_Site$Tair, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg3<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data =Ratio_NEP_GPP_Mean_Site, aes(x = Tair, y = values), inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  xlab("Annual air temperature (°C)") + ylab("Ratio NEP-GPP")+ 
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))

# 3.4. Plot

# Create an arrange plot object
pdf("Latex/Figures/Flux_Clim_Mean_Site.eps", width = 10, height = 10) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, nrow=3) # Write the grid.arrange in the file
dev.off() # Close the file




