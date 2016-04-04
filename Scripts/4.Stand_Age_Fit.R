## Script to fit ecosystem carbon flux response with stand age
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
library(gridExtra)
library (ggplot2)

#1.Function fit choice for ecosystem response vs. stand age

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
NEP<-readRDS("Output/NEP.rds")
NEP_Mean_Site<-readRDS("Output/NEP_Mean_Site.rds")
Ratio_NEP_GPP<-readRDS("Output/Ratio_NEP_GPP.rds")
Ratio_NEP_GPP_Mean_Site<- readRDS("Output/Ratio_NEP_GPP_Mean_Site.rds")
GPP<-readRDS("Output/GPP.rds")
GPP_Mean_Site<-readRDS("Output/GPP_Mean_Site.rds")

# Add GPP data to dataframe
NEP$GPP<- GPP$values
NEP_Mean_Site$GPP<- GPP_Mean_Site$values
Ratio_NEP_GPP$GPP<- GPP$values
Ratio_NEP_GPP_Mean_Site$GPP<- GPP_Mean_Site$values

# Compute stand age model performance
source("Function/NSE_NEP.R")
stat_NEP(NEP)
stat_NEP(NEP_Mean_Site)

# 2. Plot ecosystem response with stand age
source("Function/CI_Est.R")

# 2.1 NEP

# Compute the best fit function
Fun_NEP<-try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = NEP, 
                   start = list(A=192.93829, k=-0.08976, offset=-700), control = list(maxiter = 500)), silent=TRUE)

# Get CI of bootstrap
Boot_NEP<- nlsBoot(Fun_NEP, niter=50)
Boot_NEP$bootCI

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
  geom_point(data = NEP, aes(x = Stand_Age, y = values, size=MAP_CRU, colour=GPP), inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#e0f3f8", high ="#d73027", space="Lab")+
  labs(colour="GPP [gC.m-2.y-1]", size="Annual precipitation [mm.y-1]")+
  xlab("Age [years]") + ylab("NEP [gC.m-2.y-1]")+ 
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="bottom", 
        legend.box="horizontal",
        legend.text=element_text(size=9))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))
print(gg1)
ggsave("Latex/Figures/NEP_All_Years.eps", width = 12, height = 6)

# 2.2 Ratio NEP/GPP

# Compute the best fit function
Fun_Ratio_NEP_GPP<-nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = Ratio_NEP_GPP, 
                         start = list(A= 0.11795, k= -0.03746, offset= -1.5), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP_GPP<- nlsBoot(Fun_Ratio_NEP_GPP, niter=100)
Boot_NEP_GPP$bootCI
  
# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_Ratio_NEP_GPP), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(Ratio_NEP_GPP$Stand_Age),max(Ratio_NEP_GPP$Stand_Age),length=1000)
pred1 <- approx(Ratio_NEP_GPP$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(Ratio_NEP_GPP$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(Ratio_NEP_GPP$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg2<-ggplot(predVals, aes(x, lower, upper)) +
  geom_hline(yintercept=0, colour='grey', lty="dashed", size=0.8)+
  geom_point(data = subset(Ratio_NEP_GPP, !is.na(MAP_CRU)), aes(x = Stand_Age, y = values, size=MAP_CRU, colour=GPP), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_NEP_GPP, is.na(MAP_CRU)), aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#e0f3f8", high ="#d73027", space="Lab")+
  labs(colour="GPP [gC.m-2.y-1]", size="Annual precipitation [mm.y-1]")+
  xlab("Age [years]") + ylab("CUEe")+ 
  ylim(-1.5,1)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

# 2.3. Plot all carbon fluxes plots together

# Create an arrange plot object
source("Function/Legend_Grid_Arrange.R")
pdf("Latex/Figures/Annual_Flux_All_Years.eps", width = 10, height = 10) # Open a new pdf file
grid_arrange_shared_legend(gg1, gg2, nrow=3) # Write the grid.arrange in the file
dev.off() # Close the file

# 3. Plot ecosystem response with the best fit function with the mean per site
source("Function/CI_Est.R")

# 3.1 NEP

# Compute the best fit function
Fun_NEP_Mean_Site<-try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = NEP_Mean_Site, 
                             start = list(A=192.93829, k=-0.08976, offset=-700), control = list(maxiter = 500)), silent=TRUE)
# Get CI of bootstrap
Boot_NEP_Mean_Site<- nlsBoot(Fun_NEP_Mean_Site, niter=100)
Boot_NEP_Mean_Site$bootCI

# Calculate the confidence interval
predCI <- predict(as.lm.nls(Fun_NEP_Mean_Site), interval = 'confidence', level = 0.95)

# Make the predictions on our defined x
x <- seq(min(NEP_Mean_Site$Stand_Age),max(NEP_Mean_Site$Stand_Age),length=1000)
pred1 <- approx(NEP_Mean_Site$Stand_Age, predCI[, 1], xout = x) ## fitted values
pred2 <- approx(NEP_Mean_Site$Stand_Age, predCI [, 2], xout = x) ## lower CI
pred3 <- approx(NEP_Mean_Site$Stand_Age, predCI[, 3], xout = x) ## upper CI

# Put this into a data frame
predVals <- data.frame(x=x, fit=pred1$y,lower=pred2$y,upper=pred3$y)

# Plot using ggplot
gg1<-ggplot(predVals, aes(x, lower, upper)) +
  geom_point(data = NEP_Mean_Site, aes(x = Stand_Age, y = values, size=MAP_CRU, colour=GPP), inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#e0f3f8", high ="#d73027", space="Lab")+
  labs(colour="GPP [gC.m-2.y-1]", size="Annual precipitation [mm.y-1]")+
  xlab("Age [years]") + ylab("NEP [gC.m-2.y-1]")+ 
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="bottom", 
        legend.box="horizontal",
        legend.text=element_text(size=9))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))
print(gg1)
ggsave("Latex/Figures/NEP_Ave_Site.eps", width = 12, height = 6)

# 3.2 Ratio NEP/GPP

# Compute the best fit function
Fun_Ratio_NEP_GPP<-nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = Ratio_NEP_GPP_Mean_Site, 
                         start = list(A= 0.11795, k= -0.03746, offset= -1.5), control = list(maxiter = 500))

# Get CI of bootstrap
Boot_NEP_GPP_Mean_Site<- nlsBoot(Fun_Ratio_NEP_GPP, niter=50)
Boot_NEP_GPP_Mean_Site$bootCI

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
gg2<-ggplot(predVals, aes(x, lower, upper)) +
  geom_hline(yintercept=0, colour='grey', lty="dashed", size=0.8)+
  geom_point(data = subset(Ratio_NEP_GPP_Mean_Site, !is.na(MAP_CRU)),  aes(x = Stand_Age, y = values, size=MAP_CRU, colour=GPP), inherit.aes = FALSE)+
  geom_point(data = subset(Ratio_NEP_GPP_Mean_Site, is.na(MAP_CRU)),  aes(x = Stand_Age, y = values), colour='grey50', inherit.aes = FALSE)+
  geom_line(aes(y = fit), colour="#666666", size=0.8)+
  geom_line(mapping = aes(y = upper), colour="#666666", lty = "dashed", size=0.8) +
  geom_line(mapping = aes(y = lower), colour="#666666", lty = "dashed", size=0.8)+
  scale_colour_gradient(low="#e0f3f8", high ="#d73027", space="Lab")+
  labs(colour="GPP [gC.m-2.y-1]", size="Annual precipitation [mm.y-1]")+
  xlab("Age [years]") + ylab("CUEe")+ 
  ylim(-1.2,1)+
  scale_size(range = c(4, 9)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="none", 
        legend.box="horizontal",
        legend.text=element_text(size=12))+
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))

# 3.3.Plot all carbon fluxes plots together

# Create an arrange plot object
source("Function/Legend_Grid_Arrange.R")
pdf("Latex/Figures/Annual_Flux_Mean_Site.eps", width = 10, height = 10) # Open a new pdf file
grid_arrange_shared_legend(gg1, gg2, nrow=3) # Write the grid.arrange in the file
dev.off() # Close the file

