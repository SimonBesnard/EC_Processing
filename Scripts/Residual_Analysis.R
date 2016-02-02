## Script to compute residual of stand age model
## Author: Simon Besnard
## 12.10.2015
###################################
## Load the necessary packages
library(ggplot2)
library(gridExtra)

# 1. Compute residual of fitted model vs. climate variables

# 1.1 NEP

# Import dataframe
GPP<-readRDS("Output/GPP.rds")
NEP<-readRDS("Output/NEP.rds")
NEP$GPP<- GPP$values

# Compute stand age function
Fun_NEP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = NEP,
               start = list(A= 2.734e+02, B=-4.834e-03, C=-9.307e+02, D=-1.582e-01), control = list(maxiter = 500))

# Compute residual
NEP$Res<- residuals(Fun_NEP)

#Plot residuals
R2_Tair<-(cor(residuals(Fun_NEP), NEP$Tair, use="pairwise.complete.obs"))^2
R2_Preci<-(cor(residuals(Fun_NEP), NEP$Annual_Preci, use="pairwise.complete.obs"))^2
R2_GPP<-lm(residuals(Fun_NEP)~ NEP$GPP + I(NEP$GPP^2))
summary(R2_GPP)

gg1 <- ggplot(data = NEP, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual NEP model")+
  annotate("text", label = "R-squared = 0.14", x = 20, y = -275) +
  facet_wrap(~Type_Flux)

gg2 <- ggplot(data = NEP, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual NEP model")+
  annotate("text", label = "R-squared = 0.02", x = 2900, y = -275) +
  facet_wrap(~Type_Flux)

gg3<- ggplot(data = NEP, aes(x = GPP, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("GPP (gC.m-2.y-1)")+
  ylab("Residual NEP model")+
  annotate("text", label = "R-squared = 0.25", x = 3200, y = -275) +
  facet_wrap(~Type_Flux)

# 1.2 Ratio GPP-Reco

#Import dataframe
Ratio_GPP_Reco<-readRDS("Output/Ratio_GPP_Reco.rds")
levels(Ratio_GPP_Reco$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-Reco", "Ratio NEP-GPP", "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4" )
Ratio_GPP_Reco$GPP<- GPP$values

#Compute stand age function
Fun_Ratio_GPP_Reco<-nlsLM(values~A*(1-exp(k*Stand_Age)), data = Ratio_GPP_Reco, 
                          start = list(A= 1.1582, k= -0.2312), control = list(maxiter = 500))

# Compute residual
Ratio_GPP_Reco$Res<- residuals(Fun_Ratio_GPP_Reco)

#Plot residuals
R2_Tair<-lm(residuals(Fun_Ratio_GPP_Reco)~ Ratio_GPP_Reco$Tair + I(Ratio_GPP_Reco$Tair^2))
summary(R2_Tair)
R2_Preci<-lm(residuals(Fun_Ratio_GPP_Reco)~ Ratio_GPP_Reco$Annual_Preci + I(Ratio_GPP_Reco$Annual_Preci^2))
summary(R2_Preci)
R2_GPP<- lm(residuals(Fun_Ratio_GPP_Reco)~ Ratio_GPP_Reco$GPP + I(Ratio_GPP_Reco$GPP^2))
summary(R2_GPP)

gg4 <- ggplot(data = Ratio_GPP_Reco, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual GPP-ER model")+
  annotate("text", label = "R-squared = 0.25", x = 20, y = -0.4) +
  facet_wrap(~Type_Flux)

gg5 <- ggplot(data = Ratio_GPP_Reco, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual GPP-ER model")+
  annotate("text", label = "R-squared = 0.06", x = 2900, y = -0.4) +
  facet_wrap(~Type_Flux)

gg6<- ggplot(data = Ratio_GPP_Reco, aes(x = GPP, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("GPP (gC.m-2.y-1)")+
  ylab("Residual GPP-ER model")+
  annotate("text", label = "R-squared = 0.13", x = 3200, y = -0.4) +
  facet_wrap(~Type_Flux)

# 1.3 Ratio NEP-GPP

# Import dataframe
Ratio_NEP_GPP<-readRDS("Output/Ratio_NEP_GPP.rds")
levels(Ratio_NEP_GPP$Type_Flux) <- c("NEE", "GPP", "Respiration", "NEP", "mean_Uncert", "Ratio GPP-Reco", "Ratio NEP-GPP", "x1", "y1", "x2", "y2", "x3", "y3", "x4", "y4")
Ratio_NEP_GPP$GPP<- GPP$values

# Compute stand age function
Fun_Ratio_NEP_GPP<-nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = Ratio_NEP_GPP,
                         start = list(A=0.204818, B= -0.006599, C=-1.411068, D=-0.145227), control = list(maxiter = 500))

# Compute residual
Ratio_NEP_GPP$Res<- residuals(Fun_Ratio_NEP_GPP)

#Plot residuals
R2_Tair<-lm(residuals(Fun_Ratio_NEP_GPP)~ Ratio_NEP_GPP$Tair + I(Ratio_NEP_GPP$Tair^2))
summary(R2_Tair)
R2_Preci<-lm(residuals(Fun_Ratio_NEP_GPP)~ Ratio_NEP_GPP$Annual_Preci + I(Ratio_NEP_GPP$Annual_Preci^2))
summary(R2_Preci)
R2_GPP<- lm(residuals(Fun_Ratio_NEP_GPP)~ Ratio_NEP_GPP$GPP + I(Ratio_NEP_GPP$GPP^2))
summary(R2_GPP)

gg7 <- ggplot(data = Ratio_NEP_GPP, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Air temperature (°C)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.25", x = 20, y = -0.55) +
  facet_wrap(~Type_Flux)

gg8 <- ggplot(data = Ratio_NEP_GPP, aes(x = Annual_Preci, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Precipitation (mm)")+
  ylab("Residual fitted model")+
  annotate("text", label = "R-squared = 0.06", x = 2900, y = -0.55) +
  facet_wrap(~Type_Flux)

gg9<- ggplot(data = Ratio_NEP_GPP, aes(x = GPP, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x + I(x^2))+
  geom_point()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("GPP (gC.m-2.y-1)")+
  ylab("Residual CUE model")+
  annotate("text", label = "R-squared = 0.14", x = 3200, y = -0.55) +
  facet_wrap(~Type_Flux)

# 1.4 Plot all graphs together
pdf("Latex/Figures/Residual_Flux.eps", width = 17, height = 14) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, gg7, gg8, gg9, nrow=3) # Write the grid.arrange in the file
dev.off() # Close the file