## Script to compute statistical analysis annual carbon flux
## Author: Simon Besnard
## 15.05.2015
###################################
## Load the necessary packages
library (ggplot2)
library(scales)
library (dplyr)
library (plyr)
library(tidyr)
library(gridExtra)
library (car)
library (lme4)
library(relaimpo)

#1 Explain variabilty of the fluxes

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")

#1.1. Subset dataset
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]

#Subset type of flux

#NEP
NEP<-Flux_High[Flux_High$Type_Flux %in% c("NEP"),]

#GPP
GPP<-Flux_High[Flux_High$Type_Flux %in% c("GPP"),]

#Reco
Reco<-Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]

#Ratio GPP-Reco
Ratio_GPP_ER<-Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]

#Ratio NEP-GPP
Ratio_NEP_GPP<-Flux_High[Flux_High$Type_Flux %in% c("NEP_GPP"),]

# 1.2 Analysis for NEP

#Compute transform function
Fun_NEP<-nls(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data = NEP, 
             start = list(A=0.02, B=-0.6, C= 50, D=200))
coef(Fun_NEP)
f_Age_NEP<- function (x) {9.717403e-04 *x^3+-2.853405e-01*x^2+2.322953e+01*x -3.562449e+02}
f_Tair_NEP<- function (x) {3000/(1+(exp(1.315-(0.119*x))))}
f_P_NEP<- function (x) {6.313045e-04*x^3 -1.985430e-01 *x^2+1.825791e+01*x -3.266518e+02}

# Compute forward/backward stepwise regresion for climate variables type of disturbance and stand age
lm1 <- lm(values~ Annual_Preci + f_P_NEP(Annual_Preci) + Tair + f_Tair_NEP(Tair)+ Stand_Age + f_Age_NEP(Stand_Age) + Disturbance, data=NEP)
fwd.NEP<-step(lm1, direction="forward")
summary(fwd.NEP)
summary(fwd.NEP)$r.squared

# Compute relative contribution of predictor variable
NEP$Disturbance<-as.integer(NEP$Disturbance)
bootswiss <- boot.relimp(values~ Annual_Preci + f_P_NEP(Annual_Preci) + Tair + f_Tair_NEP(Tair)+ Stand_Age + f_Age_NEP(Stand_Age) +Disturbance, 
                         data=NEP, b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

pdf("Latex/Figures/VarImp_NEP.eps", width = 15, height = 12) # Open a new pdf file
plot(booteval.relimp(bootswiss, bty="perc",
                     sort = TRUE, level=c(0.8,0.9)), names.abbrev=12, sort=T)
dev.off() # Close the file

# 1.3 Analysis for GPP

#Compute transform function
Fun_GPP<-nls(values~A*Stand_Age^2+B*Stand_Age+C, data = GPP, 
             start = list(A=-0.4, B=50, C= 300))
coef(Fun_GPP)
f_Age_GPP<- function (x) {-0.1245135*x^2+21.5238161*x+419.906537}
f_Tair_GPP<- function (x) {3000/(1+(exp(1.315-(0.119*x))))}
f_P_GPP<- function (x) {6.313045e-04*x^3 -1.985430e-01 *x^2+1.825791e+01*x -3.266518e+02}

# Compute forward/backward stepwise regresion for climate variables type of disturbance and stand age
lm1 <- lm(values~ Annual_Preci + f_P_GPP(Annual_Preci) + Tair + f_Tair_GPP(Tair)+ Stand_Age + f_Age_GPP(Stand_Age) + Disturbance, data=GPP)
fwd.GPP<-step(lm1, direction="forward")
summary(fwd.GPP)
summary(fwd.GPP)$r.squared

# Compute relative contribution of predictor variable
GPP$Disturbance<-as.integer(GPP$Disturbance)
bootswiss <- boot.relimp(values~ Annual_Preci + f_P_GPP(Annual_Preci) + Tair + f_Tair_GPP(Tair)+ Stand_Age + f_Age_GPP(Stand_Age) +Disturbance, 
                         data=GPP, b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

pdf("Latex/Figures/VarImp_GPP.eps", width = 15, height = 12) # Open a new pdf file
plot(booteval.relimp(bootswiss, bty="perc",
                     sort = TRUE, level=c(0.8,0.9)), names.abbrev=12, sort=T)
dev.off() # Close the file

# 1.4 Analysis for Reco

#Compute transform function
Fun_Reco<-nls(values~A*Stand_Age^2+B*Stand_Age+C, data = Reco, 
              start = list(A=-0.4, B=50, C= 300))
coef(Fun_Reco)
f_Age_Reco<- function (x) {-0.07539582*x^2+12.36760548*x+630.93933228}
f_Tair_Reco<- function (x) {3000/(1+(exp(1.315-(0.119*x))))}
f_P_Reco<- function (x) {6.313045e-04*x^3 -1.985430e-01 *x^2+1.825791e+01*x -3.266518e+02}

# Compute forward/backward stepwise regresion for climate variables type of disturbance and stand age
lm1 <- lm(values~ Annual_Preci + f_P_Reco(Annual_Preci) + Tair + f_Tair_Reco(Tair)+ Stand_Age + f_Age_Reco(Stand_Age) + Disturbance, data=Reco)
fwd.Reco<-step(lm1, direction="forward")
summary(fwd.Reco)
summary(fwd.Reco)$r.squared

# Compute relative contribution of predictor variable
Reco$Disturbance<-as.integer(Reco$Disturbance)
bootswiss <- boot.relimp(values~ Annual_Preci + f_P_Reco(Annual_Preci) + Tair + f_Tair_Reco(Tair)+ Stand_Age + f_Age_Reco(Stand_Age) +Disturbance, 
                         data=Reco, b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

pdf("Latex/Figures/VarImp_Reco.eps", width = 15, height = 12) # Open a new pdf file
plot(booteval.relimp(bootswiss, bty="perc",
                     sort = TRUE, level=c(0.8,0.9)), names.abbrev=12, sort=T)
dev.off() # Close the file

# 1.5 Analysis for Ratio GPP-Reco

#Compute transform function
Fun_Ratio_GPP_Reco<-nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = Ratio_GPP_Reco, 
                                  start = list(A = 1000, B = 0.170, k = -0.00295))
coef(Fun_Ratio_GPP_Reco)
f_Age_Ratio_GPP_Reco<- function (x) { 0.387217066*(x^0.328086337)*(exp(-0.004177033 *x))}
f_Tair_Ratio_GPP_Reco<- function (x) {3000/(1+(exp(1.315-(0.119*x))))}
f_P_Ratio_GPP_Reco<- function (x) {6.313045e-04*x^3 -1.985430e-01 *x^2+1.825791e+01*x -3.266518e+02}

# Compute forward/backward stepwise regresion for climate variables type of disturbance and stand age
lm1 <- lm(values~ Annual_Preci + f_P_Ratio_GPP_Reco(Annual_Preci) + Tair + f_Tair_Ratio_GPP_Reco(Tair)+ Stand_Age + f_Age_Ratio_GPP_Reco(Stand_Age) + Disturbance, data=Ratio_GPP_Reco)
fwd.Ratio_GPP_Reco<-step(lm1, direction="forward")
summary(fwd.Ratio_GPP_Reco)
summary(fwd.Ratio_GPP_Reco)$r.squared

# Compute relative contribution of predictor variable
Ratio_GPP_Reco$Disturbance<-as.integer(Ratio_GPP_Reco$Disturbance)
bootswiss <- boot.relimp(values~ Annual_Preci + f_P_Ratio_GPP_Reco(Annual_Preci) + Tair + f_Tair_Ratio_GPP_Reco(Tair)+ Stand_Age + f_Age_Ratio_GPP_Reco(Stand_Age) +Disturbance, 
                         data=Ratio_GPP_Reco, b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

pdf("Latex/Figures/VarImp_Ratio_GPP_Reco.eps", width = 15, height = 12) # Open a new pdf file
plot(booteval.relimp(bootswiss, bty="perc",
                     sort = TRUE, level=c(0.8,0.9)), names.abbrev=12, sort=T)
dev.off() # Close the file

# 1.5 Analysis for Ratio NEP-GPP

#Compute transform function
Fun_Ratio_NEP_GPP<-nls(values ~ SSweibull(Stand_Age, Asym, Drop, lrc, pwr), data=Ratio_NEP_GPP)
coef(Fun_Ratio_NEP_GPP)
f_Age_Ratio_NEP_GPP<- function (x) {0.1142166-2.3041145*exp(-exp(-0.9472308)*x^0.6360160)}
f_Tair_Ratio_NEP_GPP<- function (x) {3000/(1+(exp(1.315-(0.119*x))))}
f_P_Ratio_NEP_GPP<- function (x) {6.313045e-04*x^3 -1.985430e-01 *x^2+1.825791e+01*x -3.266518e+02}

# Compute forward/backward stepwise regresion for climate variables type of disturbance and stand age
lm1 <- lm(values~ Annual_Preci + f_P_Ratio_NEP_GPP(Annual_Preci) + Tair + f_Tair_Ratio_NEP_GPP(Tair)+ Stand_Age + f_Age_Ratio_NEP_GPP(Stand_Age) + Disturbance, data=Ratio_NEP_GPP)
fwd.Ratio_NEP_GPP<-step(lm1, direction="forward")
summary(fwd.Ratio_NEP_GPP)
summary(fwd.Ratio_NEP_GPP)$r.squared

# Compute relative contribution of predictor variable
Ratio_NEP_GPP$Disturbance<-as.integer(Ratio_NEP_GPP$Disturbance)
bootswiss <- boot.relimp(values~ Annual_Preci + f_P_Ratio_NEP_GPP(Annual_Preci) + Tair + f_Tair_Ratio_NEP_GPP(Tair)+ Stand_Age + f_Age_Ratio_NEP_GPP(Stand_Age) +Disturbance, 
                         data=Ratio_NEP_GPP, b = 100,  
                         type = "lmg",
                         rank = TRUE, diff = TRUE, rela = TRUE)
print(booteval.relimp(bootswiss))

pdf("Latex/Figures/VarImp_Ratio_NEP_GPP.eps", width = 15, height = 12) # Open a new pdf file
plot(booteval.relimp(bootswiss, bty="perc",
                     sort = TRUE, level=c(0.8,0.9)), names.abbrev=12, sort=T)
dev.off() # Close the file

# 1.6. Plot relative contribution output

#Subset data
df<-read.csv("Output/VarImp_Flux.csv", header = TRUE)
df_NEP<-df[df$Flux %in% c("NEP"),]
df_GPP<-df[df$Flux %in% c("GPP"),]
df_Reco<-df[df$Flux %in% c("Reco"),]
df_NEP_GPP<-df[df$Flux %in% c("Ratio NEP-GPP"),]
df_GPP_Reco<-df[df$Flux %in% c("Ratio GPP-Reco"),]

# Create plot per flux
gg1<- ggplot(df_NEP, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")

gg2<- ggplot(df_GPP, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")

gg3<- ggplot(df_Reco, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")

gg4<- ggplot(df_GPP_Reco, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")

gg5<- ggplot(df_NEP_GPP, aes(x=reorder(Predictor, -Percentage), y=Percentage)) +
  geom_bar(stat="identity", fill="grey")+
  geom_errorbar(aes(ymin=Lower, ymax=Upper), colour="black", width=.1) +
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  facet_wrap(~Flux, scales="free_x")+
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal")+
  xlab("")

#Plot all plots together
pdf("Latex/Figures/VarImp_Flux.eps", width = 15, height = 12) # Open a new pdf file
grid.arrange(gg1, gg2, gg3, gg4, gg5, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file

# 2. Compute residual of fitted model vs. climate variables for ratio GPP-Reco

# Compute residual
Ratio_GPP_Reco$Res<- residuals(Fun_Ratio_GPP_Reco)

#Plot residuals
R2<-(cor(residuals(Fun_Ratio_GPP_Reco), Ratio_GPP_Reco$Tair))^2
p <- ggplot(data = Ratio_GPP_Reco, aes(x = Tair, y = Res)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  geom_point()+
  theme_bw(base_size = 16, base_family = "Helvetica")+
  xlab("Air temperature")+
  ylab("Residual fitted model")+
  annotate("text", label = "R2 = 0.14", x = -2, y = .4)
ggsave(filename = "Latex/Figures/Res_GPP_Reco_Tair.eps", p)



