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
library (JBTools)

#.1.Function fit choice for ecosystem response

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")

# 1.1. Subset dataframe for fitting process 
Flux_High<-dfAll_Sites[dfAll_Sites$Int_Replacement %in% c("High"),]
GPP_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
Reco_High<-Flux_High[Flux_High$Type_Flux %in% c("Respiration"),]
Ratio_High<-Flux_High[Flux_High$Type_Flux %in% c("GPP_ER"),]
NEP_High<-Flux_High[Flux_High$Type_Flux %in% c("NEP"),]

# Harvest
GPP_High_Harvest<-GPP_High[GPP_High$Disturbance %in% c("Harvest"),]
Reco_High_Harvest<-Reco_High[Reco_High$Disturbance %in% c("Harvest"),]
NEP_High_Harvest<-NEP_High[NEP_High$Disturbance %in% c("Harvest"),]
Ratio_High_Harvest<-Ratio_High[Ratio_High$Disturbance %in% c("Harvest"),]

#Fire
GPP_High_Fire<-GPP_High[GPP_High$Disturbance %in% c("Wildfire"),]
Reco_High_Fire<-Reco_High[Reco_High$Disturbance %in% c("Wildfire"),]
Ratio_High_Fire<-Ratio_High[Ratio_High$Disturbance %in% c("Wildfire"),]
NEP_High_Fire<-NEP_High[NEP_High$Disturbance %in% c("Wildfire"),]

# 1.2 Test the best fitting function to ecosystem response

# 1.2.1. Akaike Information Criterion

#Compute the AIC test
stat <- function(dat, inds) { 
  fit1 <- try(nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = dat[inds,], start = list(A = 1000, B = 0.170, k = -0.00295)), silent=TRUE); 
  Gamma <- if (inherits(fit1, "nls")) AIC(fit1) else NA;
  fit2 <- try(nls(values~A*Stand_Age^2+B*Stand_Age+C, data = dat[inds,], start = list(A=-0.4, B=50, C= 300)), silent=TRUE); 
  Second_Poly <- if (inherits(fit2, "nls")) AIC(fit2) else NA; 
  fit3 <- try(nls(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data = dat[inds,], start = list(A=0.02, B=-0.6, C= 50, D=200)), silent=TRUE); 
  Third_Poly <- if (inherits(fit3, "nls")) AIC(fit3) else NA;
  fit4 <- try(nls(values~A*(1-exp(k*Stand_Age)), data = dat[inds,], start = list(A=1500, k= -0.224)), silent=TRUE); 
  Amiro <- if (inherits(fit4, "nls")) AIC(fit4) else NA;
    c(Gamma, Second_Poly, Third_Poly, Amiro) 
} 

df.list<-list(NEP_High_Harvest, GPP_High_Harvest, Reco_High_Harvest, Ratio_High_Harvest, NEP_High_Fire, GPP_High_Fire, Reco_High_Fire, Ratio_High_Fire)
res<-lapply(df.list, function(x) boot(x, stat, R=200))

# Create a dataframe with the output results
df.stat<-list()
for (i in 1:length(res)){
  df.stat[[i]]<-res[[i]]$t
}
names(df.stat)<-c("NEP-Harvest","GPP-Harvest", "Reco-Harvest", "Ratio-Harvest","NEP-Fire","GPP-Fire", "Reco-Fire", "Ratio-Fire")
df.stat<-lapply(df.stat, "colnames<-", paste0(c("Gamma", "Second Polymonial", "Third Polymonial", "Amiro")))
df.stat<-ldply(df.stat, .id="Flux")
df.stat<-gather(df.stat, values, Function, -Flux)

#Plot this output of the test
gg1<-ggplot(df.stat, aes(x =variable, y = value)) +
  geom_boxplot() +
  facet_wrap(~ Flux, scales = "free", nrow=2)+
  xlab("") + ylab("AIC")+
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(gg1)
ggsave("Latex/Figures/AIC_results.eps", height = 12, width = 15)
ggsave("Latex/Figures/AIC_results.pdf", height = 12, width = 15)


# 1.2.1. Modelling efficiency

#Compute the Nash-Sutcliffe Efficiency test
stat <- function(dat, inds) { 
  fit1 <- try(nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = dat[inds,], start = list(A = 1000, B = 0.170, k = -0.00295)), silent=TRUE); 
  Gamma <- if (inherits(fit1, "nls")) MEF(prediction = predict(fit1), observation = dat$values) else NA;
  fit2 <- try(nls(values~A*Stand_Age^2+B*Stand_Age+C, data = dat[inds,], start = list(A=-0.4, B=50, C= 300)), silent=TRUE); 
  Second_Poly <- if (inherits(fit2, "nls")) MEF(prediction = predict(fit2), observation = dat$values) else NA; 
  fit3 <- try(nls(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data = dat[inds,], start = list(A=0.02, B=-0.6, C= 50, D=200)), silent=TRUE); 
  Third_Poly <- if (inherits(fit3, "nls")) MEF(prediction = predict(fit3), observation = dat$values) else NA;
  fit4 <- try(nls(values~A*(1-exp(k*Stand_Age)), data = dat[inds,], start = list(A=1500, k= -0.224)), silent=TRUE); 
  Amiro <- if (inherits(fit4, "nls")) MEF(prediction = predict(fit4), observation = dat$values) else NA;
  c(Gamma, Second_Poly, Third_Poly, Amiro) 
} 

df.list<-list(NEP_High_Harvest, GPP_High_Harvest, Reco_High_Harvest, Ratio_High_Harvest, NEP_High_Fire, GPP_High_Fire, Reco_High_Fire, Ratio_High_Fire)
res<-lapply(df.list, function(x) boot(x, stat, R=200))

# Create a dataframe with the output results
df.stat<-list()
for (i in 1:length(res)){
  df.stat[[i]]<-res[[i]]$t
}
names(df.stat)<-c("NEP-Harvest","GPP-Harvest", "Reco-Harvest", "Ratio-Harvest","NEP-Fire","GPP-Fire", "Reco-Fire", "Ratio-Fire")
df.stat<-lapply(df.stat, "colnames<-", paste0(c("Gamma", "Second Polymonial", "Third Polymonial", "Amiro")))
df.stat<-ldply(df.stat, .id="Flux")
df.stat<-gather(df.stat, values, Function, -Flux)

#Plot this output of the test
gg2<-ggplot(df.stat, aes(x =variable, y = -value)) +
  geom_boxplot() +
  facet_wrap(~ Flux, scales = "free", nrow=2)+
  xlab("") + ylab("Modelling Efficiency")+
  theme_bw(base_size = 14, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(gg2)
ggsave("Latex/Figures/MEF_results.eps", height = 12, width = 15)
ggsave("Latex/Figures/MEF_results.pdf", height = 12, width = 15)

# 2. Plot ecosystem response with the best fit function
source("Function/CI_Est.R")

# 2.1 NEP Harvest

# Compute the best fit function
Fun_NEP_Harvest<-nls(values~A*Stand_Age^2+B*Stand_Age+C, data=NEP_High_Harvest,
                     start = list(A=-0.4, B=50, C= 300))

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
  geom_line(aes(y = fit), colour="black", linetype="dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA,alpha=0.2)+
  geom_point(data = NEP_High_Harvest, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), alpha=0.7, inherit.aes = FALSE)+
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
  facet_grid(Type_Flux~Disturbance, scales="free_x")
print(gg3)

# 2.2 GPP Harvest

# Compute the best fit function
Fun_GPP_Harvest<-nls(values~A*Stand_Age^2+B*Stand_Age+C, data=GPP_High_Harvest,
                     start = list(A=-0.4, B=50, C= 300))

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
  geom_line(aes(y = fit), colour="black", linetype="dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA,alpha=0.2)+
  geom_point(data = GPP_High_Harvest, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), alpha=0.7, inherit.aes = FALSE)+
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
  facet_grid(Type_Flux~Disturbance, scales="free_x")
print(gg4)

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
  geom_line(aes(y = fit), colour="black", linetype="dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA,alpha=0.2)+
  geom_point(data = Reco_High_Harvest, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), alpha=0.7, inherit.aes = FALSE)+
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
  facet_grid(Type_Flux~Disturbance, scales="free_x")
print(gg5)

# 2.4 Ratio GPP/Reco Harvest

# Compute the best fit function
Fun_Ratio_Harvest<-nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data=Ratio_High_Harvest, 
                      start = list(A=1000, B=0.170, k= -0.00295))

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
  geom_line(aes(y = fit), colour="black", linetype="dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA,alpha=0.2)+
  geom_point(data = Ratio_High_Harvest, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), alpha=0.7, inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(0,2)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  geom_hline(yintercept=1, linetype="dashed", colour="grey", size=0.8)+
  facet_grid(Type_Flux~Disturbance, scales="free_x")
print(gg6)

# 2.5 NEP Fire

# Compute the best fit function
Fun_NEP_Fire<-nls(values ~ A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data=NEP_High_Fire, 
                  start = list(A=0.004, B=-0.7, C= 25, D=300))

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
  geom_line(aes(y = fit), colour="black", linetype="dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA,alpha=0.2)+
  geom_point(data = NEP_High_Fire, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), alpha=0.7, inherit.aes = FALSE)+
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
  facet_grid(Type_Flux~Disturbance, scales="free_x")
print(gg7)

# 2.6 GPP Fire

# Compute the best fit function
Fun_GPP_Fire<-nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data=GPP_High_Fire, 
                  start = list(A=1000, B=0.170, k= -0.00295))

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
  geom_line(aes(y = fit), colour="black", linetype="dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA,alpha=0.2)+
  geom_point(data = GPP_High_Fire, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), alpha=0.7, inherit.aes = FALSE)+
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
  facet_grid(Type_Flux~Disturbance, scales="free_x")
print(gg8)

# 2.7 Reco Fire

# Compute the best fit function
Fun_Reco_Fire<-nls(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data=Reco_High_Fire, 
                   start = list(A=1000, B=0.170, k= -0.00295))

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
  geom_line(aes(y = fit), colour="black", linetype="dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA,alpha=0.2)+
  geom_point(data = Reco_High_Fire, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), alpha=0.7, inherit.aes = FALSE)+
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
  facet_grid(Type_Flux~Disturbance, scales="free_x")
print(gg9)

# 2.8 Ratio GPP/Reco Fire

# Compute the best fit function
Fun_Ratio_Fire<-nls(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age+D, data=Ratio_High_Fire,
                    start = list(A=0.02, B=-0.6, C= 50, D=200))

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
  geom_line(aes(y = fit), colour="black", linetype="dashed")+
  geom_ribbon(aes(ymin=lower, ymax=upper), colour=NA,alpha=0.2)+
    geom_point(data = Ratio_High_Fire, aes(x = Stand_Age, y = values, size=Annual_Preci, colour=Tair), alpha=0.7, inherit.aes = FALSE)+
  scale_colour_gradient(low="#00FF33", high ="#FF0000")+
  labs(colour="Annual air temperature (°C)", size="Annual precipitation (mm.y-1)")+
  xlab("Year since disturbance") + ylab("Annual carbon flux (g.m-2.y-1)")+ 
  ylim(0,2)+
  scale_size(range = c(3, 8)) +
  theme_bw(base_size = 12, base_family = "Helvetica") + 
  theme(panel.grid.minor = element_line(colour="grey", size=0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="bottom", 
        legend.box="horizontal") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5),
         size = guide_legend(title.position="top", title.hjust = 0.5))+
  facet_grid(Type_Flux~Disturbance, scales="free_x")
print(gg10)

# 2.9.Plot all carbon fluxes plots together

# Create an arrange plot object
gA <- ggplotGrob(gg3)
gB <- ggplotGrob(gg7)
gC <- ggplotGrob(gg4)
gD <- ggplotGrob(gg8)
gE <- ggplotGrob(gg5)
gF <- ggplotGrob(gg9)
gG <- ggplotGrob(gg6)
gH <- ggplotGrob(gg10)

maxWidth = grid::unit.pmax(gA$widths[1:2], gB$widths[1:2], gC$widths[1:2], gD$widths[1:2], 
                           gE$widths[1:2], gF$widths[1:2], gG$widths[1:2], gH$widths[1:2])
gA$widths[1:2] <- as.list(maxWidth)
gB$widths[1:2] <- as.list(maxWidth)
gC$widths[1:2] <- as.list(maxWidth)
gD$widths[1:2] <- as.list(maxWidth)
gE$widths[1:2] <- as.list(maxWidth)
gF$widths[1:2] <- as.list(maxWidth)
gG$widths[1:2] <- as.list(maxWidth)
gH$widths[1:2] <- as.list(maxWidth)

# PLot 
gg11 <- arrangeGrob(
  gA, gB, gC, gD, gE, gF, gG, gH, nrow = 3, heights = c(0.5, 0.5))
print(gg11)
