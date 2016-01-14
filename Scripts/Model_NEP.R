## Script to model NEP vs. stand age
## Author: Simon Besnard
## 12.10.2015
###################################
## Load the necessary packages
library(ggplot2)
library(gridExtra)
library (ggplot2)
library(scales)
library (dplyr)
library (plyr)
library(tidyr)
library(grid)

#1. Model NEP dynamics using multiplicative model

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
GPP$GPP_GPPmax<- NEP$values/GPP$GPPmax
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


# 1.2 Create NEP model using different functions

# 1 2.1 Create daframe with Age and GPPclimax per ecosystem

# Compute GGPclimax per climate
GPPmax<-ddply(GPP, "Climate",
              summarize,
              GPPmax= mean(GPPmax, na.rm=T))

# Compute df
NEP_model<- as.data.frame(0:320)
colnames(NEP_model)<- c("Age")
NEP_model$Af<- 2956.9084
NEP_model$Am<- 3382.6741
NEP_model$Bsh<- 480.4198
NEP_model$Cfa<- 1498.7218
NEP_model$Cfb<- 1356.5680
NEP_model$Cfc<- 1083.7448
NEP_model$Csa<- 1273.2122
NEP_model$Dfb<- 1180.7800
NEP_model$Dfc<- 744.1699
NEP_model$Dwb<- 865.5914
NEP_model<-gather(NEP_model, variable, values, -Age)
colnames(NEP_model)<-c("Age", "Climate", "GPPmax")

# 1.2.2 Compute NEP models

# From NEP to GPP ratio
cols <- c("Af"="#f46d43", "Am"="#d73027", "Bsh"="#000099", "Cfa"="#fdae61", "Cfb"="#fee090", "Cfc"="#abd9e9", "Csa"="#ffffbf", "Dfb"="#e0f3f8",
          "Dfc"="#4575b4", "Dwb"="#74add1")
NEP_model$NEP<- (0.19945826*(exp(-0.00517618*NEP_model$Age)) -1.51419788*(exp(-0.17645323*NEP_model$Age)))*NEP_model$GPPmax
ggplot(data=NEP_model)+
  geom_line(aes(x=Age, y=NEP, colour=Climate))+
  geom_point(data=NEP, aes(x=Stand_Age,y=values, group= NULL),  shape=20, color="black", size=2)+
  theme_bw(base_size = 16, base_family = "Helvetica")+
  xlab("Stand age (year)")+ ylab("NEP (g.m-2.y-1)")+
  theme(axis.text.x = element_text(hjust = 1),
        legend.key = element_blank(),
        legend.position="right", 
        legend.box="horizontal",
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  scale_colour_manual("",breaks = c("Af", "Am", "Bsh", "Cfa", "Cfb", "Cfc", "Csa", "Dfb", "Dfc", "Dwb"),
                      values=cols, 
                      labels=c("Tropical rainforest", "Tropical Monsoon", "Hot steppe", "Humid sub-tropical", "Oceanic-warm summer",
                               "Oceanic-cold summer", "Dry summer", "Warm summer continental-No DS", "Continental subarctic", 
                               "Warm summer continental"))+
  ylim(-700, 1000)

#2. Compute C stocks and NEP dynamics using C pool model

# 2.1 Compute GPP steady state
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")
Flux_High<-dfAll_Sites[dfAll_Sites$Study %in% c("Yes"),]
GPP<- Flux_High[Flux_High$Type_Flux %in% c("GPP"),]
params <- c(GPP15= 1923, GPP1000=1827, a1=242, a2=0.049, k=-0.00025)
GPP$GPPmat<-with(as.list(params), GPP15*((1+exp(a1-a2*15))/(1+exp(a1-a2*GPP$Tair))))
GPP$GPPp<-with(as.list(params), GPP1000*((1-exp(-k*GPP$Annual_Preci))/(1-exp(-k*1000))))
GPP<-transform(GPP, GPPmax = pmin(GPPmat, GPPp))

# 2.2 Create GPP dataframe
source("Function/C_Model.R")
nyears=300
GPP_Flux<- as.data.frame(seq(0:300))
colnames(GPP_Flux)<- "Age"
f_Age_GPP<- function (x) {1287.1921959 *(1-exp(-0.1343273*x))}
GPP_Flux$f_Age<- f_Age_GPP(GPP_Flux$Age)
GPP_Flux<-GPP_Flux$f_Age

# 2.3 Set parameters values of the model
NPPeff <- 0.5
kcBio <- rep(0.05, nyears)
kCsoil1 <- rep(0.2, nyears)
kCsoil2 <- rep(0.01, nyears)
h=0.3
#Assume steady state before: mini spin up or analytical (for each pool input/rate_constant)
CbioSS <- mean(GPP$GPPmax)*NPPeff/mean(kcBio) #not enitrely correct if kCBio varies
Csoil1SS <- CbioSS*mean(kcBio)/mean(kCsoil1)
Csoil2SS <- Csoil1SS*h*mean(kCsoil1)/mean(kCsoil2)
CtotSS<- CbioSS + Csoil1SS +  Csoil2SS

# 2.5 Run the model
Cdyn1 <- simpleCdyn(GPP_Flux, NPPeff,kcBio , kCsoil1, kCsoil2, h, 0, Csoil1SS, Csoil2SS)

# 2.6 Plot outputs
gg1<-ggplot(data = Cdyn1, aes(x = time, y = NEP)) +
  geom_line()+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  xlab("Stand Age (year)")+
  ylab("NEP")

gg2<-ggplot(data = Cdyn1) +
  geom_line(aes(x = time, y = Ctot, colour="Ctot"))+
  geom_line(aes(x = time, y = Cbio, colour="Cbio"))+
  geom_line(aes(x = time, y = Csoil1, colour="Csoil1"))+
  geom_line(aes(x = time, y = Csoil2, colour="Csoil2"))+
  theme_bw(base_size = 12, base_family = "Helvetica")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key = element_blank(),
        legend.position="bottom", 
        legend.box="horizontal",
        legend.text=element_text(size=12),
        axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))+
  xlab("Stand Age (year)")+
  ylab("Carbon stock")+
  scale_colour_manual("",breaks = c("Ctot", "Cbio", "Csoil1", "Csoil2"),
                      values=c("darkgreen", "blue","brown", "black"), 
                      labels=c("Total Carbon", "Above-ground", "Below-ground", "Deadwood"))+ 
  geom_hline(yintercept=11413.33, colour="darkgreen", linetype="dotted")+
  geom_hline(yintercept=31386.67, colour="black", linetype="dotted")+
  geom_hline(yintercept=2853.334, colour="blue", linetype="dotted")+
  geom_hline(yintercept=17120, colour="brown", linetype="dotted")

pdf("Latex/Figures/NEP_Model.eps", width = 8.87, height = 5.48) # Open a new pdf file
grid.arrange(gg1, gg2, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file

