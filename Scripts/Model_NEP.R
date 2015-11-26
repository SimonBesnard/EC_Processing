## Script to model NEP vs. stand age
## Author: Simon Besnard
## 12.10.2015
###################################
## Load the necessary packages

# Import dataframe
dfAll_Sites<-readRDS("Output/df_Annual_Flux.rds")

#1. Subset dataset
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

#2. Model NEP and carbon stocks
source("Function/C_Model.R")
nyears=max(GPP$Stand_Age)
GPP <- GPP[order(GPP$Stand_Age),]
f_Age_GPP<- function (x) {1287.1921959 *(1-exp(-0.1343273*x))}
GPP$f_Age<- f_Age_GPP(GPP$Stand_Age)
GPP_Flux<-GPP$f_Age
plot(GPP$Stand_Age, GPP$f_Age)

##GPP[1:20] <- seq(0,1000,length.out = 20) # or ramp up GPP
NPPeff <- 0.5
kcBio <- rep(0.05, nyears)
kCsoil1 <- rep(0.2, nyears)
kCsoil2 <- rep(0.01, nyears)
h=0.3
#Assume steady state before: mini spin up or analytical (for each pool input/rate_constant)
CbioSS <- mean(GPP_Flux)*NPPeff/mean(kcBio) #not enitrely correct if kCBio varies
Csoil1SS <- CbioSS*mean(kcBio)/mean(kCsoil1)
Csoil2SS <- Csoil1SS*h*mean(kCsoil1)/mean(kCsoil2)

Cdyn1 <- simpleCdyn(GPP_Flux, NPPeff,kcBio , kCsoil1, kCsoil2, h, 0, Csoil1SS, Csoil2SS)

gg1<-ggplot(data = Cdyn1, aes(x = time, y = NEP)) +
  geom_point()+
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
  xlab("Stand Age (year)")+
  ylab("Carbon stock")+
  scale_colour_manual("",breaks = c("Ctot", "Cbio", "Csoil1", "Csoil2"),
                      values=c("darkgreen", "blue","brown", "black"), 
                      labels=c("Total Carbon", "Above-ground", "Below-ground", "Deadwood"))

pdf("Latex/Figures/NEP_Model.eps", width = 8.87, height = 5.48) # Open a new pdf file
grid.arrange(gg1, gg2, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file

