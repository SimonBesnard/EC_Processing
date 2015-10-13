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


#2. Model NEP and carbon stocs
source("Function/C_Model.R")
nyears=max(GPP$Stand_Age)
GPP <- GPP[order(GPP$Stand_Age),]
GPP_Flux<-GPP$values

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
par(mfrow=c(2,1), mar=c(2,4,1,1)+0.1)

with(Cdyn1, {
plot(time, NEP, ylim=c(-300,300), lwd=2, type="l")
points(time, NEP_flux, col="blue") ## just to check the match ==> mass balance ok...?
plot(time, Ctot, ylim=c(0, max(Ctot)), type="l", lwd=2)
lines(time, Cbio, col="green", lwd=2)
lines(time, Csoil1, col="purple", lwd=2)
lines(time, Csoil2, col="brown", lwd=2)

})


