nyears=100
GPP <- rep(1000, nyears)
##GPP[1:20] <- seq(0,1000,length.out = 20) # or ramp up GPP
NPPeff <- 0.5
kcBio <- rep(0.05, nyears)
kCsoil1 <- rep(0.2, nyears)
kCsoil2 <- rep(0.01, nyears)
h=0.3
#Assume steady state before: mini spin up or analytical (for each pool input/rate_constant)
CbioSS <- mean(GPP)*NPPeff/mean(kCBio) #not enitrely correct if kCBio varies
Csoil1SS <- CbioSS*mean(kCBio)/mean(kCsoil1)
Csoil2SS <- Csoil1SS*h*mean(kCsoil1)/mean(kCsoil2)

Cdyn1 <- simpleCdyn(GPP, NPPeff,kcBio , kCsoil1, kCsoil2, h, 0, Csoil1SS, Csoil2SS)
par(mfrow=c(2,1), mar=c(2,4,1,1)+0.1)

with(Cdyn1, {
plot(time, NEP, ylim=c(-300,300), lwd=2, type="l")
points(time, NEP_flux, col="blue") ## just to check the match ==> mass balance ok...?
plot(time, Ctot, ylim=c(0, max(Ctot)), type="l", lwd=2)
lines(time, Cbio, col="green", lwd=2)
lines(time, Csoil1, col="purple", lwd=2)
lines(time, Csoil2, col="brown", lwd=2)

})



simpleCdyn <- function (GPP, NPPeff, kCbio, kCsoil1, kCsoil2, h, CbioIni, Csoil1Ini, Csoil2Ini ) {
  
  n <- length(GPP)
  NPP <- rep(NA, n)
  Cbio <- rep(CbioIni, n)
  Csoil1 <- rep(Csoil1Ini, n)
  Csoil2 <- rep(Csoil2Ini, n)
  Ctot <- Cbio + Csoil1 + Csoil2
  NEP <- rep(NA, n)
  NEP_flux <- rep(NA, n)
  
  for (i in 1:(n-1)) {
  NPP[i] <- NPPeff * GPP[i]
  Cbio[i+1]  <- Cbio[i] + NPP[i] - kCbio[i]*Cbio[i]
  Csoil1[i+1] <- Csoil1[i] - kCsoil1[i]*Csoil1[i] +  kCbio[i]*Cbio[i]
  Csoil2[i+1] <- Csoil2[i] - kCsoil2[i]*Csoil2[i] + h*kCsoil1[i]*Csoil1[i]
  
  Ctot[i+1] <- Cbio[i+1]+Csoil1[i+1]+Csoil2[i+1]
  NEP[i] <- Ctot[i+1]-Ctot[i] #Stock based
  NEP_flux[i] <- NPP[i] -  kCsoil2[i]*Csoil2[i] - (1-h)*kCsoil1[i]*Csoil1[i] #Flux based (should be the same)
  
  }
  allinfo=data.frame(time=1:n, NEP=NEP, NEP_flux=NEP_flux, Ctot=Ctot, Cbio=Cbio, Csoil1=Csoil1, Csoil2=Csoil2)
  allinfo
}

