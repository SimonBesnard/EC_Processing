simpleCdyn <- function (GPP_Flux, NPPeff, kCbio, kCsoil1, kCsoil2, h, CbioIni, Csoil1Ini, Csoil2Ini ) {
  
  n <- length(GPP_Flux)
  NPP <- rep(NA, n)
  Cbio <- rep(CbioIni, n)
  Csoil1 <- rep(Csoil1Ini, n)
  Csoil2 <- rep(Csoil2Ini, n)
  Ctot <- Cbio + Csoil1 + Csoil2
  NEP <- rep(NA, n)
  NEP_flux <- rep(NA, n)
  
  for (i in 1:(n-1)) {
    NPP[i] <- NPPeff * GPP_Flux[i]
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
