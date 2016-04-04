stat_NEP <- function(dat) {
  dat$Gamma<-NA
  dat$Amiro<-NA
  dat$Chen<- NA
  dat$Coursolle<- NA
    for(id in unique(dat$Site_ID)){
    fit1 <- try(nlsLM(values~offset + A*(Stand_Age^B)*(exp(k*Stand_Age)), data = dat[dat$Site_ID != id,], 
                    start = list(A = -786.6097, B = 0.3989, k = -0.4694, offset=-700), control = list(maxiter = 500)), silent=TRUE);
    Gamma.pred <- predict(object = fit1, newdata = dat[dat$Site_ID == id,])
    dat$Gamma[dat$Site_ID == id] <- Gamma.pred
    dat$AIC_Gamma[dat$Site_ID == id]<- AIC(fit1)
    fit2 <- try(nlsLM(values~offset + A*(1-exp(k*Stand_Age)), data = dat[dat$Site_ID != id,], 
                    start = list(A=192.93829, k=-0.08976, offset=-700), control = list(maxiter = 500)), silent=TRUE);
    Amiro.pred <- predict(object = fit2, newdata = dat[dat$Site_ID == id,])
    dat$Amiro[dat$Site_ID == id] <- Amiro.pred
    dat$AIC_Amiro[dat$Site_ID == id]<- AIC(fit2)
    fit3 <- try(nlsLM(values~offset + A*(1+((B*((Stand_Age/C)^D)-1)/(exp(Stand_Age/C)))), data = dat[dat$Site_ID != id,], 
                       start = list(A = 2.280, B = 270.887, C = 56.420, D=1.008, offset=-700), control = list(maxiter = 500)), silent=TRUE); 
    Chen.pred <- predict(object = fit3, newdata = dat[dat$Site_ID == id,])
    dat$Chen[dat$Site_ID == id] <- Chen.pred
    dat$AIC_Chen[dat$Site_ID == id]<- AIC(fit3)
    fit4 <- try(nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = dat[dat$Site_ID != id,],
                      start = list(A= 2.867e+02, B=-4.941e-03, C=-1.014e+03, D=-1.802e-01), control = list(maxiter = 500)), silent=TRUE); 
    Coursolle.pred <- predict(object = fit4, newdata = dat[dat$Site_ID == id,])
    dat$Coursolle[dat$Site_ID == id] <- Coursolle.pred
    dat$AIC_Coursolle[dat$Site_ID == id]<- AIC(fit4)
    }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(dat$Gamma, dat$values, na.rm=TRUE), NSE(dat$Amiro, dat$values, na.rm=TRUE), 
         NSE(dat$Chen, dat$values, na.rm=TRUE ), NSE(dat$Coursolle, dat$values, na.rm=TRUE)),
       c(cor(dat$Gamma, dat$values)^2,  cor(dat$Amiro, dat$values)^2, 
         cor(dat$Chen, dat$values)^2, cor(dat$Coursolle, dat$values)^2),
       c(mean(dat$AIC_Gamma, na.rm=TRUE), mean(dat$AIC_Amiro, na.rm=TRUE),
         mean(dat$AIC_Chen, na.rm=TRUE), mean(dat$AIC_Coursolle, na.rm=TRUE)))
}