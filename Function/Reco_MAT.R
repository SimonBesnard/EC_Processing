stat_Reco_Temp <- function(dat) {
  id<-unique(dat$Site_ID)
  Lieth_Temp<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  for (i in id){
    fit1 <- try(nlsLM(values~A/(1+exp(B-C*Tair)), data=dat[dat$Site_ID != i,],
                      start = list(A= 1.966e+06, B=7.890e+00, C=4.827e-02), control = list(maxiter = 500)), silent=TRUE);
    Lieth_Temp[[i]] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[dat$Site_ID == i,]) else NA;
    fit2 <- try(nlsLM(values~A*Tair^2+B*Tair+C, data = dat[dat$Site_ID != i,], 
                      start = list(A= 0.4456, B=51.8351, C= 680.5344), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[[i]] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[dat$Site_ID == i,])else NA; 
    fit3 <- try(nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = dat[dat$Site_ID != i,], 
                      start = list(A=0.4033, B=-12.0200, C= 125.5088, D=712.5760), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[[i]]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[dat$Site_ID == i,]) else NA; 
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(melt(Lieth_Temp)$value, dat$values, na.rm=TRUE), NSE(melt(Second_Poly)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Third_Poly)$value, dat$values, na.rm=TRUE)),
       c(cor(melt(Lieth_Temp)$value, dat$values)^2, cor(melt(Second_Poly)$value, dat$values)^2, 
         cor(melt(Third_Poly)$value, dat$values)^2))
}