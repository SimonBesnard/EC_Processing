stat_Temp_NEP <- function(dat) {
  id<-unique(dat$Site_ID)
  Lieth_Temp<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  for (i in id){
    fit1 <- try(nlsLM(values~A/(1+exp(B-C*Tair)), data=dat[dat$Site_ID != i,],
                      start = list(A= 300, B=0.5, C=0.119), control = list(maxiter = 500)), silent=TRUE);
    Lieth_Temp[[i]] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[dat$Site_ID == i,]) else NA;
    fit2 <- try(nlsLM(values~A*Tair^2+B*Tair+C, data = dat[dat$Site_ID != i,], 
                      start = list(A=-0.9222 , B=32.9461, C= 6.5651), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[[i]] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[dat$Site_ID == i,])else NA; 
    fit3 <- try(nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = dat[dat$Site_ID != i,], 
                      start = list(A=-0.02505, B=-0.11885, C= 28.03883, D=4.61699), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[[i]]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[dat$Site_ID == i,]) else NA; 
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(melt(Lieth_Temp)$value, dat$values, na.rm=TRUE), NSE(melt(Second_Poly)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Third_Poly)$value, dat$values, na.rm=TRUE)),
       c(cor(melt(Lieth_Temp)$value, dat$values)^2, cor(melt(Second_Poly)$value, dat$values)^2, 
         cor(melt(Third_Poly)$value, dat$values)^2))
}