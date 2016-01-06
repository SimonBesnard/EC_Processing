stat_NEP_GPPmax_Temp <- function(dat) {
  id<-unique(dat$Site_ID)
  Lieth_Temp<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  for (i in id){
    fit1 <- try(nlsLM(values~A/(1+exp(B-C*Tair)), data=dat[dat$Site_ID != i,],
                      start = list(A= 5.913e+04, B=1.345e+01, C= 4.809e-02), control = list(maxiter = 500)), silent=TRUE);
    Lieth_Temp[[i]] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[dat$Site_ID == i,]) else NA;
    fit2 <- try(nlsLM(values~A*Tair^2+B*Tair+C, data = dat[dat$Site_ID != i,], 
                      start = list(A= -0.001161, B=0.035291, C= -0.018313), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[[i]] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[dat$Site_ID == i,])else NA; 
    fit3 <- try(nlsLM(values~A*Tair^3+B*Tair^2+C*Tair+D, data = dat[dat$Site_ID != i,], 
                      start = list(A=-1.384e-05, B=-7.172e-04, C= 3.258e-02, D=-1.939e-02), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[[i]]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[dat$Site_ID == i,]) else NA; 
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(melt(Lieth_Temp)$value, dat$values, na.rm=TRUE), NSE(melt(Second_Poly)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Third_Poly)$value, dat$values, na.rm=TRUE)),
       c(cor(melt(Lieth_Temp)$value, dat$values)^2, cor(melt(Second_Poly)$value, dat$values)^2, 
         cor(melt(Third_Poly)$value, dat$values)^2))
}