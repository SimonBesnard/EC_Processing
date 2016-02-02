stat_P_NEP <- function(dat) {
  id<-unique(dat$Site_ID)
  Lieth_P<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  for (i in id){
    fit1 <- try(nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = dat[dat$Site_ID != i,],
                      start = list(A= 300, B=0.000664), control = list(maxiter = 500)), silent=TRUE);
    Lieth_P[[i]] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[dat$Site_ID == i,]) else NA;
    fit2 <- try(nlsLM(values~A*Annual_Preci^2+B*Annual_Preci+C, data = dat[dat$Site_ID != i,], 
                      start = list(A=-8.514e-05 , B=3.341e-01, C= -4.485e+01), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[[i]] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[dat$Site_ID == i,])else NA; 
    fit3 <- try(nlsLM(values~A*Annual_Preci^3+B*Annual_Preci^2+C*Annual_Preci+D, data = dat[dat$Site_ID != i,], 
                      start = list(A=3.442e-08, B=-2.414e-04, C= 5.153e-01, D=-9.730e+01), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[[i]]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[dat$Site_ID == i,]) else NA; 
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(melt(Lieth_P)$value, dat$values, na.rm=TRUE), NSE(melt(Second_Poly)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Third_Poly)$value, dat$values, na.rm=TRUE)),
       c(cor(melt(Lieth_P)$value, dat$values)^2, cor(melt(Second_Poly)$value, dat$values)^2, 
         cor(melt(Third_Poly)$value, dat$values)^2))
}