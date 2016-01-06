stat_NEP_GPPmax_P <- function(dat) {
  id<-unique(dat$Site_ID)
  Lieth_P<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  for (i in id){
    fit1 <- try(nlsLM(values~A*(1-exp(-B*Annual_Preci)), data = dat[dat$Site_ID != i,],
                      start = list(A= 0.153612, B= 0.002419), control = list(maxiter = 500)), silent=TRUE);
    Lieth_P[[i]] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[dat$Site_ID == i,]) else NA;
    fit2 <- try(nlsLM(values~A*Annual_Preci^2+B*Annual_Preci+C, data = dat[dat$Site_ID != i,], 
                      start = list(A=-5.779e-08, B=1.882e-04, C= 1.905e-02), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[[i]] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[dat$Site_ID == i,])else NA; 
    fit3 <- try(nlsLM(values~A*Annual_Preci^3+B*Annual_Preci^2+C*Annual_Preci+D, data = dat[dat$Site_ID != i,], 
                      start = list(A= 2.840e-11, B=-1.867e-07, C=  3.377e-04, D=-2.424e-02), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[[i]]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[dat$Site_ID == i,]) else NA; 
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(melt(Lieth_P)$value, dat$values, na.rm=TRUE), NSE(melt(Second_Poly)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Third_Poly)$value, dat$values, na.rm=TRUE)),
       c(cor(melt(Lieth_P)$value, dat$values)^2, cor(melt(Second_Poly)$value, dat$values)^2, 
         cor(melt(Third_Poly)$value, dat$values)^2))
}