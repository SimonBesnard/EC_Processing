stat_NEP_GPP_Photo <- function(dat) {
  id<-unique(dat$Site_ID)
  Gamma<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  Amiro<-c()
  Chen<-c()
  for (i in id){
    fit1 <- try(nlsLM(values~A*(GPP^B)*(exp(k*GPP)), data = dat[dat$Site_ID != i,], 
                      start = list(A = 3.852e-09, B =2.667e+00, k =-1.278e-03), control = list(maxiter = 500)), silent=TRUE);
    Gamma[[i]] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[dat$Site_ID == i,]) else NA;
    fit2 <- try(nlsLM(values~A*GPP^2+B*GPP +C, data = dat[dat$Site_ID != i,], 
                      start = list(A=-1.337e-07, B= 5.752e-04, C=-3.631e-01), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[[i]] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[dat$Site_ID == i,])else NA; 
    fit3 <- try(nlsLM(values~A*GPP^3+B*GPP^2+C*GPP, data = dat[dat$Site_ID != i,], 
                      start = list(A= -4.880e-11, B=2.013e-07, C= -8.399e-05), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[[i]]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[dat$Site_ID == i,]) else NA; 
    fit4 <- try(nlsLM(values~A*(1-exp(k*GPP)), data = dat[dat$Site_ID != i,], 
                      start = list(A=0.3152035, k=-0.0003991), control = list(maxiter = 500)), silent=TRUE);
    Amiro[[i]]<- if (inherits(fit4, "nls")) sim = predict(fit4, newdata=dat[dat$Site_ID == i,]) else NA;
    fit5 <- try(nlsLM(values~A*(1+((B*((GPP/C)^D)-1)/(exp(GPP/C)))), data = dat[dat$Site_ID != i,], 
                      start = list(A = 0.188, B =  -4.871, C =  365.122, D=  -0.393), control = list(maxiter = 500)), silent=TRUE); 
    Chen[[i]]<- if (inherits(fit5, "nls")) sim = predict(fit5, newdata=dat[dat$Site_ID == i,]) else NA;
    
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  list(c(NSE(melt(Gamma)$value, dat$values, na.rm=TRUE), NSE(melt(Second_Poly)$value, dat$values, na.rm=TRUE), 
         NSE(melt(Third_Poly)$value, dat$values, na.rm=TRUE), NSE(melt(Amiro)$value, dat$values, na.rm=TRUE ), 
         NSE(melt(Chen)$value, dat$values, na.rm=TRUE )),
       c(cor(melt(Gamma)$value, dat$values)^2, cor(melt(Second_Poly)$value, dat$values)^2, 
         cor(melt(Third_Poly)$value, dat$values)^2, cor(melt(Amiro)$value, dat$values)^2, 
         cor(melt(Chen)$value, dat$values)^2))
}