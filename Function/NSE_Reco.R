stat_Reco <- function(dat) {
  id<-nrow(dat)
  Gamma<-c()
  Second_Poly<-c()
  Third_Poly<-c()
  Amiro<-c()
  Weibull<-c()
  Coursolle<-c()
  for (i in 1:id){
    fit1 <- try(nlsLM(values~A*(Stand_Age^B)*(exp(k*Stand_Age)), data = dat[-i,], 
                    start = list(A = 631.614933, B = 0.154252, k = -0.001269), control = list(maxiter = 500)), silent=TRUE);
    Gamma[i] <- if (inherits(fit1, "nls")) sim = predict(fit1, newdata=dat[i,]) else NA;
    fit2 <- try(nlsLM(values~A*Stand_Age^2+B*Stand_Age+C, data = dat[-i,], 
                    start = list(A=-0.01064, B=3.47055, C=895.44346), control = list(maxiter = 500)), silent=TRUE);
    Second_Poly[i] <- if (inherits(fit2, "nls")) sim = predict(fit2, newdata=dat[i,])else NA; 
    fit3 <- try(nlsLM(values~A*Stand_Age^3+B*Stand_Age^2+C*Stand_Age + D, data = dat[-i,], 
                    start = list(A= 1.326e-04, B=-6.246e-02, C= 8.412e+00, D=7.970e+02), control = list(maxiter = 500)), silent=TRUE); 
    Third_Poly[i]<- if (inherits(fit3, "nls")) sim = predict(fit3, newdata=dat[i,]) else NA; 
    fit4 <- try(nlsLM(values~A*(1-exp(k*Stand_Age)), data = dat[-i,], 
                    start = list(A=1093.452, k= -0.235), control = list(maxiter = 500)), silent=TRUE); # A=1 for ratio NEP-GPP
    Amiro[i]<- if (inherits(fit4, "nls")) sim = predict(fit4, newdata=dat[i,]) else NA;
    fit5 <- try(nlsLM(values~A-(B*exp(-exp(C)*(Stand_Age^D))), data = dat[-i,],
                      start = list(A=1117.5222, B= 580.1631, C=-1.9583, D=0.8073), control = list(maxiter = 500)), silent=TRUE); 
    Weibull[i]<- if (inherits(fit5, "nls")) sim = predict(fit5, newdata=dat[i,]) else NA;
    fit6 <- try(nlsLM(values~A*(exp(B*Stand_Age)) + C*(exp(D*Stand_Age)), data = dat[-i,],
                      start = list(A=1.081e+03, B= 2.627e-04, C=-5.352e+02, D=-1.037e-01), control = list(maxiter = 500)), silent=TRUE); 
    Coursolle[i]<- if (inherits(fit6, "nls")) sim = predict(fit6, newdata=dat[i,]) else NA;
    
  }
  # list(Gamma, Second_Poly, Third_Poly, Amiro, Asympt)
  c(NSE(Gamma, dat$values, na.rm=TRUE), NSE(Second_Poly, dat$values, na.rm=TRUE), NSE(Third_Poly, dat$values, na.rm=TRUE), 
    NSE(Amiro, dat$values, na.rm=TRUE ), NSE(Weibull, dat$values, na.rm=TRUE ), NSE(Coursolle, dat$values, na.rm=TRUE ))
}