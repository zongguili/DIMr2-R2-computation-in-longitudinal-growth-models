## Pinheiro J, Bates D, R Core Team (2025). nlme: Linear and Nonlinear Mixed Effects Models. doi:10.32614/CRAN.package.nlme, R package version 3.1-168, https://CRAN.R-project.org/package=nlme. 
## The function r2lme() computes R^2 statistics from the output estimated with the `lme` function in the 'nlme' package (the function is written for version 3.1-168). 
## The function requires the 'nlme' package. 
## The function is written for models with random intercept. The model may or may not have random slopes.  
## The following arguments must be provided.
## model: lme output object. 
## effectf: selected fixed effect terms for which the R^2 is computed (if no fixed effect is selected, effectf<-NULL).
## effectr: selected random effect terms for which the R^2 is computed (if no fixed effect is selected, effectr<-NULL; use "(Intercept)" to specify the random intercept).


r2lme<-function(model,effectf,effectr) {
  if (is.null(model$na.action)==TRUE) {
    data<-model$data
  } else {
    data<-model$data
    data<-data[!(rownames(data)) %in% as.vector(model$na.action),]
  }
  data[,"(Intercept)"]<-1
  xnames<-names(model$coefficients$fixed)
  for (xname in xnames) {
    if (grepl(":",xname)==TRUE){
      x1<-sub("\\:.*", "",xname)
      x2<-sub(".*:", "",xname)
      intname<-paste(x1,":",x2,sep="")
      data[,intname]<-data[,x1]*data[,x2]
    } else {
    }
  }
  cid<-names(model$coefficients$random)
  xrnames<-colnames(model$coefficients$random[[cid]])
  ## covariance matrix of predictors (sigmax, sigmaxr) 
  ## sigmax
  if(length(xnames)==1) {
    if(xnames=="(Intercept)") {
      sigmax<-0
    } else {
      sigmax<-var(data[,xnames])
    }
  } else if (length(xnames)==2 & "(Intercept)" %in% xnames){
    sigmax<-as.matrix(var(data[,xnames[xnames!="(Intercept)"]]))
    rownames(sigmax)<-xnames[xnames!="(Intercept)"]
    colnames(sigmax)<-xnames[xnames!="(Intercept)"]
  } else {
    sigmax<-cov(data[,xnames[xnames!="(Intercept)"]])
  }
  ## sigmaxr
  if(length(xrnames)==1) {
    if(xrnames=="(Intercept)") {
      sigmaxr<-0
    } else {
      sigmaxr<-var(data[,xrnames])
    }
  } else{
    sigmaxr<-cov(data[,xrnames])
  }
  ## mean vector of predictors associated with level-2 random effects
  if(length(xrnames)==1){
    mu<-mean(data[,xrnames])
  } else {
    mu<-colMeans(data[,xrnames])
  }
  ## parameter estimates: fixed effect coefficients (gamma)
  gammaall<-model$coefficients$fixed
  gamma<-as.matrix(gammaall[names(gammaall)!="(Intercept)"])
  ## parameter estimates: covariance matrix of level-2 random effects (tau)
  tau<-getVarCov(model)
  ## parameter estimates: variance of level-1 residuals
  esigma2<-(model$sigma)^2
  ## Dummy indicator matrix for fixed effects 
  if(length(gamma)==0) {
    gamma<-0
    DIMf<-0
  } else {
    DIMf<-matrix(0,nrow=length(gamma),ncol=1)
    rownames(DIMf)<-xnames[xnames!="(Intercept)"]
    DIMf[is.element(rownames(DIMf),effectf)]<-1
    if (nrow(DIMf)>1) {
      DIMf<-diag(as.vector(DIMf))
    } else {
      rownames(DIMf)<-NULL
    }
  }
  ## Dummy indicator matrix for random effects 
  DIMr<-matrix(0,nrow=length(xrnames),ncol=1)
  rownames(DIMr)<-xrnames
  DIMr[is.element(rownames(DIMr),effectr)]<-1
  if (nrow(DIMr)>1) {
    DIMr<-diag(as.vector(DIMr))
  } else {
    rownames(DIMr)<-NULL
  }
  ## variance 
  vartotal<-t(gamma)%*%sigmax%*%gamma+sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu+esigma2
  varf<-t(gamma)%*%sigmax%*%gamma
  varr<-sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu
  vareff<-(t(DIMf%*%gamma)%*%sigmax%*%(DIMf%*%gamma)+(sum(diag((DIMr%*%tau%*%DIMr)%*%sigmaxr))+t(mu)%*%(DIMr%*%tau%*%DIMr)%*%mu))
  ## R^2 
  R2f<-varf/vartotal
  R2r<-varr/vartotal
  R2e<-esigma2/vartotal
  R2<-vareff/vartotal
  value1<-c(varf,varr,esigma2,vartotal,R2f,R2r,R2e)
  names(value1)<-c("AllFixed","AllRandom","Residual","TotalVar","R2AllFixed","R2AllRandom","R2Residual")
  value2<-c(vartotal,vareff,R2)
  names(value2)<-c("TotalVar","EffectVar","R2")
  ## warning for correlated terms
  fxrnames<-colnames(sigmax)[(colnames(sigmax) %in% effectf)]
  fxcnames<-colnames(sigmax)[!(colnames(sigmax) %in% effectf)]
  fx<-sigmax[fxrnames,fxcnames]
  rx<-tau[rownames(tau)==effectr,colnames(tau)!=effectr]
  if (!is.null(fx) && !all(fx==0)) {
    warnf<-"The selected fixed effect term(s) are correlated with other terms."
  } else {
    warnf<-NA 
  }
  if (!is.null(rx) && !all(rx==0)) {
    warnr<-"The selected random effect term(s) are correlated with other terms."
  } else {
    warnr<-NA 
  }
  warn<-c(warnf,warnr)
  output<-list(xnames,xrnames,value1,effectf,effectr,value2,warn)
  names(output)<-c("All Fixed Effects in the model",
                   "All Random Effects in the model",
                   "Variance partitioning",
                   "Selected Fixed Effects",
                   "Selected Random Effects",
                   "R2 for selected effects",
                   "Warning for correlated terms")
  return(output)
}
