## The function r2lme computes R^2 statistics from the output estimated with lme function in nlme package (the function is written for version 3.1-155). 
## The function requires nlme package. 
## The function is written for models with random intercept. The model may or may not have random slopes.  
## The function is written for models with constant level-1 residual variance, i.e., e~N(0, sigma^2). 
## The following arguments must be provided.
## model: lme output object. 
## effectf: selected fixed effects for which the effect size is computed (if no fixed effect is selected, effectf<-NULL).
## effectr: selected random effects for which the effect size is computed (if no fixed effect is selected, effectr<-NULL; use "(Intercept)" to specify the random intercept).

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
  ## covariance matrix of predictors (sigma, sigmar) 
  if(length(xnames)==1) {
    if(xnames=="(Intercept)") {
      sigmax<-0
    } else {
      sigmax<-var(data[,xnames])
    }
  }  else if (length(xnames)==2){
    sigmax<-var(data[,xnames[xnames!="(Intercept)"]])
  }  else {
    sigmax<-cov(data[,xnames[xnames!="(Intercept)"]])
  }
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
  value1<-matrix(c(esigma2,varf,varr,vartotal,R2f,R2r,R2e),nrow=1)
  colnames(value1)<-c("Residual","AllFixed","AllRandom","TotalVar","R2AllFixed","R2AllRandom","R2Residual")
  value2<-matrix(c(vartotal,vareff,R2),nrow=1)
  colnames(value2)<-c("TotalVar","EffectVar","R2")
  output<-list(xnames,xrnames,value1,effectf,effectr,value2)
  names(output)<-c("All Fixed Effects in the model","All Random Effects in the model",
                   "Variance partitioning",
                   "Selected Fixed Effects","Selected Random Effects",
                   "R2 for selected effects")
  return(output)
}


