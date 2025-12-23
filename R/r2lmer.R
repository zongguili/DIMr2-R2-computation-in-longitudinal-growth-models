## Bates D, Mächler M, Bolker B, Walker S (2015). “Fitting Linear Mixed-Effects Models Using lme4.” Journal of Statistical Software, 67(1), 1–48. doi:10.18637/jss.v067.i01. 
## The function r2lmer() computes R^2 statistics from the output estimated with the `lmer` function in the 'lme4' package (the function is written for version 1.1-38). 
## The function requires the 'lme4' package. 
## The function is written for models with random intercept. The model may or may not have random slopes.  
## The following arguments must be provided.
## model: lmer output object. 
## effectf: selected fixed effect terms for which the R^2 is computed (if no fixed effect is selected, effectf<-NULL).
## effectr: selected random effect terms for which the R^2 is computed (if no fixed effect is selected, effectr<-NULL; use "(Intercept)" to specify the random intercept).


r2lmer<-function(model,effectf,effectr) {
  data<-as.matrix(getME(model,"X"))
  datar<-as.matrix(getME(model,"mmList")[[1]])
  ## fixed effects in the model
  xnames<-colnames(model@pp$X)
  ## random effects in the model
  randomeff<-VarCorr(model)
  cid<-names(randomeff)
  xrnames<-colnames(randomeff[[cid]])
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
      sigmaxr<-var(datar[,xrnames])
    }
  } else {
    sigmaxr<-cov(datar[,xrnames])
  }
  ## mean vector of predictors associated with level-2 random effects
  if(length(xrnames)==1) {
    mu<-mean(datar[,xrnames])
  } else {
    mu<-colMeans(datar[,xrnames])
  }
  ## parameter estimates: fixed effect coefficients (gamma)
  gammaall<-as.matrix(getME(model,"beta"))
  gamma<-as.matrix(gammaall[-1])
  ## parameter estimates: covariance matrix of level-2 random effects (tau)
  allrandom<-as.data.frame(VarCorr(model))
  allrandom$var2<-ifelse(is.na(allrandom$var2)==TRUE,allrandom$var1,
                         allrandom$var2)
  tau<-matrix(nrow=length(xrnames),ncol=length(xrnames))
  for (i in 1:nrow(tau)) {
    for (j in 1:ncol(tau)) {
      if (i>=j) {
        tau[i,j]<-allrandom[allrandom$grp==names(model@cnms) & allrandom$var1==xrnames[j] & allrandom$var2==xrnames[i],"vcov"]
      }
      else {
        tau[i,j]<-allrandom[allrandom$grp==names(model@cnms) & allrandom$var1==xrnames[i] & allrandom$var2==xrnames[j],"vcov"]
      }
    }
  }
  rownames(tau)<-xrnames
  colnames(tau)<-xrnames
  ## parameter estimates: variance of level-1 residuals
  esigma2<-allrandom[allrandom$grp=="Residual","vcov"]
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
  names(output)<-c("All Fixed Effects in the model","All Random Effects in the model",
                   "Variance partitioning",
                   "Selected Fixed Effects","Selected Random Effects",
                   "R2 for selected effects",
                   "Warning for correlated terms")
  return(output)
}
