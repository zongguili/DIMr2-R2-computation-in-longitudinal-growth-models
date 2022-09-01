## Li, Z. & Ryu, E. (2022). Effect size measures for longitudinal growth models. Unpublished manuscript. 
## The function r2lmer computes R^2 statistics from the output estimated with lmer function in lme4 package (the function is written for version 1.1-28). 
## The function requires lme4 package. 
## The function is written for models with random intercept. The model may or may not have random slopes.  
## The following arguments must be provided.
## model: lmer output object. 
## effectf: selected fixed effects for which the effect size is computed (if no fixed effect is selected, effectf<-NULL).
## effectr: selected random effects for which the effect size is computed (if no fixed effect is selected, effectr<-NULL; use "(Intercept)" to specify the random intercept).


r2lmer<-function(model,effectf,effectr) {
  data<-as.matrix(getME(model,"X"))
  ## fixed effects in the model
  xnames<-colnames(model@pp$X)
  ## random effects in the model
  randomeff<-VarCorr(model)
  cid<-names(randomeff)
  xrnames<-colnames(randomeff[[cid]])
  ## covariance matrix of predictors (sigmax, sigmaxr) 
  if(length(xnames)==1) {
    if(xnames=="(Intercept)") {
      sigmax<-0
    } else {
      sigmax<-var(data[,xnames])
    }
  }  else if (length(xnames)==2){
    sigmax<-var(data[,xnames[xnames!="(Intercept)"]])
    } else {
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
  if(length(xrnames)==1) {
    mu<-mean(data[,xrnames])
  } else {
    mu<-colMeans(data[,xrnames])
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


## Examples

library(lme4)

## Example 1

m1<-lmer(formula=reading~time+WM+sex+time*sex+(1+time|CHILDID),data=dat_wmnona)

m2<-lmer(formula=reading~time+(1+time|CHILDID),data=dat_wmnona)


## example 1-1
effectf<-c("time")
effectr<-c("time")

r2lmer(m2,effectf,effectr)

## example 1-2
effectf<-c("time:sex")
effectr<-NULL

r2lmer(m1,effectf,effectr)

## example 1-3
effectf<-NULL
effectr<-c("(Intercept)")

r2lmer(m1,effectf,effectr)

## example 1-4
effectf<-c("time","WM")
effectr<-NULL

r2lmer(m1,effectf,effectr)


## Example 2

m2<-lmer(formula=reading~time+WM+(1|CHILDID),data=df)

## example 2-1
effectf<-c("time")
effectr<-NULL

r2lmer(m2,effectf,effectr)

## example 2-2
effectf<-NULL
effectr<-c("(Intercept)")

r2lmer(m2,effectf,effectr)


## Example 3 (random intercept model)

m3<-lmer(formula=reading~1+(1|CHILDID),data=df)

effectf<-NULL
effectr<-NULL

r2lmer(m3,effectf,effectr)