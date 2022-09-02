## The function r2lmeRmat computes R^2 statistics from the output estimated with lme function in nlme package (the function is written for version 3.1-155). 
## The function requires nlme package. 
## The function is written for models with random intercept. The model may or may not have random slopes.  
## The function may be used with the following specifications for the level-1 residual covariance structure (R matrix): Constant, Diag, CS, CSH, AR1, ARH1. 
## Constant: constant variance
## Diag: diagonal structure with heterogeneous variances
## CS: Compound symmetry
## CSH: Compound symmetry with heterogeneous variances 
## AR1: First-order autoregressive (1)
## ARH1: First-order autoregressive (1) with heterogeneous variances
## The following arguments must be provided.
## model: lme output object. 
## effectf: selected fixed effects for which the effect size is computed (if no fixed effect is selected, effectf<-NULL).
## effectr: selected random effects for which the effect size is computed (if no random effect is selected, effectr<-NULL; use "(Intercept)" to specify the random intercept).

r2lmeRmat<-function(model,effectf,effectr) {
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
  if(length(xnames)==1) {
    if(xnames=="(Intercept)") {
      sigmax<-0
    } else {
      sigmax<-var(data[,xnames])
    }
  } else if (length(xnames)==2){
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
  ## level-1 residual covariance structure
  corSatt<-attributes(model$modelStruct$corStruct)$class[1]
  varSatt<-attributes(model$modelStruct$varStruct)$class[1]
  if (is.null(corSatt)==TRUE && is.null(varSatt)==TRUE) {
    Rmat<-"Constant"
    } else if (is.null(corSatt)==TRUE && varSatt=="varIdent") {
        Rmat<-"Diag"
    } else if (corSatt=="corCompSymm" && is.null(varSatt)==TRUE) {
        Rmat<-"CS"
    } else if (corSatt=="corCompSymm" && varSatt=="varIdent") {
        Rmat<-"CSH"
    } else if (corSatt=="corARMA" && is.null(varSatt)==TRUE) {
        Rmat<-"AR1"
    } else if (corSatt=="corARMA" && varSatt=="varIdent") {
        Rmat<-"ARH1"
    } else {
        Rmat<-NULL
    }
  ## parameter estimates: covariance matrix of level-1 residuals
  if (Rmat=="Constant"){
    esigma2<-(model$sigma)^2
  } else if (Rmat=="Diag") {
    formula<-as.character(attributes(model$modelStruct$varStruct)$formula)[2]
    timevar<-sub(".*\\| ","",formula)
    sigmasqtt<-(coef(model$modelStruct$varStruct,unconstrained=F,allCoef=T)*(model$sigma))^2
    esigma<-diag(sigmasqtt)
  } else if (Rmat=="CS") {
    formula<-as.character(attributes(model$modelStruct$corStruct)$formula)[2]
    timevar<-sub(" \\|.*","",formula)
    rho<-coef(model$modelStruct$corStruct,unconstrained=F)
    csmat<-matrix(rho,nrow=nlevels(factor(data[,timevar])),ncol=nlevels(factor(data[,timevar])))
    diag(csmat)<-1
    esigma<-csmat*(model$sigma)^2
  } else if (Rmat=="CSH") {
    formula<-as.character(attributes(model$modelStruct$varStruct)$formula)[2]
    timevar<-sub(".*\\| ","",formula)
    sigmatt<-coef(model$modelStruct$varStruct,unconstrained=F,allCoef=T)*(model$sigma)
    rho<-coef(model$modelStruct$corStruct,unconstrained=F)
    csmat<-matrix(rho,nrow=nlevels(factor(data[,timevar])),ncol=nlevels(factor(data[,timevar])))
    diag(csmat)<-1
    esigma<-diag(sigmatt)%*%csmat%*%diag(sigmatt)
  } else if (Rmat=="AR1") {
    formula<-as.character(attributes(model$modelStruct$corStruct)$formula)[2]
    timevar<-sub(" \\|.*","",formula)
    phi<-coef(model$modelStruct$corStruct,unconstrained=F)
    ar1cor<-function(n,phi) {
      exponent<-abs(matrix(1:n-1,nrow=n,ncol=n,byrow=TRUE)-(1:n-1))
      phi^exponent
    }
    armat<-ar1cor(nlevels(factor(data[,timevar])),phi)
    esigma<-armat*(model$sigma)^2
  } else if (Rmat=="ARH1") {
    formula<-as.character(attributes(model$modelStruct$varStruct)$formula)[2]
    timevar<-sub(".*\\| ","",formula)
    sigmatt<-coef(model$modelStruct$varStruct,unconstrained=F,allCoef=T)*(model$sigma)
    phi<-coef(model$modelStruct$corStruct,unconstrained=F)
    ar1cor<-function(n,phi) {
      exponent<-abs(matrix(1:n-1,nrow=n,ncol=n,byrow=TRUE)-(1:n-1))
      phi^exponent
    }
    armat<-ar1cor(nlevels(factor(data[,timevar])),phi)
    esigma<-diag(sigmatt)%*%armat%*%diag(sigmatt)
  } else {
    esigma<-NULL
  }
  ## level-1 residual variance unweighted and weighted
  if (Rmat=="Constant"){
    residvar<-esigma2
    residvarw<-esigma2
  } else if (is.null(Rmat)==FALSE) {
    nt<-table(model$data[,timevar])
    nt<-as.vector(nt)
    residvar<-sum(diag(esigma))*(1/nlevels(factor(data[,timevar])))
    ntw<-nt/sum(nt)
    esigmaw<-diag(sqrt(ntw))%*%esigma%*%diag(sqrt(ntw))
    residvarw<-sum(diag(esigmaw))
  }
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
  vartotal<-t(gamma)%*%sigmax%*%gamma+sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu+residvar
  vartotalw<-t(gamma)%*%sigmax%*%gamma+sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu+residvarw
  varf<-t(gamma)%*%sigmax%*%gamma
  varr<-sum(diag(tau%*%sigmaxr))+t(mu)%*%tau%*%mu
  vareff<-(t(DIMf%*%gamma)%*%sigmax%*%(DIMf%*%gamma)+(sum(diag((DIMr%*%tau%*%DIMr)%*%sigmaxr))+t(mu)%*%(DIMr%*%tau%*%DIMr)%*%mu))
  ## R^2 
  R2f<-varf/vartotal
  R2r<-varr/vartotal
  R2e<-residvar/vartotal
  R2<-vareff/vartotal
  R2fw<-varf/vartotalw
  R2rw<-varr/vartotalw
  R2ew<-residvarw/vartotalw
  R2w<-vareff/vartotalw
  value1<-matrix(c(residvar,varf,varr,vartotal,R2f,R2r,R2e),nrow=1)
  colnames(value1)<-c("Residual","AllFixed","AllRandom","TotalVar","R2AllFixed","R2AllRandom","R2Residual")
  value2<-matrix(c(vartotal,vareff,R2),nrow=1)
  colnames(value2)<-c("TotalVar","EffectVar","R2")
  value1w<-matrix(c(residvarw,varf,varr,vartotalw,R2fw,R2rw,R2ew),nrow=1)
  colnames(value1w)<-c("Residual","AllFixed","AllRandom","TotalVar","R2AllFixed","R2AllRandom","R2Residual")
  value2w<-matrix(c(vartotalw,vareff,R2w),nrow=1)
  colnames(value2w)<-c("TotalVar","EffectVar","R2")
  if (Rmat=="Constant"|Rmat=="CS"|Rmat=="AR1") {
    outputuw<-list(xnames,xrnames,Rmat,value1,effectf,effectr,value2)
    names(outputuw)<-c("All Fixed Effects in the model","All Random Effects in the model",
                       "Level-1 residual covariance structure",
                       "Variance partitioning",
                       "Selected Fixed Effects","Selected Random Effects",
                       "R2 for selected effects")
    output<-outputuw
  } else if (Rmat=="Diag"|Rmat=="CSH"|Rmat=="ARH1") {
    outputw<-list(xnames,xrnames,Rmat,nt,value1,value1w,effectf,effectr,value2,value2w)
    names(outputw)<-c("All Fixed Effects in the model","All Random Effects in the model",
                      "Level-1 residual covariance structure",
                      "Number of observations at each time point",
                      "Variance partitioning",
                      "Variance partitioning (weighted residual variance)",
                      "Selected Fixed Effects","Selected Random Effects",
                      "R2 for selected effects",
                      "R2 for selected effects (weighted residual variance)")
    output<-outputw
    }
  return(output)
}


