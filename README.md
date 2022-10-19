# DIMr2

The functions *r2lmer*, *r2lme*, and *r2lmeRmat* compute $R^2$ statistics using the DIM (dummy indicator matrices) method in longitudinal growth models and multilevel models (Li and Ryu, manuscript under review). The functions are used along with *lme4* or *nlme* functions in R.

Li, Z. and Ryu, E. (manuscript under review). Effect Size Measures for Longitudinal Growth Models.

# Description
The model should be estimated using *lmer* function from *lme4* package (lme4::lmer) before using the function *r2lmer*. The model should be estimated using *lme* function from *nlme* package (nlme::lme) before using the functions *r2lme* and *r2lmeRmat*. The function *r2lmeRmat* may be used for models with one of the following level-1 residual covariance structure: Constant (single constant variance $\sigma^2$), Diag (diagonal structure with heterogeneous variances), CS (compound symmetry), CSH (compound symmetry with heterogeneous variances), AR(1) (first-order autoregressive), ARH(1) (first-order autoregressive with heterogeneous variances). 

All three functions require three input arguments:\
(a) model: *lmer* or *lme* output object for the estimated model,\
(b) effectf: selected fixed effects for which the $R^2$ is computed (effectf <- NULL if no fixed effect is selected), and\
(c) effectr: selected random effects for which the $R^2$ is computed (effectr <- NULL if no random effect is selected; effectr <- “(Intercept)” to select the random intercept). 

# Example

Following are examples of using *r2lmer*, *r2lme*, and *r2lmeRmat*.
## *r2lmer*
```{r,r2lmer, echo=T}
library(lme4)
# Estimate the model: 
m<-lmer(formula=y ~ time+w2+(1+time|id),data=data)

# R^2 for the fixed and random effect of time
effectf<-c("time")
effectr<-c("time")
r2lmer(m,effectf,effectr)

# R^2 for the fixed effect of level-1 covariate w2
effectf<-c("w2")
effectr<-NULL
r2lmer(m,effectf,effectr)
```
## *r2lme*
```{r,r2lme, echo=T}
library(nlme)
# Estimate the model:
m<-lme(y ~ time+w2+z1+time*z1,random = ~ 1+time| id, data =data)

# R^2 for the level-2 random effects (random intercept and random slope of time)
effectf<-NULL
effectr<-c("(Intercept)","time")
r2lme(m,effectf,effectr)

# R^2 for the effects of cross-level interaction time and z1
effectf<-c("time:z1")
effectr<-NULL
r2lme(m,effectf,effectr)
```

## *r2lmeRmat*
```{r,r2lme, echo=T}
library(nlme)
# Estimate the model with constant level-1 residual covariance structure:
m<-lme(y ~ time+w2+z1+time*z1,random = ~ 1+time| id, data =data)

# R^2 for the level-2 random effects (random intercept and random slope of time)
effectf<-NULL
effectr<-c("(Intercept)","time")
r2lmeRmat(m,effectf,effectr)

# Estimate the model with CS level-1 residual covariance structure: 
m_CS<-lme(y ~ time+w2, random = ~ 1+time| id,
    correlation = corCompSymm(form = ~time|id),
    data =data)
    
# R^2 for the fixed effect of time
effectf<-c("time")
effectr<-NULL
r2lmeRmat(m_CS,effectf,effectr)

# Estimate the model with ARH(1) level-1 residual covariance structure:
m_ARH<-lme(y ~ time+w2, random = ~ 1+time| id,
    correlation = corAR1(form = ~ time|id),
    weights = varIdent(form = ~ 1 |time),
    data = data)
    
# R^2 for the fixed effect of time
effectf<-c("time")
effectr<-NULL
r2lmeRmat(m_ARH,effectf,effectr)
```
