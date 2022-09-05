*r2lmer*, *r2lme*, and *r2lmeRmat* are functions to compute the $R^2$ statistic in longitidual growth models (LGMs) or multilevel models (MLMs). These three functions are used along with model estimation functions *lme4* or *nlme* in R.

# Description
Note that to acquire the $R^2$, researchers need to estimate the longitudinal growth model by either *lmer* function from *lme4* package (lme4::lmer) or *lme* function from *nlme* package (nlme::lme) priorly. *r2lmer* and *r2lme* are functions used for computing $R^2$ incorporate constant level-1 residual variance with priorly applied model estimation functions *lmer* and *lme*, respectively. *r2lmeRmat* is used for computing weighted and unweighted $R^2$s incorporate different level-1 residual variance structures with model estimation function *lme*.

The input of the three functions are straightforward which are: a) **model**: model specification that used in the prior model estimation function, b) **effectf**: the variables that corresponding to interested fixed effect, and c) **effectr**: the variables that corresponding to interested random effect. The output of *r2lmer* and *r2lme* functions give the variance partitioning and $R^2$ for selected effects. The function *r2lmeRmat* output the number of observations at each time point, the variance partitioning with and without weights applied, and the $R^2$ for selected effects with and without weights applied. Note that, if homogeneous variance is specified across time (i.e., single constant variance, CS, and AR1), there is no output of weighted variance partitioning and weighted $R^2$.

# Example

Following are two examples of using *r2lmer* and *r2lme*.

```{r,r2lmer, echo=T}
library(lme4)
m<-lmer(formula=y~time+(1+time|id),data=data)
effectf<-c("time")
effectr<-c("time")
r2lmer(m,effectf,effectr)

```

```{r,r2lme, echo=T}
library(nlme)
m<-lme(y~time,random=~1+time|id,data=data)
effectf<-c("time")
effectr<-c("(Intercept)")
r2lme(m,effectf,effectr)

```

