## ----prelims,echo=F,cache=F---------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(pomp)
options(stringsAsFactors=FALSE)
set.seed(1350254336)
knitr::opts_chunk$set(highlight=FALSE)


## ----Q6-----------------------------------------------------------------------------------------------------------------------------------------
  ll <- c(-2446,-2444,-2443,-2442,-2440)
  mean(ll)
  sd(ll)
  sd(ll)/sqrt(length(ll))
  library(pomp)
  logmeanexp(ll,se=TRUE)


## ----ebolaModel,echo=TRUE,eval=FALSE------------------------------------------------------------------------------------------------------------
## d <- dacca(deltaI=0.08)


## ----Q9-----------------------------------------------------------------------------------------------------------------------------------------
  d <- dacca(deltaI=0.08)
  library(doFuture)
  plan(multisession)

  bake(file="Q9.rds",{
    foreach(i=1:32,.combine=c,
      .options.future=list(seed=TRUE)
    ) %dofuture% {
      library(pomp)
      logLik(pfilter(d,Np=10000))
    }
  }) -> cholera_loglik
  logmeanexp(cholera_loglik,se=TRUE)

