params <-
list(prefix = "pfnm")

## ----prelims,include=FALSE,cache=FALSE----------------------------------------------------------------------------------------------------------
library(tidyverse)
library(pomp)
set.seed(594709947L)


## ----model-construct----------------------------------------------------------------------------------------------------------------------------
source("model_measSIR.R")


## ----partrans-----------------------------------------------------------------------------------------------------------------------------------
measSIR |>
  pomp(partrans=parameter_trans(log=c("Beta","Gamma"),logit=c("Rho","Eta")),
    paramnames=c("Beta","Gamma","Eta","Rho")) -> measSIR


## ----ref-params---------------------------------------------------------------------------------------------------------------------------------
coef(measSIR)


## ----like-optim-1-------------------------------------------------------------------------------------------------------------------------------
neg.ll <- function (par, est) {
  try(
    freeze({
      allpars <- coef(measSIR,transform=TRUE)
      allpars[est] <- par
      theta <- partrans(measSIR,allpars,dir="fromEst")
      pfilter(measSIR,params=theta,Np=2000)
    },
    seed=915909831
    )
  ) -> pf
  if (inherits(pf,"try-error")) 1e10 else -logLik(pf)
}



## ----like-optim-2-eval,include=FALSE------------------------------------------------------------------------------------------------------------
stew(file="like_optim2.rda",{
  ## use Nelder-Mead with fixed RNG seed
  estpars <- c("Beta","Gamma","Eta")
  optim(
    par=coef(measSIR,estpars,transform=TRUE),
    est=estpars,
    fn=neg.ll,
    method="Nelder-Mead",
    control=list(maxit=400,trace=0)
  ) -> fit
  
  mle <- measSIR
  coef(mle,estpars,transform=TRUE) <- fit$par
  coef(mle)
  
  fit$val
  
  lls <- replicate(n=5,logLik(pfilter(mle,Np=20000)))
  ll <- logmeanexp(lls,se=TRUE); ll
})


## ----sims1--------------------------------------------------------------------------------------------------------------------------------------
mle |> simulate(nsim=10,format="data.frame",include.data=TRUE) -> sims


## ----sims1-plot,echo=F--------------------------------------------------------------------------------------------------------------------------
sims |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  guides(color="none")+
  geom_line()

