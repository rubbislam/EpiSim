params <-
list(prefix = "direct")

## ----prelims,include=FALSE,cache=FALSE----------------------------------------------------------------------------------------------------------
library(tidyverse)
library(pomp)
set.seed(594709947L)


## ----model-construct----------------------------------------------------------------------------------------------------------------------------
source("model_measSIR.R")


## ----bbs-mc-like-2------------------------------------------------------------------------------------------------------------------------------
measSIR |>
  simulate(nsim=5000,format="arrays") -> x
measSIR |>
  emeasure(x=x$states) -> sims
matplot(time(measSIR),t(sims[1,1:50,]),type='l',lty=1,
  xlab="time",ylab=expression(rho*H),bty='l',col='blue')
lines(time(measSIR),obs(measSIR,"reports"),lwd=2,col='black')


## ----bbs-mc-like-3,cache=TRUE-------------------------------------------------------------------------------------------------------------------
ell <- dmeasure(measSIR,y=obs(measSIR),x=x$states,times=time(measSIR),log=TRUE,
  params=coef(measSIR))
dim(ell)


## ----bbs-mc-like-4------------------------------------------------------------------------------------------------------------------------------
ell <- apply(ell,1,sum)
summary(ell)
summary(exp(ell))

