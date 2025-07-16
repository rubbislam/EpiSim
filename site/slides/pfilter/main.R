## ----prelims,echo=FALSE,cache=FALSE---------------------------------------------------------------------------------------------------
library(tidyverse)
library(pomp)
library(foreach)
library(iterators)
library(doFuture)
options(stringsAsFactors=FALSE)
set.seed(1350254336)




## ----model-construct,purl=TRUE--------------------------------------------------------------------------------------------------------
source("model_measSIR.R")
measSIR@params


## ----pfilter-1,cache=TRUE-------------------------------------------------------------------------------------------------------------
library(pomp)
pf <- measSIR |> pfilter(Np=5000)
logLik(pf)


## ----pf-diagnostic-1, fig.height=1.5, fig.width=8-------------------------------------------------------------------------------------
measSIR |>
  simulate(nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")




## ----pf-diagnostic-2-2, echo=FALSE, fig.width=4, fig.height=3.5-----------------------------------------------------------------------
plot(pf)


## ----pf-diagnostic-3,cache=TRUE-------------------------------------------------------------------------------------------------------
plan(multisession)
foreach (
  i=1:10, .combine=c, .options.future=list(seed=652643293)
) %dofuture% {
    measSIR |> pfilter(Np=5000)
} -> pf
logLik(pf) -> ll
logmeanexp(ll,se=TRUE)






## ----like-slice-eval,include=FALSE,purl=TRUE------------------------------------------------------------------------------------------
## What is this 'bake' function?
## See https://kingaa.github.io/sbied/misc/bake.html
## for an explanation.
bake(file="like-slice.rds",{
  slice_design(
    center = coef(measSIR),
    Beta = rep(seq(from=5,to=30,length=40),each=3),
    Gamma = rep(seq(from=0.2,to=2,length=40),each=3)
  ) -> param_slice
  
  dim(param_slice)
  library(iterators)
  plan(multisession)
  foreach (theta=iter(param_slice,"row"), .combine=rbind, 
    .options.future=list(seed=108028909)) %dofuture% {
    measSIR |> pfilter(params=theta,Np=5000) -> pf
    theta$loglik <- logLik(pf)
    theta
  } -> lik_slice
}) -> lik_slice


## ----like-slice-plot,echo=FALSE-------------------------------------------------------------------------------------------------------
#| fig-width: 6
#| fig-height: 3
#| out-width: 5in
#| out-height: 2in
lik_slice |>
  pivot_longer(c(Beta,Gamma)) |>
  filter(name==slice) |>
  ggplot(aes(x=value,y=loglik,color=name)) +
  geom_point() +
  facet_grid(
    ~tolower(name),scales="free_x",
    labeller=label_parsed
  ) +
  guides(color="none") + labs(x="parameter value",color="")






## ----pfilter-grid1-eval,include=FALSE,purl=TRUE---------------------------------------------------------------------------------------
bake(file="pfilter-grid1.rds",{
  expand.grid(
    Beta = rep(seq(from=10,to=30,length=40), each=3),
    Gamma = rep(seq(from=0.4,to=1.5,length=40), each=3),
    Rho = 0.5, k=10, Eta=0.06, N=38000
  ) -> param_grid
  
  dim(param_grid)
  plan(multisession)
  foreach (theta=iter(param_grid,"row"), .combine=rbind,
    .options.future=list(seed=421776444)) %dofuture% {
    measSIR |> pfilter(params=theta,Np=5000) -> pf
    theta$loglik <- logLik(pf)
    theta
  } -> lik_grid
  lik_grid |> arrange(Beta,Gamma) -> lik_grid
})-> lik_grid


## ----pfilter-grid1-plot,echo=FALSE,purl=TRUE------------------------------------------------------------------------------------------
#| fig-width: 6
#| fig-height: 4
#| out-width: 4.8in
#| out-height: 3.2in
lik_grid |>
  group_by(Beta,Gamma) |>
  summarize(loglik=logmeanexp(loglik)) |>
  ungroup() |>
  mutate(loglik=ifelse(loglik>max(loglik)-25,loglik,NA)) |>
  ggplot(aes(x=Beta,y=Gamma,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient()+
  labs(x=expression(beta),y=expression(gamma))

