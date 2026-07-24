params <-
list(prefix = "expense")

## ----cache=FALSE--------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(pomp)
library(iterators)
library(doFuture)


## ----model-construct----------------------------------------------------------------------------------------------------------------------------
source("model_measSIR.R")




## ----comp_eval,include=FALSE--------------------------------------------------------------------------------------------------------------------
stew(file="expense1.rda",{
  expand_grid(
    Np=ceiling(10^seq(1,5,by=0.2))
  ) -> design
  
  plan(multisession)
  foreach (
    expt=iter(design,"row"),
    .combine=bind_rows,
    .options.future=list(seed=TRUE)
  ) %dofuture% {
    system.time(measSIR |> pfilter(Np=expt$Np))[3] -> expt$time
    expt
  } -> resultA
})


## -----------------------------------------------------------------------------------------------------------------------------------------------
resultA |>
  ggplot(aes(x=Np,y=time))+
  geom_point()+
  geom_smooth(method="lm",formula=y~x)+
  expand_limits(x=0,y=0)

lm(time~Np,data=resultA) -> fit
summary(fit)



## ----comp2_eval,include=FALSE-------------------------------------------------------------------------------------------------------------------
stew(file="expense2.rda",{
  measSIR |>
    window(end=21) -> shortMeasSIR
  
  plan(multisession)
  foreach (
    expt=iter(design,"row"),
    .combine=bind_rows,
    .options.future=list(seed=TRUE)
  ) %dofuture% {
    system.time(shortMeasSIR |> pfilter(Np=expt$Np))[3] -> expt$time
    expt
  } -> resultB
  
  bind_rows(
    long=resultA,
    short=resultB,
    .id="length"
  ) |>
    mutate(
      n=case_when(
        length=="short"~length(time(shortMeasSIR)),
        length=="long"~length(time(measSIR))
      )
    ) -> result
})


## -----------------------------------------------------------------------------------------------------------------------------------------------
result |>
  ggplot(aes(x=Np,y=time,group=n,color=factor(n)))+
  geom_point()+
  labs(color="n")+
  geom_smooth(method="lm",formula=y~x)+
  expand_limits(x=0,y=0)

lm(time~n*Np,data=result) -> fit
summary(fit)

