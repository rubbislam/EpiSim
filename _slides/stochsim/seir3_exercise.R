params <-
list(prefix = "seir3_exercise")

## ----prelims,include=FALSE,purl=TRUE,cache=FALSE-------------------------------------------------------------------
library(tidyverse)
library(pomp)
theme_set(theme_bw())
options(stringsAsFactors=FALSE)
set.seed(1221234211)




## ----seir3_model---------------------------------------------------------------------------------------------------
library(tidyverse)
library(pomp)

bsflu |>
  select(day,B) |>
  pomp(times="day",t0=-6,
    rprocess=euler(Csnippet("
      double N = 763;
      double t1 = rbinom(S,1-exp(-Beta*I/N*dt));
      double t2 = rbinom(E,1-exp(-Mu_E*dt));
      double t3 = rbinom(I,1-exp(-Mu_I*dt));
      double t4 = rbinom(R1,1-exp(-Mu_R1*dt));
      double t5 = rbinom(R2,1-exp(-Mu_R2*dt));
      S  -= t1;
      E  += t1 - t2;
      I  += t2 - t3;
      R1 += t3 - t4;
      R2 += t4 - t5;"),
      delta.t=1/5),
    rinit=Csnippet("
      S  = 762;
      E  = 0;
      I  = 1;
      R1 = 0;
      R2 = 0;"),
    dmeasure=Csnippet("
      lik = dpois(B,Rho*R1+1e-6,give_log);"),
    rmeasure=Csnippet("
      B = rpois(Rho*R1+1e-6);"),
    statenames=c("S","E","I","R1","R2"),
    paramnames=c("Beta","Mu_E","Mu_I","Mu_R1","Mu_R2","Rho")
  ) -> flu


## ----fixedparams---------------------------------------------------------------------------------------------------
with(bsflu,c(Mu_R1=1/(sum(B)/512),Mu_R2=1/(sum(C)/512)))


## ----simulations---------------------------------------------------------------------------------------------------
coef(flu) <- c(Beta=5,Mu_E=0.5,Mu_I=1,Mu_R1=0.33,Mu_R2=0.55,Rho=0.95)

flu |>
  simulate(nsim=20,format="data.frame",include.data=TRUE) -> simdat

simdat |>
  select(day,B,.id) |>
  ggplot(aes(x=day,y=B,color=(.id=="data"),size=(.id=="data"),group=.id))+
  geom_line()+
  scale_color_manual(values=c(`TRUE`="red",`FALSE`="black"))+
  scale_size_manual(values=c(`TRUE`=1,`FALSE`=0.5))+
  guides(color="none",size="none")

