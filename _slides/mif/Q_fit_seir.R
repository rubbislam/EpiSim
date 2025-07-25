params <-
list(prefix = "Q_fit_seir")

library(tidyverse)
library(pomp)
library(iterators)
library(doFuture)
if (.Platform$OS.type=="windows") 
  options(pomp_cdir="./tmp")

source("model_measSIR.R")

seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-Sigma*dt));
  double dN_IR = rbinom(I,1-exp(-Gamma*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

seir_rinit <- Csnippet("
  S = nearbyint(Eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-Eta)*N);
  H = 0;
")

measSIR |>
  pomp(
    rprocess=euler(seir_step,delta.t=1/7),
    rinit=seir_rinit,
    partrans=parameter_trans(
      log=c("Beta","Sigma"),
      logit=c("Eta","Rho")
    ),
    paramnames=c("N","Beta","Sigma","Gamma","Eta","k","Rho"),
    statenames=c("S","E","I","R","H")
  ) -> measSEIR

read_csv("measles_params.csv") |>
  filter(abs(Gamma-2)<0.001) |>
  filter(loglik==max(loglik)) |>
  select(-loglik,-loglik.se) -> coef(measSEIR)

coef(measSEIR,"Sigma") <- 0.8
fixed_params <- coef(measSEIR,c("N","Gamma","k"))
coef(measSEIR)

set.seed(1014406)
measSEIR |>
  simulate(nsim=20,format="data.frame",include=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=(.id=="data")))+
  geom_line()+
  guides(color="none")+
  theme_bw()

measSEIR |>
  pfilter(Np=1000) -> pf1
logLik(pf1)
plot(pf1)

plan(multisession)

ncpu <- nbrOfWorkers()
bake(file="Q_fit_seir_local_mifs.rds",{
  foreach(i=seq_len(ncpu),.combine=c,
    .options.future=list(seed=482947940)
  ) %dofuture% {
    measSEIR |>
      mif2(
        Np=1000, Nmif=50,
        cooling.fraction.50=0.5,
        rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02), Sigma=0.02)
      )
  }
}) -> local_mifs

local_mifs |>
  traces(pars=c("loglik","Beta","Sigma","Rho","Eta")) |>
  melt() |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")

bake(file="Q_fit_seir_local_logliks.rds",{
  foreach(mf=local_mifs,.combine=rbind,
    .options.future=list(seed=901242057)
  ) %dofuture% {
    evals <- replicate(10, logLik(pfilter(mf,Np=2000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf |> coef() |> bind_rows() |>
      bind_cols(loglik=ll[1],loglik.se=ll[2],Np=2000,nfilt=10)
  }
}) -> local_logliks

local_logliks$loglik
local_logliks$loglik.se

etime <- attr(local_mifs,"system.time")+attr(local_logliks,"system.time")
mif_work <- sum(sapply(local_mifs,slot,"Nmif")*apply(sapply(local_mifs,slot,"Np"),2,mean))
pfilter_work <- with(local_logliks,sum(Np*nfilt))
efactor <- unname(etime[3]*ncpu/(mif_work+pfilter_work)*1000)
efactor

unit_cost <- (100*1000+10*2000)/1000*efactor
budget <- 60*ncpu
budget/unit_cost

freeze(
  runif_design(
    lower=c(Beta=5,Rho=0.2,Eta=0,Sigma=1/3),
    upper=c(Beta=80,Rho=0.9,Eta=1,Sigma=3),
    nseq=45
  ),
  seed=2062379496
)-> guesses

bake(file="Q_fit_seir_global1.rds",{
  foreach(guess=iter(guesses,"row"), .combine=rbind,
    .options.future=list(seed=1270401374)
  ) %dofuture% {
    measSEIR |>
      mif2(
        Nmif=100, Np=1000,
        cooling.fraction.50=0.5,
        rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02), Sigma=0.02),
        params=c(unlist(guess),fixed_params)
      ) -> mf
    replicate(
      10,
      mf |> pfilter(Np=2000) |> logLik()
    ) |>
      logmeanexp(se=TRUE) -> ll
    mf |> coef() |> bind_rows() |>
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> global1
  attr(global1,"ncpu") <- nbrOfWorkers()
  global1
}) -> global1

pairs(
  ~loglik+Beta+Eta+Rho+Sigma,
  data=global1,
  pch=16
)

unit_cost <- (200*1000+10*2000)/1000*efactor
budget <- 5*60*ncpu
budget/unit_cost

freeze(
  runif_design(
    lower=c(Beta=5,Rho=0.2,Eta=0,Sigma=1/3),
    upper=c(Beta=80,Rho=0.9,Eta=1,Sigma=3),
    nseq=150
  ),
  seed=2062379496
)-> guesses

bake(file="Q_fit_seir_global2.rds",{
  foreach(guess=iter(guesses,"row"), .combine=rbind,
    .options.future=list(seed=1270401374)
  ) %dofuture% {
    measSEIR |>
      mif2(
        Nmif=100, Np=1000,
        cooling.fraction.50=0.5,
        rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02), Sigma=0.02),
        params=c(unlist(guess),fixed_params)
      ) |>
      continue(
        cooling.fraction=0.1
      ) -> mf
    replicate(
      10,
      mf |> pfilter(Np=2000) |> logLik()
    ) |>
      logmeanexp(se=TRUE) -> ll
    mf |> coef() |> bind_rows() |>
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> global2
  attr(global2,"ncpu") <- nbrOfWorkers()
  global2
}) -> global2

pairs(
  ~loglik+Beta+Eta+Rho+Sigma,
  data=global2,
  pch=16
)

unit_cost <- (200*1000+10*2000)/1000*efactor
budget <- 30*60*ncpu
budget/unit_cost



freeze(
  runif_design(
    lower=c(Beta=5,Rho=0.2,Eta=0,Sigma=1/3),
    upper=c(Beta=80,Rho=0.9,Eta=1,Sigma=3),
    nseq=1000
  ),
  seed=2062379496
)-> guesses

## ## Uncomment this code chunk to perform the run-level 3 calculation!
## foreach(guess=iter(guesses,"row"), .combine=rbind,
##   .options.future=list(seed=1270401374)
## ) %dofuture% {
##   measSEIR |>
##     mif2(
##       Nmif=100, Np=1000,
##       cooling.fraction.50=0.5,
##       rw.sd=rw_sd(Beta=0.02, rho=0.02, eta=ivp(0.02), mu_EI=0.02),
##       params=c(unlist(guess),fixed_params)
##     ) |>
##     continue(
##       cooling.fraction=0.1
##     ) -> mf
##   replicate(
##     10,
##     mf |> pfilter(Np=2000) |> logLik()
##   ) |>
##     logmeanexp(se=TRUE) -> ll
##   mf |> coef() |> bind_rows() |>
##     bind_cols(loglik=ll[1],loglik.se=ll[2])
## } -> global3



global3 |>
  count(
    finite=is.finite(loglik),
    precise=loglik.se>0.2,
    high=loglik>max(loglik,na.rm=TRUE)-10
  ) |>
  as.data.frame()

global3 |>
  filter(
    is.finite(loglik),
    loglik>max(loglik,na.rm=TRUE)-10
  ) -> global3

pairs(
  ~loglik+Beta+Eta+Rho+Sigma,
  data=global3,
  pch=16,cex=0.5
)

bind_rows(
  `1`=global1,
  `2`=global2,
  `3`=global3,
  .id="run-level"
) -> all

all |>
  filter(loglik>max(loglik)-20) |>
  ggplot(aes(x=Sigma,y=loglik,color=`run-level`))+
  geom_point(alpha=0.3,size=1)+
  guides(color="none")+
  labs(x=expression(sigma),y=expression(log~L)) -> pl1

all |>
  ggplot(aes(x=Sigma,y=loglik,color=`run-level`))+
  geom_point(alpha=0.3,size=1)+
  labs(x=expression(sigma),y=expression(log~L)) -> pl2

pl2+annotation_custom(grob=ggplotGrob(pl1),xmin=5,xmax=18,ymin=-430,ymax=-200)

bind_rows(
  global1,
  global2,
  global3
) |>
  filter(
    is.finite(loglik),
    loglik.se < 0.2
  ) |>
  filter(loglik==max(loglik)) -> mle_seir

read_csv("measles_params.csv") |>
  filter(abs(Gamma-2)<0.001) |>
  filter(loglik==max(loglik)) -> mle_sir


all |>
  mutate(
    bin=cut(Rho,breaks=50,include.lowest=TRUE)
  ) |>
  group_by(bin) |>
  filter(rank(-loglik)<=1) |>
  ungroup() |>
  filter(loglik>max(loglik)-10) -> poorman_prof

poorman_prof |>
  select(Sigma,loglik,Rho,Eta,Beta) |>
  pivot_longer(c(loglik,Eta,Sigma,Beta)) |>
  mutate(
    name=factor(name,levels=c("loglik","Eta","Sigma","Beta"))
  ) |>
  ggplot(aes(x=Rho,y=value))+
  geom_point()+
  labs(x=expression(rho),y="")+
  facet_wrap(~name,scales="free_y",ncol=1,strip.position="left")+
  theme(
    strip.placement="outside",
    strip.background=element_rect(fill=NA,color=NA)
  )

poorman_prof |>
  select(-loglik,-loglik.se) -> guesses

unit_cost <- (100*2000+10*2000)/1000*efactor
nrow(guesses)*unit_cost/ncpu/60

## ## Uncomment this code chunk to perform a profile computation
## foreach(guess=iter(guesses,"row"), .combine=rbind,
##   .options.future=list(seed=1270401374)
## ) %dofuture% {
##   measSEIR |>
##     mif2(
##       Nmif=100, Np=2000,
##       cooling.fraction.50=0.1,
##       rw.sd=rw_sd(Beta=0.02, eta=ivp(0.02), mu_EI=0.05),
##       params=c(unlist(guess),fixed_params)
##     ) -> mf
##   replicate(
##     10,
##     mf |> pfilter(Np=2000) |> logLik()
##   ) |>
##     logmeanexp(se=TRUE) -> ll
##   mf |> coef() |> bind_rows() |>
##     bind_cols(loglik=ll[1],loglik.se=ll[2])
## } -> profile_rho



profile_rho |>
  select(Sigma,loglik,Rho,Eta,Beta) |>
  pivot_longer(c(loglik,Eta,Sigma,Beta)) |>
  mutate(
    name=factor(name,levels=c("loglik","Eta","Sigma","Beta"))
  ) |>
  ggplot(aes(x=Rho,y=value))+
  geom_point()+
  labs(x=expression(rho),y="")+
  facet_wrap(~name,scales="free_y",ncol=1,strip.position="left")+
  theme(
    strip.placement="outside",
    strip.background=element_rect(fill=NA,color=NA)
  )
