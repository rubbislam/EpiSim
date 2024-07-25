## ----prelims,echo=F,cache=F---------------------------------------------------------------
library(foreach)
library(iterators)
library(doFuture)
library(tidyverse)
library(pomp)
set.seed(1350254336)


## ----model-construct----------------------------------------------------------------------
source("model_measSIR.R")




## ----fixed_params-------------------------------------------------------------------------
fixed_params <- c(N=38000, Gamma=2, k=10)
coef(measSIR,names(fixed_params)) <- fixed_params
coef(measSIR)


## ----sanity, echo=FALSE, fig.width=6, fig.height=4, dpi=300-------------------------------
#| out-height: 70%
measSIR |>
  simulate(nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")


## ----parallel-setup,cache=FALSE-----------------------------------------------------------
library(foreach)
library(doFuture)
plan(multisession)




## ----pf,echo=FALSE------------------------------------------------------------------------
tic <- Sys.time()
foreach(i=1:10,.combine=c,
  .options.future=list(seed=TRUE)) %dofuture% {
  measSIR |> pfilter(Np=5000)
} -> pf
pf |> logLik() |> logmeanexp(se=TRUE) -> L_pf
L_pf
toc <- Sys.time()


## ----init_csv,cache=FALSE-----------------------------------------------------------------
pf[[1]] |> coef() |> bind_rows() |>
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2]) |>
  write_csv("measles_params.csv")




## ----local_search_eval,echo=FALSE---------------------------------------------------------
## What is this 'bake' function?
## See https://kingaa.github.io/sbied/pfilter/bake.html
## for an explanation.
bake(file="local_search.rds",{
  foreach(i=1:20,.combine=c,
    .options.future=list(seed=482947940)
  ) %dofuture% {
    measSIR |>
      mif2(
        Np=2000, Nmif=50, cooling.fraction.50=0.5,
        rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02)),
        partrans=parameter_trans(
          log="Beta",logit=c("Rho","Eta")
        ),
        paramnames=c("Beta","Rho","Eta")
      )
  } -> mifs_local
  attr(mifs_local,"ncpu") <- nbrOfWorkers()
  mifs_local
}) -> mifs_local
t_loc <- attr(mifs_local,"system.time")
ncpu_loc <- attr(mifs_local,"ncpu")


## ----local_search_plot,fig.width=6, fig.height=4, dpi=300---------------------------------
#| out-height: 80%
mifs_local |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")




## ----lik_local_eval,include=FALSE---------------------------------------------------------
bake(file="lik_local.rds",{
  foreach(mf=mifs_local,.combine=rbind,
    .options.future=list(seed=900242057)
  ) %dofuture% {
    evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
    ll <- logmeanexp(evals,se=TRUE)
    mf |> coef() |> bind_rows() |>
      bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results
  attr(results,"ncpu") <- nbrOfWorkers()
  results
}) -> results
t_local <- attr(results,"system.time")
ncpu_local <- attr(results,"ncpu")


## ----mle_local----------------------------------------------------------------------------
results |> filter(loglik==max(loglik))


## ----pairs_local,fig.width=8, fig.height=6, dpi=300---------------------------------------
#| out-height: 80%
pairs(~loglik+Beta+Eta+Rho,data=results,pch=16)


## ----local_database,cache=FALSE-----------------------------------------------------------
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


## ----cluster_setup,include=FALSE,purl=TRUE------------------------------------------------
if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
}


## ----global_search1-----------------------------------------------------------------------
set.seed(2062379496)

runif_design(
  lower=c(Beta=5,Rho=0.2,Eta=0),
  upper=c(Beta=80,Rho=0.9,Eta=1),
  nseq=400
) -> guesses

mf1 <- mifs_local[[1]]




## ----global_search_eval,include=FALSE-----------------------------------------------------
bake(file="global_search.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=1270401374)
    ) %dofuture% {
      mf1 |>
        mif2(params=c(guess,fixed_params)) |>
        mif2(Nmif=100) -> mf
      replicate(
        10,
        mf |> pfilter(Np=5000) |> logLik()
      ) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) |>
  filter(is.finite(loglik)) -> results
t_global <- attr(results,"system.time")
ncpu_global <- attr(results,"ncpu")


## ----cache=FALSE,include=FALSE------------------------------------------------------------
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


## ----mle_global---------------------------------------------------------------------------
results |> filter(loglik==max(loglik))


## ----pairs_global1,fig.width=8, fig.height=6, dpi=300-------------------------------------
#| out-height: 80%
read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-50) |>
  bind_rows(guesses) |>
  mutate(type=if_else(is.na(loglik),"guess","result")) |>
  arrange(type) -> all

pairs(~loglik+Beta+Eta+Rho, data=all, pch=16, cex=0.3,
  col=ifelse(all$type=="guess",grey(0.5),"red"))


## ----pairs_global2, fig.width=8, fig.height=6, dpi=300------------------------------------
#| out-height: 80%
all |>
  filter(type=="result") |>
  filter(loglik>max(loglik)-10) |>
  ggplot(aes(x=Eta,y=loglik))+
  geom_point()+
  labs(
    x=expression(Eta),
    title="poor man's profile likelihood"
  )


## ----eta_profile1a------------------------------------------------------------------------
read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-20,loglik.se<2) |>
  sapply(range) -> box
box


## ----eta_profile1b, fig.width=6, fig.height=4, dpi=300------------------------------------
#| out-height: 80%
freeze(seed=1196696958,
  profile_design(
    Eta=seq(0.01,0.95,length=40),
    lower=box[1,c("Beta","Rho")],
    upper=box[2,c("Beta","Rho")],
    nprof=15, type="runif"
  )) -> guesses
plot(guesses)




## ----eta_profile2_eval,include=FALSE------------------------------------------------------
mf1 <- mifs_local[[1]]
bake(file="eta_profile.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=830007657)
    ) %dofuture% {
      mf1 |>
        mif2(params=c(guess,fixed_params),
          rw.sd=rw_sd(Beta=0.02,Rho=0.02)) |>
        mif2(Nmif=100,cooling.fraction.50=0.3) -> mf
      replicate(
        10,
        mf |> pfilter(Np=5000) |> logLik()) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) -> results
t_eta <- attr(results,"system.time")
ncpu_eta <- attr(results,"ncpu")


## ----eta_profile_database,cache=FALSE-----------------------------------------------------
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


## ----eta_profile_pairs, fig.width=8, fig.height=6, dpi=300--------------------------------
#| out-height: 80%
read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-10) -> all

pairs(~loglik+Beta+Eta+Rho,data=all,pch=16)


## ----eta_profile_plot1, fig.width=6, fig.height=4, dpi=300--------------------------------
#| out-height: 80%
results |>
  ggplot(aes(x=Eta,y=loglik))+
  geom_point() + xlab(expression(eta))


## ----eta_profile_plot2, fig.width=6, fig.height=4, dpi=300--------------------------------
#| out-height: 80%
results |>
  filter(is.finite(loglik)) |>
  group_by(round(Eta,5)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  filter(loglik>max(loglik)-20) |>
  ggplot(aes(x=Eta,y=loglik))+
  geom_point() + xlab(expression(eta))


## ----eta_profile_plot3, fig.width=6, fig.height=4, dpi=300--------------------------------
#| out-height: 80%
maxloglik <- max(results$loglik,na.rm=TRUE)
ci.cutoff <- maxloglik-0.5*qchisq(df=1,p=0.95)

results |>
  filter(is.finite(loglik)) |>
  group_by(round(Eta,5)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  ggplot(aes(x=Eta,y=loglik))+
  geom_point() + xlab(expression(eta)) +
  geom_smooth(method="loess",span=0.25)+
  geom_hline(color="red",yintercept=ci.cutoff)+
  lims(y=maxloglik-c(5,0))


## ----eta_profile_eta_by_rho, fig.width=6, fig.height=4, dpi=300---------------------------
#| out-height: 80%
results |>
  filter(is.finite(loglik)) |>
  group_by(round(Eta,5)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  mutate(in_ci=loglik>max(loglik)-1.92) |>
  ggplot(aes(x=Eta,y=Rho,color=in_ci))+
  geom_point()+
  labs(
    color="inside 95% CI?",
    x=expression(eta),
    y=expression(rho),
    title="profile trace"
  )


## ----include=FALSE------------------------------------------------------------------------
results |>
  filter(is.finite(loglik)) |>
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) |>
  summarize(min=min(Rho),max=max(Rho)) -> rho_ci


## ----rho_profile1-------------------------------------------------------------------------
read_csv("measles_params.csv") |>
  group_by(cut=round(Rho,2)) |>
  filter(rank(-loglik)<=10) |>
  ungroup() |>
  arrange(-loglik) |>
  select(-cut,-loglik,-loglik.se) -> guesses




## ----rho_profile_eval,include=FALSE-------------------------------------------------------
mf1 <- mifs_local[[1]]
bake(file="rho_profile.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=2105684752)
    ) %dofuture% {
      mf1 |>
        mif2(params=guess,
          rw.sd=rw_sd(Beta=0.02,Eta=ivp(0.02))) |>
        mif2(Nmif=100,cooling.fraction.50=0.3) |>
        mif2() -> mf
      replicate(
        10,
        mf |> pfilter(Np=5000) |> logLik()) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) -> results
t_rho <- attr(results,"system.time")
ncpu_rho <- attr(results,"ncpu")


## ----include=FALSE,cache=FALSE------------------------------------------------------------
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


## ----prof-rho1, echo=FALSE, fig.width=8, fig.height=6, dpi=300----------------------------
#| out-height: 80%
results |> filter(is.finite(loglik)) -> results

pairs(~loglik+Beta+Eta+Rho,data=results,pch=16)


## ----echo=FALSE, fig.width=6, fig.height=4, dpi=300---------------------------------------
#| out-height: 80%
results |>
  filter(loglik>max(loglik)-10,loglik.se<1) |>
  group_by(round(Rho,2)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  ggplot(aes(x=Rho,y=loglik))+
  geom_point() + xlab(expression(rho)) +
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )


## ----rho_ci-------------------------------------------------------------------------------
results |>
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) |>
  summarize(min=min(Rho),max=max(Rho)) -> rho_ci


## ----exp_global_search1-------------------------------------------------------------------
freeze(seed=55266255,
  runif_design(
    lower=c(Beta=5,Gamma=0.2,Eta=0),
    upper=c(Beta=80,Gamma=5,Eta=0.99),
    nseq=1000
  )) |>
  mutate(
    Rho=0.6, k=10, N=38000
  ) -> guesses






## ----exp_global_search_eval,include=FALSE-------------------------------------------------
bake(file="global_search2.rds",
  dependson=guesses,{
    foreach(
      guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=610408798)
    ) %dofuture% {
      measSIR |>
        mif2(params=guess, Np=2000, Nmif=100,
          cooling.fraction.50=0.5,
          partrans=parameter_trans(
            log=c("Beta","Gamma"),
            logit="Eta"), paramnames=c("Beta","Gamma","Eta"),
          rw.sd=rw_sd(Beta=0.02,Gamma=0.02,Eta=ivp(0.02))) -> mf
      mf |>
        mif2(
          Nmif=100,rw.sd=rw_sd(Beta=0.01,Gamma=0.01,Eta=ivp(0.01))
        ) |>
        mif2(
          Nmif=100,
          rw.sd=rw_sd(Beta=0.005,Gamma=0.005,Eta=ivp(0.005))
        ) -> mf
      replicate(
        10,
        mf |> pfilter(Np=5000) |> logLik()
      ) |> logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) |>
  filter(is.finite(loglik)) -> results
t_expglob <- attr(results,"system.time")
ncpu_expglob <- attr(results,"ncpu")


## ----include=FALSE,cache=FALSE------------------------------------------------------------
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


## ----another-search, fig.width=8, fig.height=6, dpi=300-----------------------------------
#| out-height: 80%
read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-20) -> all

pairs(~loglik+Rho+Gamma+Beta+Eta,data=all,pch=16,cex=0.3,
  col=if_else(round(all$Rho,3)==0.6,1,4))


## ----poorman-gamma2, fig.width=8, fig.height=6, dpi=300-----------------------------------
#| out-height: 80%
results |>
  filter(loglik>max(loglik)-20,loglik.se<1) |>
  ggplot(aes(x=Gamma,y=loglik))+
  geom_point()+
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )


## ----gamma_profile1a----------------------------------------------------------------------
read_csv("measles_params.csv") |>
  filter(
    loglik>max(loglik)-20,
    loglik.se<2,
    abs(Rho-0.6)<0.01
  ) |>
  sapply(range) -> box


## ----gamma_profile1b----------------------------------------------------------------------
freeze(seed=610408798,
  profile_design(
    Gamma=seq(0.2,2,by=0.1),
    lower=box[1,c("Beta","Eta")],
    upper=box[2,c("Beta","Eta")],
    nprof=100, type="runif"
  )) |>
  mutate(
    N=38000,
    Rho=0.6,
    k=10
  ) -> guesses




## ----gamma_profile_eval,include=FALSE-----------------------------------------------------
bake(file="gamma_profile1.rds",
  dependson=guesses,{
    foreach(guess=iter(guesses,"row"), .combine=rbind,
      .options.future=list(seed=610408798)
    ) %dofuture% {
      measSIR |>
        mif2(params=guess, Np=2000, Nmif=100,
          partrans=parameter_trans(log="Beta",logit="Eta"),
          paramnames=c("Beta","Eta"), cooling.fraction.50=0.5,
          rw.sd=rw_sd(Beta=0.02,Eta=ivp(0.02))
        ) |>
        mif2(Nmif=100) |>
        mif2(Nmif=100,rw.sd=rw_sd(Beta=0.01,Eta=ivp(0.01))) |>
        mif2(Nmif=100,rw.sd=rw_sd(Beta=0.005,Eta=ivp(0.005))) -> mf
      replicate(10,mf |> pfilter(Np=5000) |> logLik()) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) |>
  filter(is.finite(loglik)) -> results
t_muIR <- attr(results,"system.time")
ncpu_muIR <- attr(results,"ncpu")


## ----include=FALSE,cache=FALSE------------------------------------------------------------
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


## ----gamma_profile_plot, fig.height=4, fig.width=6, dpi=300-------------------------------
#| out-height: 80%
results |>
  group_by(round(Gamma,2)) |>
  filter(rank(-loglik)<=1) |>
  ungroup() |>
  ggplot(aes(x=Gamma,y=loglik))+
  geom_point() + xlab(expression(gamma)) +
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )


## -----------------------------------------------------------------------------------------
read_csv("measles_params.csv") |>
  filter(loglik == max(loglik)) |>
  select(-loglik, -loglik.se) -> best.params

measSIR |>
  simulate(
    params=unlist(best.params),
    nsim=1000, format="data.frame", include.data=TRUE
  ) -> sims


## ----viz-predictions, fig.height=4, fig.width=6, dpi=300----------------------------------
#| out-height: 80%
sims |>
  mutate(data=.id=="data") |>
  group_by(week,data) |>
  reframe(
    p=c(0.025,0.5,0.975),
    value=wquant(reports,probs=p),
    name=c("lo","med","up")
  ) |>
  select(-p) |> pivot_wider() |> ungroup() |>
  ggplot(aes(x=week,y=med,color=data,fill=data,ymin=lo,ymax=up))+
  geom_line()+ geom_ribbon(alpha=0.2,color=NA) +
  labs(y="reports")+
  theme_bw() + guides(color="none",fill="none")

