## ----prelims,echo=F,cache=F----------------------------------------------------------------------------------------------------------------------
library(foreach)
library(iterators)
library(doFuture)
library(tidyverse)
library(microbenchmark)
library(pomp)
plan(multisession)


## ----meas-data1----------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
read_csv("Measles_Consett_1948.csv") |>
  select(week,reports=cases) |>
  filter(week<=42) -> meas
meas |> as.data.frame() |> head(n=3)


## ----meas-data2, fig.width=5, fig.height=1.5, dpi=300--------------------------------------------------------------------------------------------
meas |>
  ggplot(aes(x=week,y=reports)) +
  geom_line() + geom_point()


## ----sir_pomp1-----------------------------------------------------------------------------------------------------------------------------------
sir_stoch <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-Gamma*dt));
  S -= dN_SI; I += dN_SI - dN_IR;
  R += dN_IR; H += dN_IR;"
)

sir_init <- Csnippet("
  S = nearbyint(Eta*N);
  I = 1;
  R = nearbyint((1-Eta)*N);
  H = 0;"
)


## ----sir_pomp2-----------------------------------------------------------------------------------------------------------------------------------
dmeas <- Csnippet("lik = dnbinom_mu(reports,k,Rho*H,give_log);")
rmeas <- Csnippet("reports = rnbinom_mu(k,Rho*H);")

meas |>
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_stoch,delta.t=1/7),
    rinit=sir_init, rmeasure=rmeas,
    dmeasure=dmeas, accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","Gamma","Eta","Rho","k","N"),
    params=c(Beta=15,Gamma=0.5,Rho=0.5,k=10,Eta=0.06,N=38000)
  ) -> measSIR


## ----fig.width=8, fig.height=2-------------------------------------------------------------------------------------------------------------------
measSIR |>
  simulate(nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")


## ----fig.width=5, fig.height=2-------------------------------------------------------------------------------------------------------------------
pf <- measSIR |> pfilter(Np=5000)
min(pf@eff.sample.size)


## ----fig.width=7, fig.height=4-------------------------------------------------------------------------------------------------------------------
#| out-height: 80%
plot(pf)


## ----pf2-----------------------------------------------------------------------------------------------------------------------------------------
foreach(i=1:10,.combine=c,
  .options.future=list(seed=TRUE)) %dofuture% {
  measSIR |> pfilter(Np=5000)
} -> pf
pf |> logLik() |> logmeanexp(se=TRUE) -> L_pf
L_pf


## ----init_csv,cache=FALSE------------------------------------------------------------------------------------------------------------------------
pf[[1]] |> coef() |> bind_rows() |>
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2]) |>
  write_csv("measles_params.csv")




## ----local_search_eval,echo=FALSE----------------------------------------------------------------------------------------------------------------
bake(file="local_search.rds",{
  foreach(i=1:20,.combine=c,
    .options.future=list(seed=482947940)) %dofuture% {
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


## ----local_search_plot,fig.width=6, fig.height=4, dpi=300----------------------------------------------------------------------------------------
#| out-height: 80%
mifs_local |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")


## ----cluster_setup,include=FALSE,purl=TRUE-------------------------------------------------------------------------------------------------------
if (file.exists("CLUSTER.R")) {
  source("CLUSTER.R")
}


## ----global_search1------------------------------------------------------------------------------------------------------------------------------
set.seed(2062379496) 
runif_design(
  lower=c(Beta=5,Rho=0.2,Eta=0),
  upper=c(Beta=80,Rho=0.9,Eta=1),
  nseq=400
) -> guesses
fixed_params <- c(N=38000, Gamma=2, k=10)
mf1 <- mifs_local[[1]]




## ----global_search_eval,include=FALSE------------------------------------------------------------------------------------------------------------
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


## ----cache=FALSE,include=FALSE-------------------------------------------------------------------------------------------------------------------
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


## ----mle_global----------------------------------------------------------------------------------------------------------------------------------
results |> filter(loglik==max(loglik))


## ----pairs_global1,fig.width=8, fig.height=6, dpi=300--------------------------------------------------------------------------------------------
#| out-height: 80%
read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-50) |>
  bind_rows(guesses) |>
  mutate(type=if_else(is.na(loglik),"guess","result")) |>
  arrange(type) -> all


## ----pairs_global11,fig.width=8, fig.height=6, dpi=300-------------------------------------------------------------------------------------------
#| out-height: 80%
pairs(~loglik+Beta+Eta+Rho, data=all, pch=16, cex=0.3,
  col=ifelse(all$type=="guess",grey(0.5),"red"))


## ----pairs_global2, fig.width=8, fig.height=6, dpi=300-------------------------------------------------------------------------------------------
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


## ----eta_profile1a-------------------------------------------------------------------------------------------------------------------------------
read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-20,loglik.se<2) |>
  sapply(range) -> box
box


## ----eta_profile1b, fig.width=4, fig.height=2, dpi=300-------------------------------------------------------------------------------------------
freeze(seed=1196696958,
  profile_design(
    Eta=seq(0.01,0.95,length=40),
    lower=box[1,c("Beta","Rho")],
    upper=box[2,c("Beta","Rho")],
    nprof=15, type="runif"
  )) -> guesses


## ----fig.width=6, fig.height=4, dpi=300----------------------------------------------------------------------------------------------------------
#| out-height: 60%
plot(guesses)




## ----eta_profile2_eval,include=FALSE-------------------------------------------------------------------------------------------------------------
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
      replicate(10, mf |> pfilter(Np=5000) |> logLik()) |>
        logmeanexp(se=TRUE) -> ll
      mf |> coef() |> bind_rows() |>
        bind_cols(loglik=ll[1],loglik.se=ll[2])
    } -> results
    attr(results,"ncpu") <- nbrOfWorkers()
    results
  }) -> results
t_eta <- attr(results,"system.time")
ncpu_eta <- attr(results,"ncpu")


## ----eta_profile_database, echo=FALSE, cache=FALSE-----------------------------------------------------------------------------------------------
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


## ----eta_profile_pairs, fig.width=8, fig.height=6, dpi=300---------------------------------------------------------------------------------------
#| out-height: 70%
read_csv("measles_params.csv") |> 
  filter(loglik>max(loglik)-10) -> all

pairs(~loglik+Beta+Eta+Rho,data=all,pch=16)


## ----eta_profile_plot2, fig.width=6, fig.height=5, dpi=300---------------------------------------------------------------------------------------
#| out-height: 80%
results |>
  filter(is.finite(loglik)) |>
  group_by(round(Eta,5)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  filter(loglik>max(loglik)-20) |>
  ggplot(aes(x=Eta,y=loglik))+
  geom_point() + xlab(expression(eta))


## ----eta_profile_plot3, fig.width=6, fig.height=4, dpi=300---------------------------------------------------------------------------------------
#| out-height: 80%
maxloglik <- max(results$loglik,na.rm=TRUE)
ci.cutoff <- maxloglik-0.5*qchisq(df=1,p=0.95)

results |>
  filter(is.finite(loglik)) |>
  group_by(round(Eta,5)) |>
  filter(rank(-loglik)<3) |> ungroup() |>
  ggplot(aes(x=Eta,y=loglik))+
  geom_point() + xlab(expression(eta)) +
  geom_smooth(method="loess",span=0.25)+
  geom_hline(color="red",yintercept=ci.cutoff)+
  lims(y=maxloglik-c(5,0))


## ----eta_ci--------------------------------------------------------------------------------------------------------------------------------------
results |>
  filter(is.finite(loglik)) |>
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) |> 
  summarize(min=min(Eta),max=max(Eta)) -> Eta_ci


## ----rho_profile1--------------------------------------------------------------------------------------------------------------------------------
read_csv("measles_params.csv") |>
  group_by(cut=round(Rho,2)) |>
  filter(rank(-loglik)<=10) |>
  ungroup() |>
  arrange(-loglik) |>
  select(-cut,-loglik,-loglik.se) -> guesses




## ----rho_profile_eval,include=FALSE--------------------------------------------------------------------------------------------------------------
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


## ----include=FALSE,cache=FALSE-------------------------------------------------------------------------------------------------------------------
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


## ----prof-rho1, echo=FALSE, fig.width=8, fig.height=6, dpi=300-----------------------------------------------------------------------------------
#| out-height: 80%
results |> filter(is.finite(loglik)) -> results

pairs(~loglik+Beta+Eta+Rho,data=results,pch=16)


## ----echo=FALSE, fig.width=6, fig.height=4, dpi=300----------------------------------------------------------------------------------------------
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


## ----rho_ci--------------------------------------------------------------------------------------------------------------------------------------
results |>
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) |>
  summarize(min=min(Rho),max=max(Rho)) -> rho_ci


## ------------------------------------------------------------------------------------------------------------------------------------------------
read_csv("measles_params.csv") |>
  filter(loglik == max(loglik)) |>
  select(-loglik, -loglik.se) -> best.params

measSIR |>
  simulate(
    params=unlist(best.params),
    nsim=1000, format="data.frame", include.data=TRUE
  ) -> sims


## ----viz-predictions, fig.height=4, fig.width=6, dpi=300-----------------------------------------------------------------------------------------
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

