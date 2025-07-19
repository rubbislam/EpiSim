
# Setup -------------------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)
library(iterators)
options(stringsAsFactors=FALSE)
set.seed(1350254336)


# Particle filtering in pomp ----------------------------------------------

source("model_measSIR.R")
measSIR@params

# NOTE: IF The code above didn't work, try the following to build the model 
# from scratch:

# sir_stoch <- Csnippet("
#   double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
#   double dN_IR = rbinom(I,1-exp(-Gamma*dt));
#   S -= dN_SI;
#   I += dN_SI - dN_IR;
#   R += dN_IR;
#   H += dN_IR;"
# )
# 
# sir_init <- Csnippet("
#   S = nearbyint(Eta*N);
#   I = 1;
#   R = nearbyint((1-Eta)*N);
#   H = 0;"
# )
# 
# dmeas <- Csnippet("
#   lik = dnbinom_mu(reports,k,Rho*H,give_log);"
# )
# 
# rmeas <- Csnippet("
#   reports = rnbinom_mu(k,Rho*H);"
# )
# 
# emeas <- Csnippet("
#   E_reports = Rho*H;"
# )
# 
# read_csv(".../data/Measles_Consett_1948.csv") |>
#   select(week,reports=cases) |>
#   filter(week<=42) |>
#   pomp(
#     times="week",t0=0,
#     rprocess=euler(sir_stoch,delta.t=1/7),
#     rinit=sir_init,
#     rmeasure=rmeas,
#     dmeasure=dmeas,
#     emeasure=emeas,
#     accumvars="H",
#     statenames=c("S","I","R","H"),
#     paramnames=c("Beta","Gamma","Eta","Rho","k","N"),
#     params=c(Beta=15,Gamma=0.5,Rho=0.5,k=10,Eta=0.06,N=38000)
#   ) -> measSIR
# 
# invisible(gc())

# Particle Filtering in POMP II -------------------------------------------

measSIR |>
  simulate(nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")


# Particle Filtering in POMP III ------------------------------------------

pf <- measSIR |> pfilter(Np=5000)
logLik(pf)

# Particle Filtering in POMP IV -------------------------------------------

plot(pf)

# Particle Filtering in POMP VI -------------------------------------------

cores <- parallel::detectCores()
registerDoParallel(cores-1)  # MacOS / Linux
registerDoRNG(seed = 123456)

foreach (
  i=1:10, .combine=c
) %dopar% {
  measSIR |> pfilter(Np=5000)
} -> pf

logLik(pf) -> ll
logmeanexp(ll,se=TRUE)

# Slicing Measles SIR Likelihood ------------------------------------------

slice_design(
  center = coef(measSIR),
  Beta = rep(seq(from=5,to=30,length=40),each=3),
  Gamma = rep(seq(from=0.2,to=2,length=40),each=3)
) -> param_slice

dim(param_slice)

# Slicing Measles SIR Likelihood -----------------------------------------

library(iterators)
# Doesn't need to happen again, but included as reminder. 
registerDoParallel(cores-1)  # For MacOS / Linux
registerDoRNG(seed = 654321)

foreach (theta=iter(param_slice,"row"), .combine=rbind) %dopar% {
  measSIR |> pfilter(params=theta,Np=5000) -> pf
  theta$loglik <- logLik(pf)
  theta
} -> lik_slice

# Slice Plot --------------------------------------------------------------

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


# Two Dim Slice  ----------------------------------------------------------

# WARNING: Consider making a smaller grid, otherwise it will take a long 
# time on your laptop. E.g., set length = 10 for Beta and Gamma. 
expand.grid(
  Beta = rep(seq(from=10,to=30,length=40), each=3),
  Gamma = rep(seq(from=0.4,to=1.5,length=40), each=3),
  Rho = 0.5, k=10, Eta=0.06, N=38000
) -> param_grid

dim(param_grid)


# Two Dim Slice (WARNING) -------------------------------------------------
# ::::::WARNING:::::::
# This code requires a large number of computations, and may take a long 
# time on your laptop (1+ hours)
# 
# You can try reducing the time by making a smaller grid above. 


# Doesn't need to happen again, but included as reminder. 
registerDoParallel(cores-1)  # For MacOS / Linux
registerDoRNG(seed = 111111)

foreach (theta=iter(param_grid,"row"), .combine=rbind,
         .options.future=list(seed=421776444)) %dopar% {
           measSIR |> pfilter(params=theta,Np=5000) -> pf
           theta$loglik <- logLik(pf)
           theta
         } -> lik_grid


# Two-dim slice figure ----------------------------------------------------

lik_grid |>
  group_by(Beta,Gamma) |>
  summarize(loglik=logmeanexp(loglik)) |>
  ungroup() |>
  mutate(loglik=ifelse(loglik>max(loglik)-25,loglik,NA)) |>
  ggplot(aes(x=Beta,y=Gamma,z=loglik,fill=loglik))+
  geom_tile(color=NA)+
  scale_fill_gradient()+
  labs(x=expression(beta),y=expression(gamma))
