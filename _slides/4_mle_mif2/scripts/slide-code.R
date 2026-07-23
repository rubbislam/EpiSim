## =============================================================================
##  Lecture 4 -- slide code: finding the MLE with iterated filtering (mif2)
##
##  Open 4_mle_mif2.Rproj first, then run this top to bottom.
## =============================================================================

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)
library(reshape2)   # for melt()

source("scripts/model_measSIR.R")   # builds `measSIR`

## Build the parameter transformation ONCE, outside any foreach loop
measSIR |>
  pomp(
    partrans=parameter_trans(log="Beta", logit=c("Rho","Eta")),
    paramnames=c("Beta","Rho","Eta")
  ) -> po


# One mif2 run: the four knobs --------------------------------------------

set.seed(42)
po |>
  mif2(
    Np=2000,                     # particles
    Nmif=50,                     # iterations
    cooling.fraction.50=0.5,     # how fast perturbations shrink
    rw.sd=rw_sd(                 # how big perturbations start
      Beta=0.02, Rho=0.02, Eta=ivp(0.02)
    )
  ) -> mf

plot(mf)

## Continue a run for more iterations (re-runs from the current best):
mf |> mif2(Nmif=50) -> mf


# Twenty searches at once, in parallel ------------------------------------

cl <- makePSOCKcluster(max(1, detectCores() - 1))
registerDoParallel(cl)

foreach(i=1:20, .combine=c, .packages="pomp", .export="po",
        .options.RNG=482947940) %dorng% {
  po |> mif2(Np=2000, Nmif=50, cooling.fraction.50=0.5,
             rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02)))
} -> mifs_local

## The convergence traces (the "Twenty searches at once" figure)
mifs_local |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration, y=value, group=L1)) +
  geom_line(alpha=.3, color="#0072B2") +
  facet_wrap(~name, scales="free_y") +
  labs(x="mif2 iteration", y=NULL) +
  theme_bw()


# The number mif2 prints is NOT your likelihood: always re-filter ---------

foreach(mf=mifs_local, .combine=rbind, .packages=c("pomp","dplyr"),
        .options.RNG=900242057) %dorng% {
  evals <- replicate(10, logLik(pfilter(mf, Np=5000)))
  ll <- logmeanexp(evals, se=TRUE)
  mf |> coef() |> bind_rows() |> bind_cols(loglik=ll[1], loglik.se=ll[2])
} -> results

results |> filter(loglik==max(loglik))    # your best fit


# Cleanup -----------------------------------------------------------------

stopCluster(cl)
