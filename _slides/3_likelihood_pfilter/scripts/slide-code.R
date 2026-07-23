## =============================================================================
##  Lecture A -- slide code: the likelihood, and the particle filter
##
##  The key runnable code BEHIND the slides -- not the polished figure scripts
##  in figures/, which exist for illustration only. Open
##  3_likelihood_pfilter.Rproj first, then run this top to bottom.
## =============================================================================

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)

## The model: the Consett 1948 measles SIR
source("scripts/model_measSIR.R")   # builds `measSIR`
coef(measSIR)


# Look before you compute: simulate the model against the data ------------

set.seed(1350254336)
measSIR |>
  simulate(nsim=20, format="data.frame", include.data=TRUE) |>
  ggplot(aes(week, reports, group=.id, color=(.id=="data"))) +
  geom_line() +
  scale_color_manual(values=c("FALSE"="grey70", "TRUE"="#E69F00"), guide="none") +
  labs(x="week", y="reported cases") +
  theme_bw()


# Compute the likelihood with the particle filter -------------------------

pf <- measSIR |> pfilter(Np=1000)
logLik(pf)


# The answer is RANDOM: run it again, get a different number ---------------

replicate(5, logLik(measSIR |> pfilter(Np=1000)))

## Average on the LIKELIHOOD scale (pfilter is unbiased for L, not for log L):
replicate(10, logLik(measSIR |> pfilter(Np=1000))) |> logmeanexp(se=TRUE)


# Run many filters in parallel, then combine ------------------------------

cl <- makePSOCKcluster(max(1, detectCores() - 1))
registerDoParallel(cl)
registerDoRNG(1234)

foreach(i=1:20, .combine=c, .packages="pomp", .export="measSIR") %dopar% {
  measSIR |> pfilter(Np=5000)
} -> pfs

pfs |> logLik() |> logmeanexp(se=TRUE)


# Diagnostics: do I have enough particles? --------------------------------

pf <- measSIR |> pfilter(Np=1000)
pf |> eff_sample_size() |> plot(type="l", xlab="week", ylab="effective sample size")
pf |> cond_logLik()      # per-week log-likelihood contributions


# A likelihood SLICE: fix everything but Eta, vary Eta --------------------

slice_design(center=coef(measSIR),
             Eta=rep(seq(from=0.01, to=0.2, length=30), each=3)) -> p_eta

registerDoRNG(4321)
foreach(theta=iter(p_eta,"row"), .combine=rbind, .packages="pomp") %dopar% {
  measSIR |> pfilter(params=theta, Np=2000) -> pf
  theta$loglik <- logLik(pf)
  theta
} -> eta_slice

eta_slice |>
  ggplot(aes(Eta, loglik)) +
  geom_point(color="#0072B2", alpha=.6) +
  labs(x=expression(eta), y="log likelihood") +
  theme_bw()


# Cleanup -----------------------------------------------------------------

stopCluster(cl)
