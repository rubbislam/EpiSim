## =============================================================================
##  Lecture A: The particle filter -- SOLUTIONS
##
##  Try scripts/EXERCISES.R first. Open 3_likelihood_pfilter.Rproj before running this.
##  Works on macOS, Linux and Windows.
## =============================================================================


# Setup -------------------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)

source("scripts/model_measSIR.R")

cl <- makePSOCKcluster(max(1, detectCores() - 1))
registerDoParallel(cl)


# Exercise 1: does the model look anything like the data? ------------------

set.seed(1350254336)
measSIR |>
  simulate(nsim=20, format="data.frame", include.data=TRUE) |>
  ggplot(aes(x=week, y=reports, group=.id, color=(.id=="data"))) +
  geom_line() +
  scale_color_manual(values=c("FALSE"="grey70", "TRUE"="#E69F00")) +
  guides(color="none") +
  labs(x="week", y="reports", title="20 simulations (grey) vs the Consett data (amber)")

## The simulations are in the right ballpark but the epidemic peaks too early
## and too low. Our starting parameters are a guess, not a fit. That is the
## whole problem: we need a way to SCORE parameters so we can improve them.


# Exercise 2: your first likelihood ----------------------------------------

set.seed(1)
pf <- measSIR |> pfilter(Np=1000)
logLik(pf)
## about -136

## A single log-likelihood is meaningless in isolation. It becomes meaningful
## as a DIFFERENCE: this parameter set vs that one.


# Exercise 3: run it again ------------------------------------------------

replicate(5, logLik(measSIR |> pfilter(Np=1000)))
## Five different numbers. The pfilter is a Monte Carlo algorithm.

replicate(5, logLik(measSIR |> pfilter(Np=100)))    # jumps around a LOT, and is biased LOW
replicate(5, logLik(measSIR |> pfilter(Np=10000)))  # much tighter, and higher

## KEY POINT: too few particles doesn't just make the answer noisy -- it makes
## it systematically TOO LOW. More particles buys accuracy, not just precision.


# Exercise 4: read the filter's diagnostics -------------------------------

set.seed(1350254336)
pf <- measSIR |> pfilter(Np=5000)
plot(pf)                              # data / ess / cond.logLik, per week

## (a) the worst weeks by effective sample size
ess <- eff_sample_size(pf)
tibble(week = time(measSIR),
       reports = as.data.frame(measSIR)$reports,
       ess = ess) |>
  arrange(ess) |>
  head(5)
## Lowest at week 14 (ess ~ 700 of 5000 -- and remarkably STABLE run to run)
## and in the tail weeks 27-30 (the crash to near-zero cases). The tail ess is
## itself very noisy between runs -- it can land anywhere from ~100 to ~1500 --
## which is a hint that those weeks are exactly where the filter is on the edge.

## (b) what is the epidemic doing there?
## Week 14 is the EARLY RISE -- the data climb to 34 while the model, at these
## not-yet-fitted parameters, mostly hasn't taken off, so few particles predict
## it. Weeks 27-30 are the CRASH to 0-2 cases, which surprises particles that
## expect an ongoing epidemic.
## The big week-17->18 jump (18 -> 75) is NOT a problem: ess there is ~3200.
## The filter handles the surge fine; it struggles at the turning points.

## (c) raise Np and watch the worst week recover
essat <- function(np, seed=1350254336) {
  set.seed(seed); eff_sample_size(measSIR |> pfilter(Np=np))
}
w <- time(measSIR)
tibble(
  Np = c(5000, 20000),
  ess_wk14 = c(essat(5000)[w==14], essat(20000)[w==14])
)
## ess at week 14 goes from ~700 to ~2800 -- it scales with Np, but stays about
## the same FRACTION (~14%). More particles gives you more of them everywhere;
## the hard weeks stay proportionally hard.
##
## And the total log-likelihood creeps UP by a log unit or two -- the low-Np
## bias from Exercise 3 shrinking, exactly as expected. (The exact numbers
## depend on the seed; it is the DIRECTION that is reliable.)

## Q: A colleague reports a loglik at Np=1000 with ess under 200 at several
##    weeks. Do you trust it? NO -- with so few effective particles it is
##    almost certainly biased low and very noisy. Tell them to RAISE Np (and
##    report a logmeanexp over several runs, with a standard error -- Ex 5).
##
## This is a FILTER health check (do I have enough particles?). Using the same
## plot to criticise the MODEL -- which weeks it fits badly -- comes later, in
## the workflow lecture.


# Exercise 5: so run it many times, in parallel ----------------------------

registerDoRNG(1234)
foreach(i=1:20, .combine=c, .packages="pomp", .export="measSIR") %dopar% {
  measSIR |> pfilter(Np=5000)
} -> pfs

pfs |> logLik() |> logmeanexp(se=TRUE)
## est ~ -130, se ~ 0.9  (your numbers will differ slightly by machine)

## Why logmeanexp and not mean?
## The particle filter is UNBIASED for the likelihood L, not for log L.
## So we average on the likelihood scale and then take the log:
##   log( mean( exp(loglik) ) )
## logmeanexp does this without overflowing.


# Exercise 6: what does accuracy cost you? --------------------------------

Nps <- c(100, 500, 1000, 5000, 10000)
foreach(np=Nps, .combine=c, .packages="pomp", .export="measSIR") %dopar% {
  system.time(measSIR |> pfilter(Np=np))[3]
} -> times

tibble(Np=Nps, seconds=times) |>
  ggplot(aes(x=Np, y=seconds)) +
  geom_point(color="#0072B2", size=2) + geom_line(color="#0072B2") +
  labs(x="number of particles", y="elapsed seconds",
       title="Particle filter cost is linear in Np")

## Cost is LINEAR in Np: each particle is simulated and weighted once per
## observation. Doubling the particles doubles the work.
##
## The planning consequence: halving your Monte Carlo standard error costs
## roughly 4x the compute (se scales like 1/sqrt(Np)). Budget accordingly.


# Exercise 7 (stretch): a slice through the likelihood ---------------------

slice_design(
  center = coef(measSIR),
  Eta = rep(seq(from=0.01, to=0.2, length=30), each=3)
) -> p_eta

registerDoRNG(4321)
foreach(theta=iter(p_eta,"row"), .combine=rbind, .packages="pomp", .export="measSIR") %dopar% {
  measSIR |> pfilter(params=theta, Np=2000) -> pf
  theta$loglik <- logLik(pf)
  theta
} -> eta_slice

eta_slice |>
  ggplot(aes(x=Eta, y=loglik)) +
  geom_point(color="#0072B2", alpha=.6) +
  labs(x=expression(eta), y="log likelihood",
       title="Likelihood slice in the initial-susceptible fraction")

## The slice peaks at Eta ~ 0.049, with a log-likelihood of about -119.
## At our starting guess (Eta ~ 0.062) it is about -131.
##
## Look at what that means: a TINY move in Eta -- from 0.062 to 0.049 -- buys
## 12 log units of likelihood. The likelihood is extremely sensitive to Eta,
## and there is no way you could have eyeballed that from the simulation plot
## in Exercise 1. This is exactly why we need to fit rather than guess.
##
## CAUTION: a slice holds every other parameter FIXED. It is not a profile,
## and it should not be used to build a confidence interval.


# Cleanup -----------------------------------------------------------------

stopCluster(cl)
