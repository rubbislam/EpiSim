## =============================================================================
##  Lecture 4: Iterated filtering (mif2) -- EXERCISES
##
##  SISMID: Simulation-based Inference for Epidemiological Dynamics
##
##  HOW TO USE THIS FILE
##    1. Open 4_mle_mif2.Rproj first. That makes all the paths below work.
##    2. Run the Setup section once.
##    3. Work through the exercises in order. Answers are in
##       scripts/SOLUTIONS.R -- try first, peek after.
##
##  In the last session you learned to EVALUATE the likelihood with pfilter,
##  and saw that handing it to a standard optimizer fails. Now you MAXIMIZE it.
##
##  This script works on macOS, Linux and Windows without changes.
## =============================================================================


# Setup -------------------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)
library(reshape2)   # for melt() -- gives the column L1 used below

source("scripts/model_measSIR.R")

cl <- makePSOCKcluster(max(1, detectCores() - 1))
registerDoParallel(cl)

## IMPORTANT: build the parameter transformation ONCE, here, outside any
## foreach loop. Beta must stay positive; Rho and Eta must stay in [0,1].
## Doing this inside a parallel loop forces a recompile on every worker and
## is a common source of grief.
measSIR |>
  pomp(
    partrans=parameter_trans(log="Beta", logit=c("Rho","Eta")),
    paramnames=c("Beta","Rho","Eta")
  ) -> po


# Exercise 1: your first mif2 run -----------------------------------------
#
# Run ONE iterated filtering search and look at what it did.
#
# HINT:
#   po |> mif2(
#     Np=2000,                                   # particles
#     Nmif=50,                                   # iterations
#     cooling.fraction.50=0.5,                   # perturbation decay
#     rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02))
#   ) -> mf
#   plot(mf)
#
# Q: The top panel is the log-likelihood. Is it going up?
# Q: The other panels are parameters. Have they settled down, or are they
#    still drifting at iteration 50?
#
# NOTE: ivp() marks Eta as an "initial value parameter" -- it only affects
#       the start of the epidemic, so it is perturbed only at the beginning.




# Exercise 2: turn the knobs ----------------------------------------------
#
# You can feed a finished mif2 run straight back into mif2 to continue it,
# changing only what you want to change:
#
#   mf |> mif2(cooling.fraction.50=0.3) -> mf2
#   plot(mf2)
#
# Try each of these on its own and see what changes:
#   (a) cooling.fraction.50 = 0.1   vs  0.9
#   (b) rw.sd Beta = 0.002          vs  0.2
#   (c) Np = 200                    vs  5000
#
# Q: What does a too-small rw.sd look like? A too-large one?
# Q: cooling.fraction.50=0.5 means "after FIFTY iterations the perturbation
#    is half its starting size" -- not "half every iteration". Does the
#    behaviour you see match that?




# Exercise 3: many searches at once ---------------------------------------
#
# One search tells you very little. Run 20 independent searches in parallel
# and plot all their traces together.
#
# HINT:
#   foreach(i=1:20, .combine=c, .packages="pomp", .export="po",
#           .options.RNG=482947940) %dorng% {
#     po |> mif2(Np=2000, Nmif=50, cooling.fraction.50=0.5,
#                rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02)))
#   } -> mifs_local
#
#   mifs_local |> traces() |> melt() |>
#     ggplot(aes(x=iteration, y=value, group=L1)) +
#     geom_line(alpha=.3) +
#     facet_wrap(~name, scales="free_y")
#
# Q: Do all 20 agree on the log-likelihood?
# Q: Do all 20 agree on Beta? On Eta?
# Q: If the answers to those two questions are different -- why?
#    (Think back to the ridge you saw in the last session.)




# Exercise 4: the number mif2 prints is NOT your likelihood ----------------
#
# During a mif2 run the parameters are being deliberately perturbed. So the
# likelihood mif2 reports belongs to a NOISIER model than yours. To get the
# real number you must re-run a plain pfilter at the final parameters.
#
# HINT:
#   foreach(mf=mifs_local, .combine=rbind, .packages=c("pomp","dplyr"),
#           .options.RNG=900242057) %dorng% {
#     evals <- replicate(10, logLik(pfilter(mf, Np=5000)))
#     ll <- logmeanexp(evals, se=TRUE)
#     mf |> coef() |> bind_rows() |> bind_cols(loglik=ll[1], loglik.se=ll[2])
#   } -> results
#
#   results |> filter(loglik==max(loglik))
#
# Now compare against what mif2 itself reported:
#
#   sapply(mifs_local, logLik)
#
# Q: Which is higher? Is the gap in a consistent direction, or is it random?
# Q: Which parameter set is the best one you found?
#
# WORTH KNOWING: look at traces(mf) for any single run. It has Nmif+1 rows,
# and the loglik in the LAST row is NA. That is not a bug. mif2 records the
# likelihood of the parameters at the START of each iteration, so the final
# parameter vector -- the one you actually keep in coef(mf) -- has never been
# filtered at all. logLik(mf) is the likelihood of the previous, perturbed
# parameters. It is not the likelihood of your answer. That is the whole
# reason this exercise exists.




# Exercise 5 (stretch): does the starting point matter? -------------------
#
# Start several searches from DIFFERENT, deliberately bad starting points
# and see whether they all arrive at the same place.
#
# HINT: build a starting-parameter data frame with runif_design(), then
#       loop over its rows, passing params=<row> to mif2.
# HINT: ?runif_design
#
# Q: Do they converge to the same log-likelihood? To the same parameters?




# Exercise 6 (lunch break): does the latent period earn its keep? ----------
#
# In Lecture 2 you BUILT and SIMULATED an SEIR model -- the SIR with an extra
# "exposed" (E) compartment, so infected individuals spend a latent period
# before becoming infectious. Here you FIT it with mif2, exactly as you fit
# the SIR, and ask whether that latent period improves the fit enough to
# justify its extra parameter, Sigma.
#
# This launches 20 searches plus a re-filter, so start it before lunch.

source("scripts/model_measSEIR.R")   # -> measSEIR (Sigma = 1, a 1-week latent period)

# Build the transform ONCE. Sigma is a positive rate, so log-transform it too:
measSEIR |> pomp(
  partrans = parameter_trans(log = c("Beta","Sigma"), logit = c("Rho","Eta")),
  paramnames = c("Beta","Sigma","Rho","Eta")
) -> po_seir

# 20 local searches, now estimating Sigma alongside Beta, Rho, Eta:
foreach(i=1:20, .combine=c, .packages="pomp", .export="po_seir",
        .options.RNG=482947940) %dorng% {
  po_seir |> mif2(Np=2000, Nmif=50, cooling.fraction.50=0.5,
                  rw.sd=rw_sd(Beta=0.02, Sigma=0.02, Rho=0.02, Eta=ivp(0.02)))
} -> mifs_seir

# Honest re-estimate of each endpoint (as in Exercise 4), then keep the best:
foreach(mf=mifs_seir, .combine=rbind, .packages=c("pomp","dplyr"),
        .options.RNG=900242057) %dorng% {
  evals <- replicate(10, logLik(pfilter(mf, Np=5000)))
  ll <- logmeanexp(evals, se=TRUE)
  mf |> coef() |> bind_rows() |> bind_cols(loglik=ll[1], loglik.se=ll[2])
} -> seir_results

seir_results |> filter(loglik == max(loglik))    # your best SEIR fit

# Q: What is the best SEIR log-likelihood? The best SIR fit (Exercise 4) was
#    about -105.4. Is the SEIR better, and by how much?
# Q: The SEIR cost you ONE extra parameter (Sigma). A rough rule (AIC): an
#    extra parameter earns its keep only if it buys MORE than ~2 log-likelihood
#    units. Does the latent period clear that bar here?
# Q: Look at Sigma across the 20 searches: sort(seir_results$Sigma). Is it
#    pinned down, or does it wander? What does that say about whether these
#    data can actually see the latent period?




# Cleanup -----------------------------------------------------------------

stopCluster(cl)
