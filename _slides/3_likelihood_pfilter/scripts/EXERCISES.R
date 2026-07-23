## =============================================================================
##  Lecture 3: The particle filter -- EXERCISES
##
##  SISMID: Simulation-based Inference for Epidemiological Dynamics
##
##  HOW TO USE THIS FILE
##    1. Open 3_likelihood_pfilter.Rproj first. That makes all the paths below work.
##    2. Run the Setup section once.
##    3. Work through the exercises in order. Each one is meant to take a few
##       minutes. Answers are in scripts/SOLUTIONS.R -- try first, peek after.
##
##  This script works on macOS, Linux and Windows without changes.
## =============================================================================


# Setup -------------------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)

## Build the Consett 1948 measles model. This gives you `measSIR`.
source("scripts/model_measSIR.R")

## What parameters is it currently using?
coef(measSIR)

## Start a parallel cluster. makePSOCKcluster works on every platform.
## Leave one core free so your laptop stays usable.
cl <- makePSOCKcluster(max(1, detectCores() - 1))
registerDoParallel(cl)

## NOTE: run stopCluster(cl) at the very end of the session (see the bottom).


# Exercise 1: does the model look anything like the data? ------------------
#
# Before you compute any likelihood, LOOK. Simulate 20 outbreaks from the
# model and plot them against the real Consett data.
#
# Q: Do the simulations look like the data? Where do they get it wrong?
#
# HINT: measSIR |> simulate(nsim=20, format="data.frame", include.data=TRUE)
#       gives you a data frame with a .id column; .id=="data" is the real data.
# HINT: ggplot(aes(x=week, y=reports, group=.id, color=(.id=="data")))




# Exercise 2: your first likelihood ----------------------------------------
#
# Run ONE particle filter and get a log-likelihood.
#
# HINT: pf <- measSIR |> pfilter(Np=1000)
#       logLik(pf)
#
# Q: You should get something around -130. Is that number good or bad?
#    (Trick question -- a log-likelihood on its own means nothing. It only
#     means something COMPARED to another one.)




# Exercise 3: run it again ------------------------------------------------
#
# Run the exact same line from Exercise 2 a second time. And a third.
#
# Q: Do you get the same number? Why not?
# Q: Try Np=100 and then Np=10000. What happens to how much the answer
#    jumps around? What happens to the number itself?
#
# This is the single most important thing to understand here: the particle
# filter does not return THE likelihood. It returns a random ESTIMATE of it.




# Exercise 4: read the filter's diagnostics -------------------------------
#
# WHERE does that randomness come from? Every pfilter object carries per-week
# diagnostics. The quickest look is one command:
#
#   pf <- measSIR |> pfilter(Np=5000)
#   plot(pf)
#
# You get three stacked panels sharing the week axis:
#   - reports     : the data
#   - ess         : the EFFECTIVE SAMPLE SIZE -- how many particles really
#                   count that week (out of Np)
#   - cond.logLik : that week's own contribution to the log-likelihood
#
# (a) Pull the effective sample size out and find the worst weeks:
#
#       ess <- eff_sample_size(pf)
#
#     HINT: build a tibble with the week, the reports and the ess, then
#           arrange() by ess. Which 3-4 weeks are lowest?
#
# (b) Look at the DATA at those weeks (Exercise 1's plot, or just the numbers).
#     What is the epidemic doing there? Is it where you would have guessed the
#     filter would struggle? (Check the big week-17->18 jump too -- is IT low?)
#
# (c) The effective sample size is your "do I have enough particles?" gauge.
#     Re-run with Np=20000 and look at the ess at the worst week from part (a).
#     Did it recover? Did the TOTAL log-likelihood move?
#
# Q: A colleague reports a log-likelihood computed with Np=1000, and its ess
#    is under 200 at several weeks. Do you trust the number? What one thing
#    would you tell them to change?




# Exercise 5: so run it many times, in parallel ----------------------------
#
# Since each pfilter is independent, run 20 of them across your cores and
# combine them into one estimate with a standard error.
#
# HINT:
#   registerDoRNG(1234)                      # reproducible parallel seeds
#   foreach(i=1:20, .combine=c, .packages="pomp", .export="measSIR") %dopar% {
#     measSIR |> pfilter(Np=5000)
#   } -> pfs
#   pfs |> logLik() |> logmeanexp(se=TRUE)
#
# Q: Why logmeanexp() and not just mean()?
#    (Read the name backwards: log(mean(exp(x))). We average the LIKELIHOODS,
#     then take the log -- because the pfilter is unbiased for the likelihood,
#     not for the log-likelihood.)




# Exercise 6: what does accuracy cost you? --------------------------------
#
# How much processing time does a particle filter take, and how does it
# scale with the number of particles?
#
# Form a guess FIRST, then measure it.
#
# HINT: system.time( <your code> )[3] returns elapsed seconds.
# HINT: try Np = c(100, 500, 1000, 5000, 10000) and plot time against Np.
#
# Q: Is the relationship linear? Does that match your guess?
# Q: You have a 90-minute coffee break and 8 cores. How many particle
#    filters can you afford?




# Exercise 7 (stretch): a slice through the likelihood ---------------------
#
# Fix every parameter except Eta, vary Eta, and compute the likelihood at
# each value. This is a likelihood "slice".
#
# HINT: slice_design(center=coef(measSIR), Eta=seq(0.01, 0.2, length=30))
#       then loop over the rows with foreach, calling pfilter on each.
# HINT: ?slice_design
#
# Q: Where is the peak? Is it where you expected?




# Cleanup -----------------------------------------------------------------

## Always shut the cluster down when you are done.
stopCluster(cl)
