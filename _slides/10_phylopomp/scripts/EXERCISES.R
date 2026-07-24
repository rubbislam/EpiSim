## =============================================================================
##  Lecture C: Genealogies and phylodynamics with phylopomp -- EXERCISES
##
##  SISMID: Simulation-based Inference for Epidemiological Dynamics
##
##  HOW TO USE THIS FILE
##    1. Open 10_phylopomp.Rproj first, so the paths below work.
##    2. Install the latest phylopomp (only needed once):
##         # install.packages("pak")
##         pak::pak("kingaa/phylopomp")
##    3. Run the Setup section once, then work through the exercises in order.
##       Each is meant to take a few minutes. Answers are in SOLUTIONS.R --
##       try first, peek after.
##
##  API NOTE (phylopomp >= 0.19.4): simulate with the per-model run...()
##  functions (runSIR, runSIIR, ...). The old simulate("SIR", ...) string
##  form is gone. Initial states (S0, I0, ...) are FRACTIONS of `pop`.
## =============================================================================


# Setup -------------------------------------------------------------------

library(phylopomp)
library(pomp)
library(tidyverse)

packageVersion("phylopomp")     # should be 0.19.4.2 or newer


# Exercise 1: simulate a genealogy, and look at it ------------------------
#
# Simulate an SIR genealogy, then view it two ways: the tree, and its
# lineages-through-time (how many branches are alive over time).

set.seed(2024)
runSIR(
  time = 4,
  Beta = 3, gamma = 1,
  psi  = 0.4,             # sampling rate: bigger psi -> more tips
  pop  = 500,
  S0 = 0.97, I0 = 0.03
) -> G

plot(G, points = TRUE)     # the genealogy
plot(lineages(G))          # lineages through time
getInfo(G, nsample = TRUE)$nsample   # how many tips did you get?

# TODO: raise psi to 1.5, then lower it to 0.1, re-simulate, and re-plot.
# Q: How does the number of sampled tips change? Does the OVERALL SHAPE
#    (the rise and fall of the lineages-through-time) change with psi, or
#    just how densely the tree is sampled?



# Exercise 2: two strains in one tree -------------------------------------
#
# runSIIR simulates two competing strains. With obscure = FALSE, phylopomp
# colours each lineage by the strain it belongs to.

set.seed(7)
runSIIR(
  time = 3,
  Beta1 = 4, Beta2 = 8,    # strain 2 spreads faster
  gamma = 1.5,
  psi1 = 0.6, psi2 = 0.6,
  pop = 300,
  S_0 = 0.95, I1_0 = 0.03, I2_0 = 0.02
) -> GI

plot(GI, obscure = FALSE)              # colour by strain
plot(lineages(GI, obscure = FALSE))    # lineages of each strain over time

# Q: Which strain leaves more (and larger) clades in the tree? Why?
# TODO: set Beta1 = Beta2 = 6 (equal strains) and re-simulate. Does either
#       strain still dominate, or is it now a coin toss?



# Exercise 3: the likelihood of a genealogy -------------------------------
#
# Inference is the SAME particle filter you used in Lectures A and B. Turn
# the tree into a pomp object with sir_pomp(), then pfilter() it. Because
# each pfilter is a noisy estimate, take several and combine with logmeanexp.

# A small helper: the log-likelihood of the SIR model at transmission rate b,
# given the genealogy G from Exercise 1.
loglik_at <- function(b) {
  G |>
    sir_pomp(Beta = b, gamma = 1, psi = 0.4,
             S0 = 0.97, I0 = 0.03, R0 = 0, pop = 500) |>
    pfilter(Np = 1000) |>
    replicate(n = 8) |> concat() -> pf
  logmeanexp(sapply(pf, logLik), se = TRUE)
}

loglik_at(3)     # at the TRUE Beta = 3

# TODO: evaluate loglik_at() across a grid of Beta and find the peak.
# HINT:
#   grid <- seq(1.5, 5, by = 0.5)
#   map_dfr(grid, \(b) { e <- loglik_at(b)
#                        tibble(Beta = b, loglik = e[1], se = e[2]) }) -> ll
#   ll |> ggplot(aes(Beta, loglik)) + geom_line() + geom_point()
# Q: Where is the maximum? Is it close to the true Beta = 3?



# Exercise 4 (stretch): fit the WRONG model, and compare ------------------
#
# A genealogy does not announce which model made it, so you can fit more than
# one and compare their likelihoods. The linear birth-death process (lbdp) is
# the "phylodynamic null": every lineage is born and dies at CONSTANT rates,
# with no susceptible pool to deplete.
#
# Simulate a tree FROM lbdp, then fit two models to it -- the true lbdp and a
# (wrong) SIR -- with the same pfilter you already know.

set.seed(12)
runLBDP(time = 5, lambda = 2, mu = 1, psi = 1, n0 = 5) -> GL   # lbdp is the TRUTH
getInfo(GL, nsample = TRUE)$nsample

# the TRUE model (lbdp):
GL |>
  lbdp_pomp(lambda = 2, mu = 1, psi = 1, n0 = 5) |>
  pfilter(Np = 2000) |> replicate(n = 8) |> concat() -> pf_true
logmeanexp(sapply(pf_true, logLik), se = TRUE)

# the WRONG model (SIR):
GL |>
  sir_pomp(Beta = 3, gamma = 1, psi = 1,
           S0 = 0.99, I0 = 0.01, R0 = 0, pop = 1000) |>
  pfilter(Np = 2000) |> replicate(n = 8) |> concat() -> pf_wrong
logmeanexp(sapply(pf_wrong, logLik), se = TRUE)

# Q: Which model has the higher log-likelihood, and by how much?
# Q: n0 is the population size at t = 0. lbdp has no susceptibles to run out,
#    so nothing slows it down. Look back at an SIR lineages-through-time plot
#    (Exercise 1): it RISES then FALLS. Could a constant-rate birth-death
#    process reproduce that fall?
#
# CAUTION: this comparison runs cleanly because the lbdp tree is simple. Going
# the OTHER way -- fitting lbdp (or SEIRS) to an SIR EPIDEMIC tree -- is
# numerically fragile: susceptible depletion stresses the mismatched model,
# and pfilter can stop with an "invalid rate" error or return a wildly
# unstable estimate. The inference tools are still young.
