## =============================================================================
##  Lecture C: Genealogies and phylodynamics with phylopomp -- SOLUTIONS
##
##  Try scripts/EXERCISES.R first. Open 10_phylopomp.Rproj before running this.
##  Numbers below were produced with phylopomp 0.19.4.2 and the seeds shown;
##  phylopomp is stochastic, so set the seed to reproduce them exactly.
## =============================================================================


# Setup -------------------------------------------------------------------

library(phylopomp)
library(pomp)
library(tidyverse)


# Exercise 1: simulate a genealogy, and look at it ------------------------

set.seed(2024)
runSIR(time = 4, Beta = 3, gamma = 1, psi = 0.4,
       pop = 500, S0 = 0.97, I0 = 0.03) -> G

plot(G, points = TRUE)
plot(lineages(G))
getInfo(G, nsample = TRUE)$nsample          # 153 tips

## Raise / lower the sampling rate psi:
set.seed(2024)
getInfo(runSIR(time = 4, Beta = 3, gamma = 1, psi = 1.5,
               pop = 500, S0 = 0.97, I0 = 0.03), nsample = TRUE)$nsample  # ~ many more
set.seed(2024)
getInfo(runSIR(time = 4, Beta = 3, gamma = 1, psi = 0.1,
               pop = 500, S0 = 0.97, I0 = 0.03), nsample = TRUE)$nsample  # ~ far fewer

## ANSWER. psi is the SAMPLING rate: it controls how many of the infections
## end up as tips, not the epidemic itself. Bigger psi -> a denser tree;
## smaller psi -> a sparser one. The underlying rise-and-fall (driven by
## Beta, gamma, pop) is the same. psi changes how much of that shape you
## GET TO SEE, not the shape itself. This is exactly the observation model
## idea from Lecture A: the dynamics are hidden; sampling is how we observe.


# Exercise 2: two strains in one tree -------------------------------------

set.seed(7)
runSIIR(time = 3, Beta1 = 4, Beta2 = 8, gamma = 1.5,
        psi1 = 0.6, psi2 = 0.6, pop = 300,
        S_0 = 0.95, I1_0 = 0.03, I2_0 = 0.02) -> GI

plot(GI, obscure = FALSE)
plot(lineages(GI, obscure = FALSE))
getInfo(GI, nsample = TRUE)$nsample          # 102 tips

## ANSWER. Strain 2 (Beta2 = 8) spreads twice as fast as strain 1
## (Beta1 = 4), so it infects more hosts and leaves more, larger clades in
## the tree -- the olive lineages dominate the magenta ones. Sampling is
## equal (psi1 = psi2), so this is a difference in TRANSMISSION, not in
## observation: the faster strain simply has more descendants to sample.

## With equal transmission the advantage disappears and which strain "wins"
## in any single tree is down to chance (early stochastic events):
set.seed(7)
runSIIR(time = 3, Beta1 = 6, Beta2 = 6, gamma = 1.5,
        psi1 = 0.6, psi2 = 0.6, pop = 300,
        S_0 = 0.95, I1_0 = 0.03, I2_0 = 0.02) |>
  plot(obscure = FALSE)


# Exercise 3: the likelihood of a genealogy -------------------------------

loglik_at <- function(b) {
  G |>
    sir_pomp(Beta = b, gamma = 1, psi = 0.4,
             S0 = 0.97, I0 = 0.03, R0 = 0, pop = 500) |>
    pfilter(Np = 1000) |>
    replicate(n = 8) |> concat() -> pf
  logmeanexp(sapply(pf, logLik), se = TRUE)
}

grid <- seq(1.5, 5, by = 0.5)
map_dfr(grid, function(b) {
  set.seed(2024)                       # same filter randomness at each Beta
  e <- loglik_at(b)
  tibble(Beta = b, loglik = e[1], se = e[2])
}) -> ll

ll |> ggplot(aes(Beta, loglik)) +
  geom_vline(xintercept = 3, linetype = "dashed", colour = "#009E73") +
  geom_errorbar(aes(ymin = loglik - se, ymax = loglik + se), width = .08) +
  geom_line(colour = "#0072B2") + geom_point(colour = "#0072B2") +
  labs(x = expression(beta), y = "log-likelihood of the genealogy")

ll |> slice_max(loglik, n = 1)         # the MLE on the grid

## ANSWER. The likelihood peaks at Beta = 3, the value we simulated from:
##
##   Beta   loglik     se
##   1.5   -243.05   0.28
##   2.0   -231.02   0.06
##   2.5   -226.47   0.04
##   3.0   -226.44   0.03   <- maximum, and the truth
##   3.5   -229.26   0.07
##   4.0   -234.43   0.23
##   4.5   -242.10   0.59
##   5.0   -249.41   3.01
##
## Beta = 2.5 and 3.0 are within Monte Carlo error of each other -- the top
## of the ridge is flat, just as in Lecture B -- but the peak sits at the
## truth. This is the SAME particle filter and the SAME logmeanexp as
## Lecture A; only the data (a tree, not a case count) has changed. To turn
## this grid search into a real fit, hand the pomp object to mif2 (Lecture B).


# Exercise 4 (stretch): fit the WRONG model, and compare ------------------

## lbdp -- the linear birth-death process -- is the "phylodynamic null":
## constant per-capita birth (lambda), death (mu) and sampling (psi), and no
## susceptible pool. Simulate a tree from it, then fit the true lbdp and a
## wrong SIR, and compare their likelihoods.
set.seed(12)
runLBDP(time = 5, lambda = 2, mu = 1, psi = 1, n0 = 5) -> GL
getInfo(GL, nsample = TRUE)$nsample                 # 92 tips

## the TRUE model:
GL |>
  lbdp_pomp(lambda = 2, mu = 1, psi = 1, n0 = 5) |>
  pfilter(Np = 2000) |> replicate(n = 8) |> concat() -> pf_true
logmeanexp(sapply(pf_true, logLik), se = TRUE)      # about -63.0 (se 0.06)

## the WRONG model:
GL |>
  sir_pomp(Beta = 3, gamma = 1, psi = 1,
           S0 = 0.99, I0 = 0.01, R0 = 0, pop = 1000) |>
  pfilter(Np = 2000) |> replicate(n = 8) |> concat() -> pf_wrong
logmeanexp(sapply(pf_wrong, logLik), se = TRUE)     # about -76 (se 0.5)

## ANSWER. The true model wins: lbdp scores about -63.0, the misspecified SIR
## about -76 -- roughly 13 log units lower. A genealogy does not tell you its
## model, so you compare candidates by likelihood and prefer the higher one
## (rigorously, maximise each and compare with AIC; here we just evaluate at
## representative parameters). The SIR machinery, with its susceptible
## depletion, is simply the wrong story for a constant-rate birth-death tree.
##
## CAUTION -- direction matters. This ran cleanly because the lbdp tree is
## simple. Fitting lbdp (or SEIRS) to an SIR EPIDEMIC tree is fragile: with
## n0 below the tree's peak lineage count pfilter stops with an "invalid rate"
## error, and even above it the estimate is wildly unstable (se in the
## hundreds), because a constant-rate model cannot reproduce the rise-and-fall
## that susceptible depletion carves into an epidemic tree. Read a failed or
## unstable fit as information too: a model that cannot even evaluate the tree
## is telling you something, and the tools are still young.
