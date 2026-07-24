## =============================================================================
##  Lecture C -- slide code: genealogies and phylodynamics with phylopomp
##
##  The key runnable code BEHIND the slides. Open 10_phylopomp.Rproj first.
##  The trees here are simulated (for illustration); the figures/ scripts just
##  render polished versions of these same calls.
##
##  Install the latest phylopomp once (requires the current package):
##      # install.packages("pak")
##      pak::pak("kingaa/phylopomp")
## =============================================================================

library(phylopomp)
library(pomp)
library(tidyverse)

packageVersion("phylopomp")   # 0.19.4.2 or newer


# Simulate a genealogy from an SIR model ----------------------------------

set.seed(2024)
runSIR(
  time = 4,               # observation window, t = 0 .. 4
  Beta = 3, gamma = 1,    # transmission and recovery rates
  psi  = 0.4,             # per-lineage sampling rate (sets how many tips)
  pop  = 500,             # total host population
  S0 = 0.97, I0 = 0.03    # initial FRACTIONS of pop (not counts)
) -> G

plot(G, points = TRUE)    # the genealogy
plot(lineages(G))         # lineages through time

## Change the parameters partway through:
runSIR(time = 2, Beta = 3, gamma = 1, psi = 0.5,
       pop = 500, S0 = 0.98, I0 = 0.02) |>
  continueSIR(time = 5, Beta = 6, gamma = 2) -> G2
plot(G2)


# A richer model: two strains (SIIR) --------------------------------------

set.seed(7)
runSIIR(
  time = 3, Beta1 = 4, Beta2 = 8,   # strain 2 spreads faster
  gamma = 1.5, psi1 = 0.6, psi2 = 0.6,
  pop = 300, S_0 = 0.95, I1_0 = 0.03, I2_0 = 0.02
) -> GI

plot(GI, obscure = FALSE)           # colour each lineage by its strain


# Inference: the likelihood of a genealogy --------------------------------

## Same particle filter and logmeanexp as Lectures A and B; only the DATA
## changed, from a case-count series to a genealogy.
G |>
  sir_pomp(Beta = 3, gamma = 1, psi = 0.4,
           S0 = 0.97, I0 = 0.03, R0 = 0, pop = 500) |>
  pfilter(Np = 1000) |>
  replicate(n = 8) |> concat() -> pf

logmeanexp(sapply(pf, logLik), se = TRUE)

