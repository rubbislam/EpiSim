## The Consett 1948 measles outbreak: stochastic SEIR model -- the SIR with a
## latent (exposed) compartment E, entered at rate Sigma. Sigma = 1 gives a
## one-week (7-day) mean latent period.
##
## Parameters follow the A/B SIR convention (Gamma=0.5, Rho=0.5, Eta=0.06), so
## the SEIR nests against our SIR for a fair likelihood comparison (B's
## Exercise 6). Lecture 2 used slightly different values (Beta=7.5, Rho=0.95,
## Eta=0.1); that is a different lesson (simulation, not fitting), and the
## structure is identical, so this still works as a debugging reference.
##
## Open the .Rproj file first, then:  source("scripts/model_measSEIR.R")
## Paths below are relative to the project root, which is where the .Rproj lives.

library(tidyverse)
library(pomp)

seir_stoch <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-Sigma*dt));
  double dN_IR = rbinom(I,1-exp(-Gamma*dt));
  S -= dN_SI;
  E += dN_SI - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

seir_rinit <- Csnippet("
  S = nearbyint(Eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-Eta)*N);
  H = 0;
")

seir_dmeas <- Csnippet("
  lik = dnbinom_mu(reports, k, Rho*H, give_log);
")

seir_rmeas <- Csnippet("
  reports = rnbinom_mu(k, Rho*H);
")

read_csv("data/Measles_Consett_1948.csv",show_col_types=FALSE) |>
  select(week,reports=cases) |>
  filter(week<=42) |>
  pomp(
    times = "week", t0 = 0,
    rprocess = euler(seir_stoch, delta.t=1/7),
    rinit = seir_rinit,
    rmeasure = seir_rmeas,
    dmeasure = seir_dmeas,
    accumvars = "H",
    statenames = c("S","E","I","R","H"),
    paramnames = c("Sigma", "Beta","Gamma","N","Eta","Rho","k"),
    params = c(Sigma=1, Beta=15, Gamma=0.5, N=38000, Eta=0.06, Rho=0.5, k=10)
  ) -> measSEIR

invisible(gc())