library(tidyverse)
library(pomp)

## Read in the Measles Data
read_csv("raw-data/Measles_Consett_1948.csv") |>
  select(week,reports=cases) -> meas
meas |> as.data.frame() |> head(n=3)

# Code the SEIR model ----------------------------------------------
## Specify the process model
seir_stoch <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-Sigma*dt));
  double dN_IR = rbinom(I,1-exp(-Gamma*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

## Specify initial model conditions
seir_rinit <- Csnippet("
  S = nearbyint(Eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-Eta)*N);
  H = 0;
")

## Specify the measurement model
seir_dmeas <- Csnippet("
  lik = dnbinom_mu(reports, k, Rho*H, give_log);
")

seir_rmeas <- Csnippet("
  reports = rnbinom_mu(k, Rho*H);
")

## Create the SEIR pomp model
meas |>
  pomp(
    times = 'week',
    t0 = 0,
    rprocess = euler(seir_stoch, delta.t=1/7),
    rinit = seir_rinit,
    rmeasure = seir_rmeas,
    dmeasure = seir_dmeas,
    accumvars = "H",
    statenames = c("S","E","I","R","H"),
    paramnames = c("Sigma", "Beta","Gamma","N","Eta","Rho","k"),
    partrans = parameter_trans(
      log = c("Sigma", "Beta", "Gamma", "k"),
      logit = c("Eta", "Rho")
    )
  ) -> measSEIR

## Specify the parameters
c(Sigma = 0.5,
  Beta = 15,
  Gamma = 1,
  Rho = .95,
  k = 10,
  Eta = 0.1,
  N = 38000) -> meas_seir_params
