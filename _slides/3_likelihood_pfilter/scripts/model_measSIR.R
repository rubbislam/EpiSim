## The Consett 1948 measles outbreak: stochastic SIR model.
##
## Open the .Rproj file first, then:  source("scripts/model_measSIR.R")
## Paths below are relative to the project root, which is where the .Rproj lives.

library(tidyverse)
library(pomp)

sir_stoch <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-Gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;"
)

sir_init <- Csnippet("
  S = nearbyint(Eta*N);
  I = 1;
  R = nearbyint((1-Eta)*N);
  H = 0;"
)

sir_dmeas <- Csnippet("
  lik = dnbinom_mu(reports,k,Rho*H,give_log);"
)

sir_rmeas <- Csnippet("
  reports = rnbinom_mu(k,Rho*H);"
)

read_csv("data/Measles_Consett_1948.csv",show_col_types=FALSE) |>
  select(week,reports=cases) |>
  filter(week<=42) |>
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_stoch,delta.t=1/7),
    rinit=sir_init,
    rmeasure=sir_rmeas,
    dmeasure=sir_dmeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","Gamma","Eta","Rho","k","N"),
    params=c(Beta=15,Gamma=0.5,Rho=0.5,k=10,Eta=0.06,N=38000)
  ) -> measSIR

invisible(gc())
