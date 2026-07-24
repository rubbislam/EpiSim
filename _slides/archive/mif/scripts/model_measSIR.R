## ----model-construct----------------------------------------------------------
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

dmeas <- Csnippet("
  lik = dnbinom_mu(reports,k,Rho*H,give_log);"
)

rmeas <- Csnippet("
  reports = rnbinom_mu(k,Rho*H);"
)

emeas <- Csnippet("
  E_reports = Rho*H;"
)

read_csv("../raw-data/Measles_Consett_1948.csv") |>
  select(week,reports=cases) |>
  filter(week<=42) |>
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_stoch,delta.t=1/7),
    rinit=sir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    emeasure=emeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","Gamma","Eta","Rho","k","N"),
    params=c(Beta=15,Gamma=2,Rho=0.5,k=10,Eta=0.06,N=38000)
  ) -> measSIR
