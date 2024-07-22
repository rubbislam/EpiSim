## ----model-construct----------------------------------------------------------
library(tidyverse)
library(pomp)

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

seir_init <- Csnippet("
  S = nearbyint(Eta*N);
  E = 0; 
  I = 1;
  R = nearbyint((1-Eta)*N);
  H = 0;
")

dmeas <- Csnippet("
  lik = dnbinom_mu(reports,k,Rho*H,give_log);
")

rmeas <- Csnippet("
  reports = rnbinom_mu(k,Rho*H);
")

emeas <- Csnippet("
  E_reports = Rho*H;
")

read_csv("https://kingaa.github.io/sbied/pfilter/Measles_Consett_1948.csv") |>
  select(week,reports=cases) |>
  filter(week<=42) |>
  pomp(
    times="week",t0=0,
    rprocess=euler(seir_stoch,delta.t=1/7),
    rinit=seir_init,
    rmeasure=rmeas,
    dmeasure=dmeas,
    emeasure=emeas,
    accumvars="H",
    statenames=c("S","E","I","R","H"),
    paramnames=c("Beta","Sigma","Gamma","Rho","k","Eta","N"),
    params=c(Beta=40,Sigma=0.8,Gamma=0.5,Rho=0.5,k=10,Eta=0.06,N=38000)
  ) -> measSEIR
