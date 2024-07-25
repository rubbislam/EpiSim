## ----load-packages, echo=FALSE------------------------------------------------------------
library(tidyverse)
library(pomp)
library(ggplot2)
library(reshape2)


## ----read-data----------------------------------------------------------------------------
dat <- read_csv("1918flu_3cities_1wave.csv")
head(dat)


## ----plot, echo=FALSE, fig.width=6, fig.height=3, dpi=300---------------------------------
#| out-height: 90%
dat |>
  select(-date) |>
  reshape2::melt(id.vars = "week") |>
  ggplot(aes(x=week,y=value)) +
  geom_line() +
  scale_y_sqrt() +
  facet_wrap(variable ~ .) +
  labs(y="weekly cases",x="") +
  theme_minimal()


## ----meta-rproc---------------------------------------------------------------------------
sir_meta <- Csnippet("
  double *S = &S1, *I = &I1, *R = &R1, *D = &D1, *M = &M1;
  double *P = &P1, *H = &H1; int N[3] = {N1, N2, N3};
  double Beta, dN_SI, dN_IRD, dN_IR, dN_ID, dN_DM, dN_P;
  for (int i = 0; i < 3; i++) {
    Beta = beta0*pow(1-P[i]/N[i],kappa);
    dN_SI = rbinom(S[i],1-exp(-Beta*I[i]/N[i]*dt));
    dN_IRD = rbinom(I[i],1-exp(-gamma*dt));
    dN_IR = nearbyint((1-phi)*dN_IRD); dN_ID = nearbyint(phi*dN_IRD);
    dN_DM = rbinom(D[i],1-exp(-g*dt));
    dN_P = rbinom(P[i],1-exp(-lambda*dt));
    S[i] -= dN_SI; I[i] += dN_SI - dN_IRD; R[i] += dN_IR; 
    D[i] += dN_ID - dN_DM; M[i] += dN_DM; P[i] += dN_DM - dN_P;
    H[i] += dN_IRD;}
")


## ----meta-init----------------------------------------------------------------------------
sir_meta_rinit <- Csnippet("
  double *S = &S1, *I = &I1, *R = &R1, *D = &D1, *M = &M1;
  double *P = &P1, *H = &H1; int N[3] = {N1, N2, N3};
  double eta[3] = {eta1, eta2, eta3};
  double psi[3] = {psi1, psi2, psi3};
  for (int i = 0; i < 3; i++) {
    S[i] = nearbyint(N[i]*eta[i]); 
    I[i] = nearbyint(N[i]*psi[i]); 
    R[i] = nearbyint(N[i]*(1-eta[i]-psi[i]));
    D[i] = M[i] = P[i] = H[i] = 0;
  }
")


## ----meta-meas----------------------------------------------------------------------------
sir_meta_dmeas <- Csnippet("
  double lik1, lik2, lik3;
  lik1 = (ISNA(London)) ? dnbinom_mu(London,rho*H1,k,1) : 0;
  lik2 = (ISNA(Birmingham)) ? dnbinom_mu(Birmingham,rho*H2,k,1) : 0;
  lik3 = (ISNA(Liverpool)) ? dnbinom_mu(Liverpool,rho*H3,k,1) : 0;
  lik = lik1 + lik2 + lik3;
  lik = (give_log) ? lik : exp(lik);
")

sir_meta_rmeas <- Csnippet("
  London = rnbinom_mu(k,rho*H1);
  Birmingham = rnbinom_mu(k,rho*H2);
  Liverpool = rnbinom_mu(k,rho*H3);
")


## ----meta-pomp----------------------------------------------------------------------------
dat |> select(-date) |>
  pomp(
    times = "week", t0 = 0, 
    rprocess=euler(sir_meta,delta.t=1/7),
    rinit=sir_meta_rinit, rmeasure=sir_meta_rmeas,
    dmeasure=sir_meta_dmeas, accumvars = sprintf("H%d",1:3),
    statenames=c(sprintf("S%d",1:3),sprintf("I%d",1:3),
      sprintf("R%d",1:3), sprintf("D%d",1:3),sprintf("M%d",1:3),
      sprintf("P%d",1:3),sprintf("H%d",1:3)),
    paramnames=c(
      "beta0","kappa","gamma","phi","g","lambda","rho","k", 
      sprintf("N%d",1:3),sprintf("eta%d",1:3),sprintf("psi%d",1:3))
  ) -> pomp_meta


## ----meta-sim, eval=FALSE, purl=TRUE------------------------------------------------------
## params <- c(beta0 = 5, kappa = 5, gamma = 1/4, phi = 0.0119, g = 1/8,
##   lambda = 1/2, rho = 0.05, k = 10, N1 = 4484523, N2 = 919444,
##   N3 = 802940, eta1 = 0.2, eta2 = 0.2, eta3 = 0.2,
##   psi1 = 0.0005, psi2 = 0.0005, psi3 = 0.0005)
## 
## pomp_meta |>
##   simulate(
##     params = params, nsim=20,format="data.frame",include.data=TRUE
##   ) |>
##   select(week,.id,London,Birmingham,Liverpool) |>
##   reshape2::melt(id.vars = c("week",".id")) |>
##   ggplot(aes(x=week, y=value, group=.id, color=.id=="data")) +
##   geom_line() + facet_wrap(variable ~ .) +
##   theme_minimal() + guides(color="none")


## ----echo=FALSE, fig.width=6, fig.height=3, dpi=300---------------------------------------
#| out-height: 90%
params <- c(beta0 = 5, kappa = 5, gamma = 1/4, phi = 0.0119, g = 1/8, 
  lambda = 1/2, rho = 0.05, k = 10, N1 = 4484523, N2 = 919444, 
  N3 = 802940, eta1 = 0.2, eta2 = 0.2, eta3 = 0.2, 
  psi1 = 0.0005, psi2 = 0.0005, psi3 = 0.0005)

pomp_meta |> 
  simulate(
    params = params, nsim=20,format="data.frame",include.data=TRUE
  ) |>
  select(week,.id,London,Birmingham,Liverpool) |>
  reshape2::melt(id.vars = c("week",".id")) |>
  ggplot(aes(x=week, y=value, group=.id, color=.id=="data")) +
  geom_line() + facet_wrap(variable ~ .) + 
  theme_minimal() + guides(color="none")

