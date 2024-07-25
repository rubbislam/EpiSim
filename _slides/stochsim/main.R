## ----doFuture,echo=FALSE,cache=FALSE------------------------------------------------------
library(doFuture)
library(microbenchmark)
plan(multisession)
set.seed(2488820)


## ----det-skel-----------------------------------------------------------------------------
library(pomp)
Csnippet("
        DS = -Beta*S*I/N;
        DI = Beta*S*I/N - I*Gamma;
        DR = I*Gamma;
        ") -> sir_det_skel


## ----meas-data1---------------------------------------------------------------------------
library(tidyverse)
read_csv(paste0("https://kingaa.github.io/sbied/stochsim/", 
  "Measles_Consett_1948.csv")) |> 
  select(week,reports=cases) -> meas
meas |> as.data.frame() |> head(n=3)


## ----meas-data2,echo=FALSE,fig.width=6, fig.height=3--------------------------------------
library(tidyverse)
meas |>
  ggplot(aes(x=week,y=reports)) +
  geom_line() +
  geom_point()


## ----rproc1R------------------------------------------------------------------------------
sir_stoch <- function (S, I, R, N, Beta, Gamma, delta.t, ...) {
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-Gamma*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  c(S = S, I = I, R = R)
}




## ----init1R-------------------------------------------------------------------------------
sir_rinit <- function (N, Eta, ...) {
  c(S = round(N*Eta), I = 1, R = round(N*(1-Eta)))
}


## ----pomp1R-------------------------------------------------------------------------------
library(pomp)
meas |>
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_stoch,delta.t=1/7),
    rinit=sir_rinit
  ) -> measSIR


## ----rproc2R------------------------------------------------------------------------------
sir_stoch <- function (S, I, R, N, Beta, Gamma, delta.t, H, ...) {
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-Gamma*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  H <- H + dN_IR
  c(S = S, I = I, R = R, H = H)
}

sir_rinit <- function (N, Eta, ...) {
  c(S = round(N*Eta), I = 1, R = round(N*(1-Eta)), H = 0)
}


## ----zero1R-------------------------------------------------------------------------------
measSIR |>
  pomp(
    rprocess=euler(sir_stoch,delta.t=1/7),
    rinit=sir_rinit, 
    accumvars="H"
  ) -> measSIR


## ----meas-modelR--------------------------------------------------------------------------
sir_dmeas <- function (reports, H, Rho, k, log, ...) {
  dnbinom(x=reports, size=k, mu=Rho*H, log=log)
}

sir_rmeas <- function (H, Rho, k, ...) {
  c(reports=rnbinom(n=1, size=k, mu=Rho*H))
}


## ----add-meas-modelR----------------------------------------------------------------------
measSIR |>
  pomp(
    rmeasure=sir_rmeas,
    dmeasure=sir_dmeas
  ) -> measSIR


## ----csnips-------------------------------------------------------------------------------
sir_stoch <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-Gamma*dt));
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")


## ----more-csnips--------------------------------------------------------------------------
    sir_rinit <- Csnippet("
      S = nearbyint(Eta*N);
      I = 1;
      R = nearbyint((1-Eta)*N);
      H = 0;
    ")
    sir_dmeas <- Csnippet("
      lik = dnbinom_mu(reports,k,Rho*H,give_log);
    ")
    sir_rmeas <- Csnippet("
      reports = rnbinom_mu(k,Rho*H);
    ")


## ----sir_pomp-----------------------------------------------------------------------------
measSIR |>
  pomp(
    rprocess=euler(sir_stoch,delta.t=1/7),
    rinit=sir_rinit,
    rmeasure=sir_rmeas,
    dmeasure=sir_dmeas,
    accumvars="H",
    statenames=c("S","I","R","H"),
    paramnames=c("Beta","Gamma","N","Eta","Rho","k")
  ) -> measSIR_C


## ----CSnippet-compare, echo=FALSE, fig.width=4, fig.height=4, cache=TRUE------------------
fR <- function() {measSIR |>
      simulate(
        params=c(Beta=7.5,Gamma=0.5,Rho=0.5,k=10,
          Eta=0.03,N=38000),
        nsim=100,format="data.frame",include.data=TRUE
      )} 
fC <- function() {measSIR_C |>
      simulate(
        params=c(Beta=7.5,Gamma=0.5,Rho=0.5,k=10,
          Eta=0.03,N=38000),
        nsim=100,format="data.frame",include.data=TRUE
      )} 
res <- microbenchmark(fR(), fC(), times=100L)

ggplot(res, aes(x=expr, y=time/1000000, color=expr)) +
  geom_boxplot(outliers = FALSE) +
  scale_color_brewer(palette="Dark2") +
  scale_x_discrete(labels = c("fR()" = "R", "fC()" = "CSnippet")) +
  scale_y_log10() + 
  labs(y = "time (milliseconds)", x = NULL) +
  theme_bw() + guides(color="none")


## ----sir_sims1, echo=FALSE----------------------------------------------------------------
measSIR |>
  simulate(
    params=c(Beta=7.5,Gamma=0.5,Rho=0.5,k=10, Eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE
  ) -> sims


## ----sir_sim1_plot,fig.width=8,fig.height=2-----------------------------------------------
measSIR |> 
  simulate(params=c(Beta=7.5,Gamma=0.5,Rho=0.5,k=10, Eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")


## ----fig.width=8, fig.height=1.5----------------------------------------------------------
measSIR |>
  simulate(params=c(Beta=25,Gamma=0.5,Rho=0.5,k=10,Eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")


## ----fig.width=8, fig.height=1.5----------------------------------------------------------
measSIR |>
  simulate(params=c(Beta=40,Gamma=0.5,Rho=0.5,k=10,Eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")


## ----fig.width=8, fig.height=1.5----------------------------------------------------------
measSIR |>
  simulate(params=c(Beta=40,Gamma=0.2,Rho=0.5,k=10,Eta=0.03,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line() + guides(color="none")


## ----fig.width=8, fig.height=1.5----------------------------------------------------------
measSIR |>
  simulate(params=c(Beta=15,Gamma=0.5,Rho=0.5,k=10,Eta=0.06,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")


## -----------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------
seir_init <- Csnippet("
  S = nearbyint(Eta*N);
  E = 0; I = 1;
  R = nearbyint((1-Eta)*N);
  H = 0;
")

measSIR |>
  pomp(
    rprocess=euler(seir_stoch,delta.t=1/7),
    rinit=seir_init,
    paramnames=c("N","Beta","Sigma","Gamma","Rho","Eta","k"),
    statenames=c("S","E","I","R","H")
  ) -> measSEIR


## ----fig.width=8, fig.height=1.3----------------------------------------------------------
measSEIR |>
  simulate(params=c(Beta=30,Sigma=0.8,Gamma=1.3,
                    Rho=0.5,k=10,Eta=0.06,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")


## ----fig.width=8, fig.height=1.3----------------------------------------------------------
measSEIR |> 
  simulate(params=c(Beta=40,Sigma=0.8,Gamma=1.3,
                    Rho=0.5,k=10,Eta=0.06,N=38000),
    nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")

