## ----boardingplot, echo=F, fig.dim = c(6, 3), out.height="60%"-----------------------------------------------------------------------------------
bsflu |> 
  ggplot(aes(day, B)) +
  geom_line()


## ----eval = F------------------------------------------------------------------------------------------------------------------------------------
## rproc <- Csnippet("
##   double N = 2000;
##   double t1 = rbinom(S,1-exp(-Beta*I/N*dt));
##   double t2 = rbinom(I,1-exp(-mu_I*dt));
##   S  -= t1;
##   I  += t1 - t2;
##   NI += t1;
##   R += t2;
## ")
## 
## rmeas <- Csnippet("
##   B = rpois(rho*NI+1e-6);
## ")


## ----prior-ex------------------------------------------------------------------------------------------------------------------------------------
priorDens <- Csnippet("
  lik = dunif(Beta, 1, 4, 1) +
        dunif(mu_I, 0.5, 3, 1) +
        dunif(rho, 0.5, 1, 1);
  if (!give_log) lik = exp(lik);
")


## ----echo = FALSE--------------------------------------------------------------------------------------------------------------------------------
source('bayes-example/boarding-school-model-incidence.R')

# Setting up model -----------------------------------------
## Let's build the model with the data
## We use euler process and simulate 12 steps per day
bsflu %>%
  select(day,B) %>%
  pomp(
    times="day", 
    t0=0,
    rmeasure=rmeas,
    dmeasure=dmeas,
    rprocess=euler(rproc, delta.t=1/12),
    rinit=rinit,
    partrans=parameter_trans(fromEst=fromEst,toEst=toEst),
    accumvars = 'NI',
    statenames=statenames,
    paramnames=paramnames
  ) -> flu

sim_params <- c(Beta = 2, ## Makes an R0 of 2 (beta/mu_I)
                mu_I = 1, ## average infectious period of 1 day (1/1)
                rho = 0.9, ## 90% of infections reported
                mu_R1 = 1/3)


## ----pmcmc-fxn, eval=FALSE-----------------------------------------------------------------------------------------------------------------------
## flu |> ## Standard pomp object that has already been created
##   pomp(dprior = priorDens, ## Prior specified from previous slide
##        params = sim_params, ## Parameter starting point
##        paramnames=c("Beta","mu_I","rho")) |> ## Parameter names
##   pmcmc(Nmcmc = 10000, ## Number of MCMC iterations
##         Np = 200, ## Number of particles to use
##         proposal = mvn_diag_rw(rw.sd = c(Beta=0.3, mu_I=0.3, rho=0.1))
##   ) -> test_mcmc


## ----pmcmc-fxn2, eval=FALSE----------------------------------------------------------------------------------------------------------------------
## flu |> ## Standard pomp object that has already been created
##   pomp(dprior = priorDens, ## Prior specified from previous slide
##        params = sim_params, ## Parameter starting point
##        paramnames=c("Beta","mu_I","rho")) |> ## Parameter names
##   pmcmc(Nmcmc = 10000, ## Number of MCMC iterations
##         Np = 200, ## Number of particles to use
##         proposal = mvn_diag_rw(rw.sd = c(Beta=0.3, mu_I=0.3, rho=0.1))
##   ) -> test_mcmc


## ----pmcmc-fxn2-1, echo=FALSE--------------------------------------------------------------------------------------------------------------------
bake(file="pmcmc1.rds",{
  flu |> ## Standard pomp object that has already been created 
    pomp(dprior = priorDens, ## Prior specified from previous slide
         params = sim_params, ## Parameter starting point
         paramnames=c("Beta","mu_I","rho")) |> ## Parameter names
    pmcmc(Nmcmc = 10000, ## Number of MCMC iterations
          Np = 200, ## Number of particles to use
          proposal = mvn_diag_rw(rw.sd = c(Beta=0.3, mu_I=0.3, rho=0.1))
    ) -> test_mcmc
  attr(test_mcmc,"ncpu") <- nbrOfWorkers()
  test_mcmc
}) -> test_mcmc
t_loc <- attr(test_mcmc,"system.time")
ncpu_loc <- attr(test_mcmc,"ncpu")


## ----autocorr------------------------------------------------------------------------------------------------------------------------------------
library(coda)

test_mcmc |> 
  traces() |>  
  autocorr.diag(lags=c(10, 50, 100))


## ----pmcmc-fxn3, eval=FALSE----------------------------------------------------------------------------------------------------------------------
## flu |>
##   pomp(dprior = priorDens,
##        params = sim_params, ##Using the sim params as a starting spot
##        paramnames=c("Beta","mu_I","rho")) |>
##   pmcmc(Nmcmc = 10000,
##         Np = 200,
##         proposal = mvn_rw(covmat(test_mcmc, thin = 50))
##   ) -> test_mcmc2


## ----pmcmc-fxn3-1, echo=FALSE--------------------------------------------------------------------------------------------------------------------
bake(file="pmcmc2.rds",{
  flu |> 
    pomp(dprior = priorDens,
         params = sim_params, ##Using the sim params as a starting spot
         paramnames=c("Beta","mu_I","rho")) |> 
    pmcmc(Nmcmc = 10000, 
          Np = 200,
          proposal = mvn_rw(covmat(test_mcmc, thin = 50))
    ) -> test_mcmc2
  attr(test_mcmc2,"ncpu") <- nbrOfWorkers()
  test_mcmc2
}) -> test_mcmc2
t_loc <- attr(test_mcmc2,"system.time")
ncpu_loc <- attr(test_mcmc2,"ncpu")


## ----autocorr2-----------------------------------------------------------------------------------------------------------------------------------
test_mcmc2 |> 
  traces() |>  
  autocorr.diag(lags=c(10, 50, 100))

