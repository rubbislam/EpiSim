library(tidyverse)
library(pomp)
library(cowplot)
library(coda)
theme_set(theme_cowplot())

## This reads in the main components of the pomp model
source('ref materials qmd/bayes/bayes-example/boarding-school-model-incidence.R')

# Checking the data -------------------------------------------------------
## Quickly look at data to see if it looks good
## In actuality B corresponds to a prevalence metric (total number of boys bed-confined)
## In our example we assume it's incidence and make the population larger
## This is because there are statistical issues with fitting to prevalence estimates 
## (similar issues to fitting to cumulative counts)
bsflu |>
  ggplot(aes(day, B)) + 
  geom_line() +
  background_grid(major = 'xy', minor = 'xy')


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


# Validating the model ----------------------------------------------------

## Always good to simulate data to make sure model looks okay
## Specify some reasonable numbers for each of the four parameters:
## R0 of 2, (beta/mu_I = 2, and 
##           1/mu_I meaning people are infectious for 1 day before being bed-confined)
## 90% of students are reported
## People are in bed on average for 3 days 
sim_params <- c(Beta = 2, ## Makes an R0 of 2 (beta/mu_I)
                mu_I = 1, ## average infectious period of 1 day (1/1)
                rho = 0.9, ## 90% of infections reported
                mu_R1 = 1/3)

## Run 100 simulations
flu |>  
  simulate(params = sim_params,
           nsim = 100,
           format = "data.frame") -> sims

## Plot the simulations to inspect
## Also add data as red points to compare
head(sims)
sims |> 
  ggplot(aes(x=day,y=B,group=.id))+
  geom_line(alpha = .2) +
  geom_point(data = bsflu, aes(day, B), inherit.aes=F, color = 'red') 

## Simulations aren't too bad compared with data, so let's make sure pfilter works
## pfilter is main workhorse, so it's good to make sure it's working before more complex analyses
## We assess the likelihood at the parameters used for simulating
flu |> 
  pfilter(params=sim_params,
          Np=1000) -> test_pf
plot(test_pf)

## These parameters provide decent counts for the ess and loglik isn't too bad 
## Try changing some parameters to see what impact they have?


# Setting up bayesian estimation ---------------------------------------------

## We now can specify the prior distribution for the three parameters in the model
## That we will estimate
## We use uniform distributions for all, and choose reasonable minima and maxima
#### note that you will likely want to run subsequent analyses to test the sensitivity
#### of the results to these priors, but okay for now to try out
priorDens <- Csnippet("
  lik = dunif(Beta, 1, 4, 1) +
        dunif(mu_I, 0.5, 3, 1) +
        dunif(rho, 0.5, 1, 1);
  if (!give_log) lik = exp(lik);
")

## Run a short pmcmc run with simulated parameters to make sure things look okay
## Uses the flu pomp model and adds the prior density, initial parameters we used before, 
## It also changes paramnames to only those parameters that are being estimated (not fixed)
## Nmcmc is the number of mcmc samples to obtain
## Np is the number of particles to use 
## proposal indicates how mcmc parameter search happens, here we are using a
## simple random walk with standard deviations for each of the parameters
## Can play around with the random walk to see how it impacts things, but
## you are looking for a reasonable starting place for decent chain mixing

flu |> 
  pomp(dprior = priorDens,
       params = sim_params, ##Using the sim params as a starting spot
       paramnames=c("Beta","mu_I","rho")) |> 
  pmcmc(Nmcmc = 10000, 
        Np = 200,
        proposal = mvn_diag_rw(rw.sd = c(Beta=0.3, mu_I=0.3, rho=0.1))
  ) -> test_mcmc

## A few notes, here you are just looking to see some movement of chains, and that the
## loglik is increasing roughly each iteration - in our case our initial parameters 
## are pretty good so not as much loglik movement
## in general increasing the particles increases the acceptance rate
## Decreasing the standard deviation of the parameter random walk does the same.
test_mcmc |> traces() |> rejectionRate()
test_mcmc |> traces() |> traceplot() ##traces function turns the pomp output into mcmc object for coda


## These don't look "good" but they are sufficient 
## to give us confidence to move to main chain runs
## We are going to change our mixing to be based on the empirical acceptance variance/covariance matrix
## of the previous chains


# Running the bayesian estimation -----------------------------------------
library(foreach)
library(doParallel)
registerDoParallel()
library(lhs) ## Uses latin hypercube sampling to ensure reasonable parameter starting locations

set.seed(34631)
parm_lhs <- randomLHS(n = 5, k = 3) ## We are running 5 chains and have 3 parameters in the model

## Each parameter (column) we need to use the generated quantile to sample from our prior uniform distribution
## This looks at a good range of possible parameter values, so they won't necessarily have the good starting spot
## We were using before.
parm_lhs[,1] <- qunif(parm_lhs[,1], 1, 4) ## Beta
parm_lhs[,2] <- qunif(parm_lhs[,2], 0.5, 3) ## mu_I
parm_lhs[,3] <- qunif(parm_lhs[,3], 0.5, 1) ## rho

## give column names, convert to tibble, and add in the fixed parameter mu_R1
colnames(parm_lhs) <- names(sim_params)[1:3]
parm_lhs |> 
  as_tibble() -> parm_lhs


bake('ref materials qmd/bayes/bayes-example/pmcmc-chains.rds',
     seed=13434,
     {
       foreach (curr_params=iter(parm_lhs,"row"), 
                .packages='pomp',
                .inorder=FALSE) %dopar% {
                  flu |> 
                    pomp(dprior = priorDens,
                         params = curr_params,
                         paramnames=c("Beta","mu_I","rho")) |> 
                    pmcmc(
                      Nmcmc=11000,
                      Np = 300,
                      proposal=mvn_rw(covmat(test_mcmc, thin = 50))
                    ) -> flu_pmcmc ## Always save the pmcmc object so you can add samples later and simulate, etc.
                }
     }) -> pmcmc_chains

## Thinning 50 looks fine (though 100 might be better and we may want to run 
## more Nmcmc samples for each chain if this were for publication)
pmcmc_chains |> purrr::map(traces) |>  purrr::map(autocorr.diag, lags=c(10, 50, 100))
pmcmc_chains |> purrr::map(traces) |> traceplot()
## Can see in this traceplot how loglik increases over iterations
## Very quickly find the same locations though which is good!

## Remove the burn-in and thin the chains, we will have 200 total samples per chain and
## 1,000 total posterior samples across all 5 chains
pmcmc_chains |> 
  purrr::map(traces) |> 
  as.mcmc.list() |> 
  window(start=1000, end=11000, thin=50) -> posteriors

## Look at the final traceplots -- all variables look pretty consistent and solid for all chains
posteriors |> traceplot()

## Convergence check - want to see values as close to 1 as possible
## All values are good for the variables we estimated
posteriors |> gelman.diag(confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = FALSE)


## We convert into a tibble for using to plot
posteriors |> 
  purrr::map(as_tibble) |> 
  bind_rows(.id = 'chain') |> 
  group_by(chain) |> 
  mutate(samp = seq_along(loglik)) -> processed_chains


# Get posterior parameter estimates ---------------------------------------

## Posterior density plotting
## Looking for a reasonably smooth density for each of the estimated parameters
posteriors |> densplot()

## Posterior parameter estimates for common quantiles
posteriors |> summary()

## Calculate estimated reproduction number with credible intervals
posteriors |> purrr::map(as_tibble) |> bind_rows() |> 
  mutate(rnot = Beta/mu_I) |> 
  summarize(lo = quantile(rnot, probs = 0.025),
            med = quantile(rnot, probs = 0.5),
            avg = mean(rnot),
            hi = quantile(rnot, probs = 0.975))

# Compare model fit -------------------------------------------------------

## This uses the filtering process to estimate the latent state values
## We use the same chain start number (1000) and thinning value (50)
## See this for more info: https://kingaa.github.io/pomp/vignettes/getting_started.html#Sampling_the_posterior_using_particle_Markov_chain_Monte_Carlo
pmcmc_chains |>
  purrr::map(filter_traj)  |> 
  purrr::map(pomp::melt) |> 
  bind_rows(.id = 'chain') |> 
  filter(rep>1000, rep %% 50 == 0) |>  ## Removes the burn-in and thins with same parameters as above
  mutate(day = time-1) |>  ## Sets the time to the correct day for the influenza outbreak can use time() pomp function if needed
  pivot_wider(names_from =  name) |> ## pulls state counts into columns
  group_by(day) |> ## We are summarizing by the day of simulation
  reframe(label=c("lo", "med", "hi"),
            p=c(0.025, 0.5, 0.975),
            q=wquant(NI, probs=p)) |> ## Summarizes for every day the 95% credible interval for the states accounting for parameter uncertainty
  select(-p) |> 
  pivot_wider(names_from = label, values_from = q) -> in_samp_quants
in_samp_quants
  

## Plot the model alongside the data
in_samp_quants |> 
  ggplot(aes(day, med), color = 'blue') +
  geom_line() +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2, fill = 'blue') +
  geom_point(data = bsflu, aes(x = day, y = B), inherit.aes = F)

## Nice model fit including parameter uncertainty!

