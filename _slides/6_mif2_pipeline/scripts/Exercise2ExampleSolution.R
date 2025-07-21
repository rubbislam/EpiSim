# Example solution for Exercise 2. 

library(pomp)
library(tidyverse)
library(foreach)
library(doParallel)
library(doRNG)

library(iterators)

# Build the model
source("model_measSIR.R")

# Transformations to stay "in-bounds": 
#   - log: parameters must be positive.
#   - logit: parameters in (0, 1)
measSIR <- measSIR |> pomp(
  partrans = parameter_trans(
    log=c("Beta", "Gamma", "k"),logit=c("Rho","Eta")
  ),
  paramnames = c("Beta","Gamma","Eta","Rho","k","N")
)


# Create starting region for estimated parameters: Beta, Gamma, Eta. 
freeze(seed=55266255,
       runif_design(
         lower=c(Beta=5,Gamma=0.2,Eta=0),
         upper=c(Beta=80,Gamma=5,Eta=0.99),
         nseq=1000
       )) |>
  mutate(
    Rho=0.6, k=10, N=38000
  ) -> guesses


# Set up parallelization code ---------------------------------------------

## MAC / Linux 
# registerDoParallel(detectCores() - 1)


## Windows
cl <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)  # DON'T FORGET TO CLOSE CLUSTER! 


# Here we will adopt a *simulated tempering* approach 
# (following a metallurgical analogy), in which we increase the size 
# of the random perturbations some amount (i.e., "reheat"), and then continue 
# cooling.


bake(file="global_search_gamma.rds",
     dependson=guesses,{
       registerDoRNG(610408)
       foreach(
         guess=iter(guesses,"row"), .combine=rbind) %dopar% {
           
           # Initial search and cooldown
           measSIR |>
             mif2(params=guess, Np=2000, Nmif=100,
                  cooling.fraction.50=0.5,
                  rw.sd=rw_sd(Beta=0.02,Gamma=0.02,Eta=ivp(0.02))) -> mf
           
           # Restarting the search and "warm-up".
           mf |>
             mif2(
               Nmif=100,rw.sd=rw_sd(Beta=0.01,Gamma=0.01,Eta=ivp(0.01))
             ) |>
             mif2(
               Nmif=100,
               rw.sd=rw_sd(Beta=0.005,Gamma=0.005,Eta=ivp(0.005))
             ) -> mf
           
           # Calculate the likelihood of final parameter estimate
           replicate(
             10,
             mf |> pfilter(Np=5000) |> logLik()
           ) |> logmeanexp(se=TRUE) -> ll
           
           # Combine parameter estimate with likelihood values, return result.
           mf |> coef() |> bind_rows() |>
             bind_cols(loglik=ll[1],loglik.se=ll[2])
         } -> results
       results
     }) |>
  filter(is.finite(loglik)) -> results


# If we want to save to a database with parameter estimates: 
# read_csv("measles_params.csv") |>
results |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

# If we want to visualize the results: 
read_csv("measles_params.csv") |>
   filter(loglik>max(loglik)-20) -> all

pairs(~loglik+Rho+Gamma+Beta+Eta,data=all,pch=16,cex=0.3,
      col=if_else(round(all$Rho,3)==0.6,1,4))


# Profiling over Gamma ----------------------------------------------------

# Since we are interested in gamma, we can make a profile. 

# First step is to make a "search-box" for all other parameters. If we saved 
# our results, we can do it via: 
read_csv("measles_params.csv") |>
  filter(
    loglik>max(loglik)-20,
    loglik.se<2,
    abs(Rho-0.6)<0.01
  ) |>
  sapply(range) -> box

# Sample starting points for profile search from box
freeze(seed=610408798,
       profile_design(
         Gamma=seq(0.2,2,by=0.1),
         lower=box[1,c("Beta","Eta")],
         upper=box[2,c("Beta","Eta")],
         nprof=100, type="runif"
       )) |>
  mutate(
    N=38000,  # Add these fixed parameters to each row
    Rho=0.6,  # Add these fixed parameters to each row
    k=10      # Add these fixed parameters to each row
  ) -> guesses

# Conduct 
bake(file="gamma_profile1.rds",
     dependson=guesses,{
       
       foreach(guess=iter(guesses,"row"), .combine=rbind) %dopar% {
         measSIR |>
           mif2(params=guess, Np=2000, Nmif=100, cooling.fraction.50=0.5,
                rw.sd=rw_sd(Beta=0.02,Eta=ivp(0.02))
           ) |>
           mif2(Nmif=100) |>
           mif2(Nmif=100,rw.sd=rw_sd(Beta=0.01,Eta=ivp(0.01))) |>
           mif2(Nmif=100,rw.sd=rw_sd(Beta=0.005,Eta=ivp(0.005))) -> mf
         replicate(10,mf |> pfilter(Np=5000) |> logLik()) |>
           logmeanexp(se=TRUE) -> ll
         mf |> coef() |> bind_rows() |>
           bind_cols(loglik=ll[1],loglik.se=ll[2])
       } -> results
       
         results
     }) |>
  filter(is.finite(loglik)) -> results

# Add our results to the existing results. 
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

# Plot results: 
results |>
group_by(round(Gamma,2)) |>
  filter(rank(-loglik)<=1) |>
  ungroup() |>
  ggplot(aes(x=Gamma,y=loglik))+
  geom_point() + xlab(expression(gamma)) +
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )

# WINDOWS: 
# stopCluster(cl)