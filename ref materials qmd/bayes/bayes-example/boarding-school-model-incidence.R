## Code to specify the model for the influenza boarding school example
## Drawn from Seb's github repo: https://github.com/sbfnk/inference_pomp_rbi/tree/master
library(pomp)

## Assumes bed usage is poisson distributed
## around the number of newly infectious kids
## multiplied by a reporting rate (rho)
#### Added a small positive number ($10^{-6}$) to the expected incidence 
#### to prevent problems if $NI$ goes to zero
dmeas <- Csnippet("
  lik = dpois(B, rho*NI+1e-6, give_log);
")

## Simulating from the same measurement model described above
rmeas <- Csnippet("
  B = rpois(rho*NI+1e-6);
")

## Specification of the process model simulator
## Specifies transitions from Susceptible to Infectious to Recovering
## Basic idea is the people are either susceptible (S), initially infectious (I), bed-confined (R1), 
## out of bed but symptomatic (R2), or fully recovered (R3)
## In this example we are altering the model a bit so that instead of fitting to the number of
## currently bed-confined boys, we assume the data are actually the newly infectious incidence
## This is done because the statistical models assume independence across data points, but
## bed confinement is a measure of prevalence.
## since R2 and R3 do not impact transmission, we only need to track S, I, NI, and R1 for our model
## See for more info on model: https://github.com/sbfnk/inference_pomp_rbi/tree/master
rproc <- Csnippet("
  double N = 2000;
  double t1 = rbinom(S,1-exp(-Beta*I/N*dt));
  double t2 = rbinom(I,1-exp(-mu_I*dt));
  S  -= t1;
  I  += t1 - t2;
  NI += t1;
  R += t2;
")

## Assumes 1 initially infected individual of population of 2,000 boys
rinit <- Csnippet("
 S = 1999;
 I = 1;
 NI = 0;
 R = 0;
")

## We carry out log transofrmations for rate parameters Beta and mu_I, as they must be strictly positive
## Logit transformations for rho to keep it between 0 and 1
## The toEst takes from the real scale and converts to scale used for parameter estimation
#### e.g. for beta, takes from strictly positive scale and converts to -Inf to Inf scale
#### e.g. for rho, takes from 0-1 scale and converts to -Inf to Inf scale
toEst <- Csnippet("
 T_Beta = log(Beta);
 T_mu_I = log(mu_I);
 T_rho = logit(rho);
")

## This does the opposite, takes from parameter estimation scale and converts to real parameter scale
#### e.g. for beta, takes from -Inf to Inf scale and converts to strictly positive scale
#### e.g. for rho, takes from -Inf to Inf scale and converts to 0-1 scale 
fromEst <- Csnippet("
 Beta = exp(T_Beta);
 mu_I = exp(T_mu_I);
 rho = expit(T_rho);
")

## We are tracking three states, and there are four parameters in our model
statenames <- c("S", "I", "NI","R")
paramnames <- c("Beta", "mu_I", "rho")




