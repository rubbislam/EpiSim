library(tidyverse)
library(pomp)

## Read in the Measles Data
read_csv("../raw-data/Measles_Consett_1948.csv") |>
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
measSEIR <- meas |>
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
  )

## Specify the parameters
meas_seir_params <- c(Sigma = 0.5,
  Beta = 15,
  Gamma = 1,
  Rho = .95,
  k = 10,
  Eta = 0.1,
  N = 38000)

coef(measSEIR) <- meas_seir_params


# Exercises ---------------------------------------------------------------

# EXERCISE 1a:
# use one of the following functions to time how long a single particle
# filter iteration with Np = 1000 will take:
#    - system.time()   -- Single evaluation
#    - microbenchmark::microbenchmark()   -- multiple evaluations

# EXERCISE 1b:
# Based on the above calculation, approximately how long will a single
#   mif2() run with Np = 1000, Nmif = 50 take?
#   Try timing this and see how close your estimate was.

# EXERCISE 1c:
# Use the function parallel::detectCores() to see how many cores you have access to.
#   If you use all-but-one of your available cores, how many searches can you complete
#   in 20 mins on your local machine?

# EXERCISE 2a: Compute a small "global search" for the SEIR model parameters:
#    - Sigma, Beta, Eta, Rho
#
#    HINT: keep in mind your previous calculations. You may want to pick a number
#          of searches to be less than what can be done in 20 mins, for instance.
#          At home, you may want to try larger computations.

# EXERCISE 2b: Plot the "traces" of your global search. Are you converging to a
#   particular region in the likelihood surface?

# EXERCISE 3: Use pfilter to evaluate the likelihood of your global search. Do
#  you feel you have gotten close to finding the MLE? If you had more time, how
#  might you refine your search?

# EXERCISE 4 (at home): Set up a profile search for one of the parameters,
#    and compute a confidence interval.








# ==============================================================================
# IN-CLASS EXERCISES: Global Search and Profiling
# ==============================================================================
# Before starting, ensure you have loaded the pomp object 'measSIR' and the
# fixed parameters from the previous steps.
# ==============================================================================

library(pomp)
library(foreach)
library(doParallel)
library(doRNG)
library(microbenchmark)

# Assume measSIR and fixed_params are loaded in your environment.
# fixed_params <- c(N=38000, Gamma=2, k=10)

# ==============================================================================
# EXERCISE 1a: Benchmarking a single particle filter
# ==============================================================================
# Use one of the following functions to time how long a single particle
# filter iteration with Np = 1000 will take:
#    - system.time()   -- Single evaluation
#    - microbenchmark::microbenchmark()   -- multiple evaluations


# ==============================================================================
# EXERCISE 1b: Extrapolating to mif2
# ==============================================================================
# Based on the above calculation, approximately how long will a single
# mif2() run with Np = 1000, Nmif = 50 take?
# Try timing this and see how close your estimate was.

# STARTER CODE FOR MFI2:
# mif2(
#   measSIR,
#   Np = ...,
#   Nmif = ...,
#   cooling.fraction.50 = 0.5,
#   rw.sd = rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02))
# )


# ==============================================================================
# EXERCISE 1c: Estimating capacity
# ==============================================================================
# Use the function parallel::detectCores() to see how many cores you have access to.
# If you use all-but-one of your available cores, how many searches can you complete
# in 20 mins on your local machine?


# ==============================================================================
# EXERCISE 2a: Small Global Search
# ==============================================================================
# Compute a small "global search" for the SIR model parameters:
#    - Beta, Eta, Rho
#
# HINT: Keep in mind your previous calculations. Pick a number of searches (nseq)
# that can be completed well within your class time (e.g., 5-10 minutes).

# 1. Setup parallel backend
#
# FOR MAC:
#    registerDoParallel(cores = cores_to_use)
# For Windows:
#    cl <- makePSOCKcluster(cores_to_use)
#    registerDoParallel(cl)
#
# registerDoRNG(1234) # For reproducibility

# 2. Define the parameter box (lower and upper bounds)
# seir_box <- rbind(
#   Beta = c(____, ____),
#   Eta  = c(____, ____),
#   Rho  = c(____, ____)
# )

# 3. Generate starting guesses
# global_starts <- runif_design(lower = seir_box[,1], upper = seir_box[,2], nseq = _____)

# 4. Run the global search using foreach
# --- YOUR CODE HERE ---
# global_mifs <- foreach(guess = iter(_____, "row"), .combine = c, .packages = c("pomp")) %dopar% {
#   start_params <- unlist(c(_____, guess))
#
#   measSIR |> mif2(
#     params = start_params,
#     Np = _____,
#     Nmif = _____,
#     cooling.fraction.50 = 0.5,
#     rw.sd = rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02))
#   )
# }


# ==============================================================================
# EXERCISE 2b: Trace Plots
# ==============================================================================
# Plot the "traces" of your global search. Are you converging to a
# particular region in the likelihood surface?

# --- YOUR CODE HERE ---
# Hint: Use the traces() function on your list of mif2 objects (global_mifs).
# You may need to melt() the result and use ggplot2 as shown in the slides.
#
# library(ggplot2)
# global_mifs |>
#   traces() |>
#   melt() |>
#   ggplot(aes(x=iteration, y=value, group=.L1, color=factor(.L1))) +
#   geom_line() +
#   facet_wrap(~name, scales="free_y") +
#   guides(color="none")


# ==============================================================================
# EXERCISE 3: Evaluating the Global Search
# ==============================================================================
# Use pfilter to evaluate the likelihood of your global search. Do
# you feel you have gotten close to finding the MLE? If you had more time, how
# might you refine your search?

# --- YOUR CODE HERE ---
# Evaluate likelihoods for the final state of each mif2 run
# global_liks <- foreach(m = global_mifs, .combine = rbind, .packages=c("pomp")) %dopar% {
#
#   evals <- replicate(10, ...)
#   ll <- logmeanexp(evals, se=TRUE)
#
#   # Combine parameters and likelihood into a data frame row
#   m |> coef() |> bind_rows() |>
#     bind_cols(loglik=ll[1], loglik.se=ll[2])
# }


# ==============================================================================
# EXERCISE 4 (at home): Profile Search
# ==============================================================================
# Set up a profile search for one of the parameters (e.g., Eta),
# and compute a confidence interval.

# 1. Create a grid for your parameter of interest
# eta_grid <- seq(_____, _____, length=_____)

# 2. Generate profile design
# profile_guesses <- profile_design(
#   Eta = eta_grid,
#   lower = c(Beta = _____, Rho = _____),
#   upper = c(Beta = _____, Rho = _____),
#   nprof = _____
# )

# 3. Run the profile search
# --- YOUR CODE HERE ---
# profile_mifs <- foreach(guess = iter(_____, "row"), .combine = c, .packages = c("pomp")) %dopar% {
#   start_params <- unlist(c(_____, guess))
#
#   measSIR |> mif2(
#     params = start_params,
#     Np = ____, Nmif = ____, cooling.fraction.50 = 0.5,
#     # REMEMBER: Turn off the random walk for the parameter you are profiling!
#     rw.sd = rw_sd(___)
#   )
# }

# 4. Evaluate likelihoods and find the confidence interval cutoff using Wilks' theorem.
# --- YOUR CODE HERE ---


