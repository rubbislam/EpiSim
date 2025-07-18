library(foreach)
library(doParallel)
library(doRNG)                 # load doRNG for random number generation
registerDoParallel(cores=8)

## -------- load model -------- ##
source("model_measSIR.R")

## ---- compute likelihood ---- ##
foreach(
  i=1:20, .combine="c", .packages = "pomp", 
  .options.RNG = 1234          # set seed for RNG
) %dorng% {
  # codes that you would like to run in parallel
  measSIR |> pfilter(Np=5000)
} -> pfs

pfs |> logLik() |> logmeanexp(se=TRUE)
