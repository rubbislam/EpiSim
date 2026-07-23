library(foreach)
library(doParallel)
library(doRNG)                 # load doRNG for random number generation
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)

## -------- load model -------- ##
source("scripts/model_measSIR.R")

## ---- compute likelihood ---- ##
registerDoRNG(1234)            # set seed for RNG
foreach(
  i=1:20, 
  .combine=c, 
  .packages=c("pomp"), 
  .export = c("measSIR")
) %dopar% {
  # codes that you would like to run in parallel
  measSIR |> pfilter(Np=5000)
} -> pfs

stopCluster(cl)

pfs |> logLik() |> logmeanexp(se=TRUE)