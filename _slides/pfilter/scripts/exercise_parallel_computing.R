library(foreach)
library(doFuture)
plan(multisession,workers=8)    # number of cores

## -------- load model -------- ##
source("model_measSIR.R")

## ---- compute likelihood ---- ##

foreach(
  i=1:20, .combine="c", .options.future = list(seed = TRUE),
  .errorhandling = "remove"
) %dofuture% {
  # codes that you would like to run in parallel
  measSIR |> pfilter(Np=5000) |> logLik()
} -> pfs

pfs

logmeanexp(pfs,se=TRUE)
