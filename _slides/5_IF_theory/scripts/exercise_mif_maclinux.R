library(foreach)
library(doParallel)
library(doRNG)

registerDoParallel(cores=detectCores())

## -------- load model -------- ##
source("model_measSIR.R")

## ------ 1. mif basics ------- ##
measSIR |>
  mif2(
    Np=2000,                      # number of particles
    Nmif=50,                      # number of iteration
    cooling.fraction.50=0.5,      # decay of perturbation
    rw.sd=rw_sd(                  # perturbation
      Beta=0.02, Rho=0.02, Eta=ivp(0.02)
    ),
    partrans=parameter_trans(     # transformation of parameter space
      log="Beta",                 # Beta > 0
      logit=c("Rho","Eta")        # 0 <= Rho, Eta <= 1
    ),
    paramnames=c("Beta","Rho","Eta")
  ) -> mf

plot(mf)

mf |>
  mif2(
    # keep other setting, update the decay
    cooling.fraction.50 = 0.3
  ) -> mf

plot(mf)

## ------- 2. mif in parallel ------- ##

## put transformation outside
## to avoid compilation within the foreach loop
measSIR |>
  pomp(
    partrans=parameter_trans(
      log="Beta",logit=c("Rho","Eta")
    ),
    paramnames=c("Beta","Rho","Eta")
  ) -> po

foreach(
  i=1:20,
  .combine=c,                         # concatenate all results
  .packages = "pomp",
  .options.RNG = 482947940            # set random seed
) %dorng% {
  po |>
    mif2(
      Np=2000, Nmif=50, cooling.fraction.50=0.5,
      rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02))
    )
} -> mifs_local

mifs_local |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1))) +
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")

## ------- 3. likelihood after mif ------- ##
foreach(
  mf = mifs_local,                 # loop over each object in `mifs_local`
  .combine = rbind,                # row bind all results
  .package = "pomp",
  .options.RNG = 900242057
) %dorng% {
  evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
  ll <- logmeanexp(evals,se=TRUE)
  mf |> coef() |> bind_rows() |>
    bind_cols(loglik=ll[1],loglik.se=ll[2])
} -> results

results |> filter(loglik==max(loglik))
