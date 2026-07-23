## =============================================================================
##  Lecture 4: Iterated filtering (mif2) -- SOLUTIONS
##
##  Try scripts/EXERCISES.R first. Open 4_mle_mif2.Rproj before running this.
##  Works on macOS, Linux and Windows.
## =============================================================================


# Setup -------------------------------------------------------------------

library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)
library(reshape2)   # for melt()

source("scripts/model_measSIR.R")

cl <- makePSOCKcluster(max(1, detectCores() - 1))
registerDoParallel(cl)

measSIR |>
  pomp(
    partrans=parameter_trans(log="Beta", logit=c("Rho","Eta")),
    paramnames=c("Beta","Rho","Eta")
  ) -> po


# Exercise 1: your first mif2 run -----------------------------------------

set.seed(42)
po |>
  mif2(
    Np=2000, Nmif=50, cooling.fraction.50=0.5,
    rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02))
  ) -> mf

plot(mf)

## The loglik panel climbs steeply and then flattens -- the search is finding
## much better parameters than our starting guess. The parameter panels are
## still wandering at iteration 50: the search has found a good LIKELIHOOD
## but has not pinned down individual parameters. Hold that thought.


# Exercise 2: turn the knobs ----------------------------------------------

mf |> mif2(cooling.fraction.50=0.3) -> mf2
plot(mf2)

## (a) cooling.fraction.50 small (0.1): perturbations die quickly, the search
##     freezes early -- fast but it can get stuck wherever it happens to be.
##     Large (0.9): stays noisy for a long time, explores more, converges slowly.
## (b) rw.sd too small (0.002): the swarm barely moves; traces look like flat
##     lines and the loglik improves painfully slowly.
##     rw.sd too large (0.2): traces look like noise; the search never settles.
## (c) Np too small (200): the likelihood estimate is so noisy that the search
##     is partly chasing Monte Carlo error rather than real signal.
##
## cooling.fraction.50=0.5 is a promise about ITERATION 50: at iteration 50 the
## random-walk sd is half its starting value. The decay is geometric and
## continuous -- it keeps shrinking after 50 too.


# Exercise 3: many searches at once ---------------------------------------

foreach(
  i=1:20, .combine=c, .packages="pomp", .export="po", .options.RNG=482947940
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
  ggplot(aes(x=iteration, y=value, group=L1)) +
  geom_line(alpha=.3, color="#0072B2") +
  facet_wrap(~name, scales="free_y") +
  labs(x="mif2 iteration", y=NULL) +
  theme_bw()

## THE KEY OBSERVATION, and it is the punchline of both lectures:
##
##   The loglik panel converges tightly. The parameter panels do NOT.
##
## Measured over these 20 searches (your numbers will be close):
##
##   loglik : -106.7 to -105.8   -- coefficient of variation 0.24%
##   Beta   :  15.6  to  23.4    -- CV 9.6%
##   Rho    :   0.24 to   0.37   -- CV 10.7%
##   Eta    :   0.044 to  0.067  -- CV 10.3%
##
## The fit converges about 40x tighter than the parameters do. All 20 searches
## agree the epidemic is well described; they disagree by 50% about Beta.
##
## Why? Look at how the endpoints correlate with each other:
##
##   cor(Beta, Rho) = +0.96      cor(Beta, Eta) = -0.97      cor(Rho, Eta) = -0.98
##
## Those are nearly +/-1. The 20 answers are not scattered in a cloud -- they
## lie almost exactly on a LINE through (Beta, Rho, Eta) space. Slide along it
## and the likelihood barely changes: a higher transmission rate, spread through
## a smaller pool of initial susceptibles, with a higher fraction of cases
## reported, produces very nearly the same 42 numbers.
##
## This is a RIDGE -- the same phenomenon you met in Lecture A, where the
## Beta-Gamma surface was a long diagonal crest. (It is a different ridge:
## here Gamma is held fixed at 0.5 and is not being estimated at all.)
##
## And you can find WHAT the data actually pin down. Try plotting Beta*Eta:
##
##   Beta    coefficient of variation: 9.6%
##   Beta*Eta coefficient of variation: 2.1%     <- 4.5x tighter
##
## Beta*Eta is (roughly) the initial force of infection: the transmission rate
## times the fraction of the population that is susceptible to begin with.
## THAT is what these 42 numbers measure. Beta on its own, they do not.
##
## This is the single most useful habit in the whole module: when replicates
## disagree about a parameter, go looking for the COMBINATION they agree on.
## It is usually something with a clean epidemiological meaning.
##
## The practical lesson: a single mif2 run tells you almost nothing. Run many,
## and read the loglik panel and the parameter panels as answering DIFFERENT
## questions -- "is the fit good?" and "do I know this parameter?"


# Exercise 4: the number mif2 prints is NOT your likelihood ----------------

foreach(
  mf=mifs_local, .combine=rbind, .packages=c("pomp","dplyr"),
  .options.RNG=900242057
) %dorng% {
  evals <- replicate(10, logLik(pfilter(mf, Np=5000)))
  ll <- logmeanexp(evals, se=TRUE)
  mf |> coef() |> bind_rows() |> bind_cols(loglik=ll[1], loglik.se=ll[2])
} -> results

results |> filter(loglik==max(loglik))

## Compare against what mif2 reported:
sapply(mifs_local, logLik)

## The honest re-estimate is HIGHER than mif2's own number, consistently --
## across these 20 replicates it was higher in 20 out of 20, by about +0.6 log
## units on average. Two reasons stack up:
##
##   1. During mif2 the parameters are deliberately perturbed, so the
##      likelihood being reported belongs to a NOISIER model than yours.
##      Extra noise costs likelihood.
##   2. Look at traces(mf): it has Nmif+1 rows and the last loglik is NA.
##      mif2 records the likelihood of the parameters at the START of each
##      iteration, so your final parameter vector -- the one in coef(mf),
##      the one you keep -- has never been filtered at all.
##
## So logLik(mf) is the likelihood of the previous, perturbed parameters.
## It is not the likelihood of your answer. Always re-filter.
##
## For reference: the best loglik here lands around -105.5, which beats the
## -107.5 best point from the coarse Beta-Gamma grid in Lecture A. The search
## found something the grid missed.


# Exercise 5 (stretch): does the starting point matter? -------------------

runif_design(
  lower=c(Beta=5,  Rho=0.2, Eta=0.01),
  upper=c(Beta=40, Rho=0.9, Eta=0.20),
  nseq=12
) -> starts
starts$Gamma <- 0.5; starts$k <- 10; starts$N <- 38000

head(starts)

summary(starts)

plot(starts)

foreach(
  i=1:nrow(starts), .combine=c, .packages="pomp", .export=c("po","starts"), .options.RNG=531
) %dorng% {
  po |>
    mif2(
      params=unlist(starts[i,]),
      Np=2000, Nmif=50, cooling.fraction.50=0.5,
      rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02))
    )
} -> mifs_global

mifs_global |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration, y=value, group=L1)) +
  geom_line(alpha=.3, color="#0072B2") +
  facet_wrap(~name, scales="free_y") +
  labs(x="mif2 iteration", y=NULL, title="Searches from 12 scattered starting points") +
  theme_bw()

## Measured result from 12 scattered starts at these settings:
##
##   11 of 12 climbed to a loglik between -108.0 and -106.3.
##    1 of 12 got stuck at -140 and never recovered.
##
## So: MOST searches find the good region from a cold start, but not all. The
## stuck one is not a bug -- it is the normal behaviour of a search on a hard
## surface, and it is exactly why you run many starts and keep the best.
##
## Two warnings from experience, both of which you will meet on real problems:
##
##   * Turn the effort down (try Np=1000, Nmif=30) and this gets much worse --
##     in a test, only 4 of 12 reached the good region, several stalled between
##     -160 and -410, and one returned -Inf. A search that "converges" is
##     partly a statement about how much compute you gave it.
##   * A loglik of -Inf means some particle filter hit an impossible
##     observation and every particle died. Usually the starting parameters
##     were absurd. Don't panic; look at the starting values.
##
## Among the 11 successful searches, Beta ran from 6.1 to 23.1 -- a factor of
## nearly four -- for essentially the same likelihood. The ridge from
## Exercise 3 is still there, and starting from far away just spreads you
## further along it.


# Exercise 6 (lunch break): does the latent period earn its keep? ----------

## The SEIR model is the SIR with an extra "exposed" (E) compartment: infected
## individuals spend a latent period before becoming infectious. You built and
## simulated it in Lecture 2; here we FIT it and compare to the SIR.
source("scripts/model_measSEIR.R")   # -> measSEIR (Sigma = 1, a 1-week latent period)

## Sigma is a positive rate, so it is log-transformed alongside Beta:
measSEIR |> pomp(
  partrans = parameter_trans(log = c("Beta","Sigma"), logit = c("Rho","Eta")),
  paramnames = c("Beta","Sigma","Rho","Eta")
) -> po_seir

## 20 local searches, estimating Sigma as well as Beta, Rho, Eta:
foreach(i=1:20, .combine=c, .packages="pomp", .export="po_seir",
        .options.RNG=482947940) %dorng% {
  po_seir |> mif2(Np=2000, Nmif=50, cooling.fraction.50=0.5,
                  rw.sd=rw_sd(Beta=0.02, Sigma=0.02, Rho=0.02, Eta=ivp(0.02)))
} -> mifs_seir

## Honest re-estimate of each endpoint, then keep the best:
foreach(mf=mifs_seir, .combine=rbind, .packages=c("pomp","dplyr"),
        .options.RNG=900242057) %dorng% {
  evals <- replicate(10, logLik(pfilter(mf, Np=5000)))
  ll <- logmeanexp(evals, se=TRUE)
  mf |> coef() |> bind_rows() |> bind_cols(loglik=ll[1], loglik.se=ll[2])
} -> seir_results

seir_results |> filter(loglik == max(loglik))    # about -104.3
sort(seir_results$Sigma)                          # ranges ~1.2 to 3.4

## ANSWER. The best SEIR fit lands at loglik = -104.3 (se 0.02), versus -105.4
## for the SIR: the latent period buys about 1.1 log-likelihood units.
##
## Is that worth an extra parameter? Barely. The AIC rule of thumb wants MORE
## than ~2 units per added parameter, and 1.1 does not clear it:
##   AIC_SIR  = -2(-105.4) + 2*3 = 216.8   (Beta, Rho, Eta estimated)
##   AIC_SEIR = -2(-104.3) + 2*4 = 216.5   (Beta, Rho, Eta, Sigma estimated)
## a difference of about 0.3 -- essentially a tie. A likelihood-ratio test
## agrees: 2*1.1 = 2.2, below the chi-square(1) cutoff of 3.84.
##
## The deeper reason is identifiability. Across the 20 searches Sigma ran from
## about 1.2 to 3.4 -- a factor of three -- for essentially the same
## likelihood. The data barely constrain the latent rate; Sigma is another
## ridge direction, just like Beta in Exercise 3. A single outbreak like
## Consett's cannot "see" the latent period well enough to prefer the SEIR.
##
## The honest lesson: a richer model is not automatically a better one. Adding
## structure the data cannot resolve buys almost nothing -- here, roughly one
## log unit, spent on a parameter you cannot even pin down.


# Cleanup -----------------------------------------------------------------

stopCluster(cl)
