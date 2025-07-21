# CODE FOR LECTURE 6 SLIDES 
# 
# WARNING: This code takes several hours to run on a machine with 24 cores. 
#          you may not be able to run this code on your own computer in a 
#          reasonable time frame.

library(foreach)
library(doParallel)
library(doRNG)
library(iterators)
library(tidyverse)
library(pomp)
set.seed(1350254336)


# Applying IF2 to the Consett measles outbreak ----------------------------

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

plot(measSIR)


# Setting up Estimation Problem -------------------------------------------

fixed_params <- c(N=38000, Gamma=2, k=10)
coef(measSIR,names(fixed_params)) <- fixed_params
coef(measSIR)

# Sanity Check ------------------------------------------------------------

measSIR |>
  simulate(nsim=20,format="data.frame",include.data=TRUE) |>
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data")) +
  geom_line() + guides(color="none")

# Parallel Computing ------------------------------------------------------

library(foreach)
library(doParallel)
library(doRNG)

cores <- parallel::detectCores()
# registerDoParallel(cores)

# Windows: 
cl <- makePSOCKcluster(cores)  # Windows
registerDoParallel(cl)  # Windows


# Running a particle filter I ----------------------------------------------

registerDoRNG(123)
tic <- Sys.time()
foreach(i=1:10,.combine=c, .packages = 'pomp') %dopar% {
  measSIR |> pfilter(Np=5000)
} -> pf
pf |> logLik() |> logmeanexp(se=TRUE) -> L_pf
L_pf
toc <- Sys.time()


# Building up a picture of the likelihood surface I -----------------------

pf[[1]] |>  # 10 pfilterd objects, all have same parameters
  coef() |>  # Extract parameters 
  bind_rows() |>  # Convert to data.frame
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2]) |>  # Add evaluation.
  write_csv("measles_params.csv")  # Save to "database"


# A local search of the likelihood surface II -----------------------------

registerDoRNG(829479)
foreach(i=1:20,.combine=c,.packages = 'pomp') %dopar% {
  measSIR |>
    mif2(
      Np=2000, Nmif=50, cooling.fraction.50=0.5,
      rw.sd=rw_sd(Beta=0.02, Rho=0.02, Eta=ivp(0.02)),
      partrans=parameter_trans(
        log="Beta",logit=c("Rho","Eta")
      ),
      paramnames=c("Beta","Rho","Eta")
    )
} -> mifs_local


# Iterated Filtering Diagnostics ------------------------------------------

mifs_local |>
  traces() |>
  melt() |>
  ggplot(aes(x=iteration,y=value,group=.L1,color=factor(.L1)))+
  geom_line()+
  guides(color="none")+
  facet_wrap(~name,scales="free_y")


# Estimating the likelihood II --------------------------------------------

registerDoRNG(900242)
foreach(mf=mifs_local,.combine=rbind,.packages='pomp') %dopar% {
  evals <- replicate(10, logLik(pfilter(mf,Np=5000)))
  ll <- logmeanexp(evals,se=TRUE)
  mf |> coef() |> bind_rows() |>
    bind_cols(loglik=ll[1],loglik.se=ll[2])
} -> results

# Estimating the likelihood III -------------------------------------------

results |> filter(loglik==max(loglik))

# Estimating the likelihood IV --------------------------------------------

pairs(~loglik+Beta+Eta+Rho,data=results,pch=16)

# Building up a picture of the likelihood surface -------------------------

read_csv("measles_params.csv") |>
  bind_rows(results) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


# Global search of the likelihood surface III -----------------------------

set.seed(2062379496)

runif_design(
  lower=c(Beta=5,Rho=0.2,Eta=0),
  upper=c(Beta=80,Rho=0.9,Eta=1),
  nseq=400
) -> guesses

mf1 <- mifs_local[[1]]


# Global search of the likelihood surface IV ------------------------------

registerDoRNG(127040)
foreach(guess=iter(guesses,"row"), .combine=rbind, packages='pomp') %dopar% {
  mf1 |>
    mif2(params=c(guess,fixed_params)) |>
    mif2(Nmif=100) -> mf
  replicate(
    10,
    mf |> pfilter(Np=5000) |> logLik()
  ) |>
    logmeanexp(se=TRUE) -> ll
  mf |> coef() |> bind_rows() |>
    bind_cols(loglik=ll[1],loglik.se=ll[2])
} -> results

# Save results 
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


# Show our maximum likelihood result (so far):
results |> filter(loglik==max(loglik))


# Global search for the likelihood surface VII -----------------------------

read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-50) |>
  bind_rows(guesses) |>
  mutate(type=if_else(is.na(loglik),"guess","result")) |>
  arrange(type) -> all

pairs(~loglik+Beta+Eta+Rho, data=all, pch=16, cex=0.3,
      col=ifelse(all$type=="guess",grey(0.5),"red"))


# Global search for the likelihood surface X ------------------------------

all |>
  filter(type=="result") |>
  filter(loglik>max(loglik)-10) |>
  ggplot(aes(x=Eta,y=loglik))+
  geom_point()+
  labs(
    x=latex2exp::TeX('$\\eta$'),
    title="poor man's profile likelihood"
  )


# Profile likelihood over \eta II -----------------------------------------

# Create profile search box
read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-20,loglik.se<2) |>
  sapply(range) -> box
box


# Profile likelihood over \eta III ----------------------------------------

# Randomly sample from search box.
# Total number of starting points
#    - 40 unique values of Eta (Eta = ..., length = 40)
#    - 15 random values of all other parameters, for each fixed ETA (nprof = 15)
#    - Total: 40 X 15 = 600 unique starting points.
freeze(seed=1196696958,
       profile_design(
         Eta=seq(0.01,0.95,length=40),
         lower=box[1,c("Beta","Rho")],
         upper=box[2,c("Beta","Rho")],
         nprof=15, type="runif"
       )) -> guesses
plot(guesses)


# Profile over \eta VI ----------------------------------------------------

registerDoRNG(830007)
foreach(guess=iter(guesses,"row"), .combine=rbind, packages='pomp') %dopar% {
  mf1 |>
    mif2(params=c(guess,fixed_params),
         rw.sd=rw_sd(Beta=0.02,Rho=0.02)) |>
    mif2(Nmif=100,cooling.fraction.50=0.3) -> mf
  replicate(
    10,
    mf |> pfilter(Np=5000) |> logLik()) |>
    logmeanexp(se=TRUE) -> ll
  mf |> coef() |> bind_rows() |>
    bind_cols(loglik=ll[1],loglik.se=ll[2])
} -> results


# Visualizing the profile likelihood --------------------------------------

# Add results to database 
read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")


read_csv("measles_params.csv") |>
  filter(loglik>max(loglik)-10) -> all

pairs(~loglik+Beta+Eta+Rho,data=all,pch=16)


# Visualizing the profile likelihood II -----------------------------------

results |>
  ggplot(aes(x=Eta,y=loglik))+
  geom_point() + xlab(expression(eta))


# Visualizing the profile likelihood III ----------------------------------

results |>
  filter(is.finite(loglik)) |>
  group_by(round(Eta,5)) |>  # Group by rounded values of Eta
  filter(rank(-loglik)<3) |>  # Only take top two results for each Eta value
  ungroup() |>
  filter(loglik>max(loglik)-20) |>  # Only look at top of likelihood results
  ggplot(aes(x=Eta,y=loglik))+
  geom_point() + xlab(expression(eta))

# Visualizing the profile likelihood VIII ---------------------------------

maxloglik <- max(results$loglik,na.rm=TRUE)
ci.cutoff <- maxloglik-0.5*qchisq(df=1,p=0.95)

results |>
  filter(is.finite(loglik)) |>
  group_by(round(Eta,5)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  ggplot(aes(x=Eta,y=loglik))+
  geom_point() + xlab(expression(eta)) +
  geom_smooth(method="loess",span=0.25)+
  geom_hline(color="red",yintercept=ci.cutoff)+
  lims(y=maxloglik-c(5,0))


# Visualizing the profile likelihood XI -----------------------------------

results |>
  filter(is.finite(loglik)) |>
  group_by(round(Eta,5)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  mutate(in_ci=loglik>max(loglik)-1.92) |>
  ggplot(aes(x=Eta,y=Rho,color=in_ci))+
  geom_point()+
  labs(
    color="inside 95% CI?",
    x=expression(eta),
    y=expression(rho),
    title="profile trace"
  )

# Profile over rho I ------------------------------------------------------

results |>
  filter(is.finite(loglik)) |>
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) |>
  summarize(min=min(Rho),max=max(Rho)) -> rho_ci

read_csv("measles_params.csv") |>
  group_by(cut=round(Rho,2)) |>
  filter(rank(-loglik)<=10) |>
  ungroup() |>
  arrange(-loglik) |>
  select(-cut,-loglik,-loglik.se) -> guesses


# Profile over rho II -----------------------------------------------------

registerDoRNG(210568)
foreach(guess=iter(guesses,"row"), .combine=rbind, packages='pomp') %dopar% {
  mf1 |>
    mif2(params=guess,
         rw.sd=rw_sd(Beta=0.02,Eta=ivp(0.02))) |>
    mif2(Nmif=100,cooling.fraction.50=0.3) |>
    mif2() -> mf
  replicate(
    10,
    mf |> pfilter(Np=5000) |> logLik()) |>
    logmeanexp(se=TRUE) -> ll
  mf |> coef() |> bind_rows() |>
    bind_cols(loglik=ll[1],loglik.se=ll[2])
} -> results


read_csv("measles_params.csv") |>
  bind_rows(results) |>
  filter(is.finite(loglik)) |>
  arrange(-loglik) |>
  write_csv("measles_params.csv")

# Profile over rho: results I + II -----------------------------------------


results |> filter(is.finite(loglik)) -> results

pairs(~loglik+Beta+Eta+Rho,data=results,pch=16)

results |>
  filter(loglik>max(loglik)-10,loglik.se<1) |>
  group_by(round(Rho,2)) |>
  filter(rank(-loglik)<3) |>
  ungroup() |>
  ggplot(aes(x=Rho,y=loglik))+
  geom_point() + xlab(expression(rho)) +
  geom_hline(
    color="red",
    yintercept=max(results$loglik)-0.5*qchisq(df=1,p=0.95)
  )


# Profile over rho: results III -------------------------------------------

results |>
  filter(loglik>max(loglik)-0.5*qchisq(df=1,p=0.95)) |>
  summarize(min=min(Rho),max=max(Rho)) -> rho_ci


# Visualizing Simulations I -----------------------------------------------

read_csv("measles_params.csv") |>
  filter(loglik == max(loglik)) |>
  select(-loglik, -loglik.se) -> best.params

measSIR |>
  simulate(
    params=unlist(best.params),
    nsim=1000, format="data.frame", include.data=TRUE
  ) -> sims


# Visualizing Simulations II ----------------------------------------------

sims |>
  mutate(data=.id=="data") |>
  group_by(week,data) |>
  reframe(
    p=c(0.025,0.5,0.975),
    value=wquant(reports,probs=p),
    name=c("lo","med","up")
  ) |>
  select(-p) |> pivot_wider() |> ungroup() |>
  ggplot(aes(x=week,y=med,color=data,fill=data,ymin=lo,ymax=up))+
  geom_line()+ geom_ribbon(alpha=0.2,color=NA) +
  labs(y="reports")+
  theme_bw() + guides(color="none",fill="none")


stopCluster(cl)
