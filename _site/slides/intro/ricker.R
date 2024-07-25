params <-
list(prefix = "ricker")

## ----prelims,include=FALSE,purl=TRUE,cache=FALSE-------------------------------------------------------------------
library(pomp)


## ----parus-data----------------------------------------------------------------------------------------------------
dat <- read.csv("parus.csv")
head(dat)
plot(pop~year,data=dat,type='o')


## ----parus-pomp1---------------------------------------------------------------------------------------------------
library(pomp)
parus <- pomp(dat,times="year",t0=1959)


## ----parus-plot1---------------------------------------------------------------------------------------------------
plot(parus)


## ----parus-sim-defn------------------------------------------------------------------------------------------------
stochStep <- Csnippet("
  e = rnorm(0,sigma);
  N = r*N*exp(-c*N+e);
")
pomp(parus,rprocess=discrete_time(step.fun=stochStep,delta.t=1),
     paramnames=c("r","c","sigma"),statenames=c("N","e")) -> parus


## ----ricker-first-sim----------------------------------------------------------------------------------------------
sim <- simulate(parus,params=c(N_0=1,e_0=0,r=12,c=1,sigma=0.5),
                format="data.frame")
plot(N~year,data=sim,type='o')


## ----parus-rmeas-defn----------------------------------------------------------------------------------------------
rmeas <- Csnippet("pop = rpois(phi*N);")


## ----parus-dmeas-defn----------------------------------------------------------------------------------------------
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")


## ----parus-add-meas------------------------------------------------------------------------------------------------
pomp(parus,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> parus


## ----ricker-add-params---------------------------------------------------------------------------------------------
coef(parus) <- c(N_0=1,e_0=0,r=20,c=1,sigma=0.1,phi=200)

## ----ricker-second-sim,results='markup'----------------------------------------------------------------------------
library(ggplot2)
sims <- simulate(parus,nsim=3,format="d",include.data=TRUE)
ggplot(data=sims,mapping=aes(x=year,y=pop))+geom_line()+
  facet_wrap(~.id,ncol=1,scales="free_y")


## ----parus-skel-defn-----------------------------------------------------------------------------------------------
skel <- Csnippet("DN = r*N*exp(-c*N);")


## ----parus-add-skel------------------------------------------------------------------------------------------------
parus <- pomp(parus,skeleton=map(skel),paramnames=c("r","c"),statenames=c("N"))


## ----parus-first-traj,results='markup'-----------------------------------------------------------------------------
traj <- trajectory(parus,params=c(N_0=1,r=12,c=1),format="data.frame")
plot(N~year,data=sim,type='o')
lines(N~year,data=traj,type='l',col='red')


## ----plot-ricker---------------------------------------------------------------------------------------------------
plot(parus)


## ----sim-ricker1---------------------------------------------------------------------------------------------------
x <- simulate(parus)


## ------------------------------------------------------------------------------------------------------------------
class(x)
plot(x)


## ------------------------------------------------------------------------------------------------------------------
y <- as.data.frame(parus)
head(y)
head(simulate(parus,format="data.frame"))


## ------------------------------------------------------------------------------------------------------------------
x <- simulate(parus,nsim=10)
class(x)
sapply(x,class)
x <- simulate(parus,nsim=10,format="data.frame")
head(x)
str(x)


## ----fig.height=8--------------------------------------------------------------------------------------------------
library(ggplot2)
x <- simulate(parus,nsim=9,format="data.frame",include.data=TRUE)
ggplot(data=x,aes(x=year,y=pop,group=.id,color=(.id=="data")))+
  geom_line()+guides(color="none")+
  facet_wrap(~.id,ncol=2)


## ----traj-ricker---------------------------------------------------------------------------------------------------
y <- trajectory(parus,format="data.frame")
plot(N~year,data=y,type="l")


## ----coef-ricker---------------------------------------------------------------------------------------------------
coef(parus)


## ------------------------------------------------------------------------------------------------------------------
theta <- coef(parus)
theta[c("r","N_0")] <- c(5,3)
y <- trajectory(parus,params=theta)
plot(N~year,data=as.data.frame(y),type='l')
x <- simulate(parus,params=theta)
plot(x,var="pop")


## ------------------------------------------------------------------------------------------------------------------
coef(parus,c("r","N_0","sigma")) <- c(39,0.5,1)
coef(parus)
plot(simulate(parus),var="pop")

