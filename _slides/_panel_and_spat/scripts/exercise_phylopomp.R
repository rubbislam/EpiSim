library(phylopomp)

## ---------- Exercise 1: simulation ----------- ##

# 1. simulate from a two-strain SIR model
??siir              # check the model setting

phylopomp::simulate(
  "SIIR",            # model name
  time = 5,          # simulation interval t=0 to t=5
  Beta1=2,           # transmission rate of strain 1
  Beta2=50,          # transmission rate of strain 2
  gamma=1,           # recovery rate of the infected
  psi1=2,            # sampling rate for individual infected by strain 1
  psi2=1,            # sampling rate for individual infected by strain 2
  S0=300,I1_0=20,I2_0=2,
) -> model.siir

model.siir |> plot(obscure=FALSE)
model.siir |> lineages(obscure=FALSE) |> plot()

# 2. simulate from a two-class SIR model with super-spreader
??si2r

phylopomp::simulate(
  "SI2R",            # model name
  time = 5,          # simulation interval
  Beta=5,            # transmission rate
  mu=2,              # mean superspreading-event cluster size
  gamma=1,
  psi1=1,
  psi2=0,
  sigma12=1,         # transition rate from normal-spreader to super-spreader
  sigma21=3,         # transition rate from normal-spreader to super-spreader
  S0=500,I0=10
) -> model.si2r

model.si2r |> plot(obscure=FALSE)
model.si2r |> lineages(obscure=FALSE) |> plot()

## ---------- Exercise 2: inference ----------- ##
library(pomp)

# 1. fit an SIR model to a tree simulated from an SEIR model
??seirs

phylopomp::simulate(
    "SEIRS",
    time = 5,
    Beta=4,sigma=1,gamma=1,psi=1,omega=1,
    S0=100,E0=3,I0=5,R0=100,
) -> G.seir
plot(G.seir)

# tune the parameters to see what happens
G.seir |>
  sir_pomp(
    Beta=4,gamma=2,psi=1,omega=1,
    S0=100,I0=8,R0=100
  ) -> po

po |> 
  pfilter(Np=1000) |> 
  replicate(n=20) |> 
  concat() -> pf

pf |> logLik() |> logmeanexp(se=TRUE)

# 2. fit an SEIR model to a tree simulated from an SIR model
??sir

phylopomp::simulate(
  "SIR",time=10,
  Beta=3,gamma=1,psi=2,omega=1,S0=100,I0=5
) -> G.sir
plot(G.sir)

## This can take a very long time!!!
G.sir |>
  seirs_pomp(
    Beta=3,sigma=1,gamma=1,psi=2,omega=1,
    S0=100,E0=0,I0=5,R0=0
  ) |> pfilter(Np=1000) |> 
  replicate(n=20) |> 
  concat() -> pf

pf |> logLik() |> logmeanexp(se=TRUE)



