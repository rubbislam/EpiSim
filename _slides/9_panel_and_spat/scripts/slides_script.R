## ----load-packages, echo=FALSE---------------------------------------------------------------------------------------
library(tidyverse)
library(pomp)
library(panelPomp)
library(phylopomp)
library(ggplot2)
library(reshape2)


## ----ppompBuild, echo=TRUE, eval=FALSE-------------------------------------------------------------------------------
library(panelPomp)

mod1 <- pomp(..., params = c('p1' = 0.1, 'p2' = 1.2, 'p3' = 0.9))
mod2 <- pomp(..., params = c('p1' = 0.1, 'p2' = 1.1, 'p3' = 0.6))
mod3 <- pomp(..., params = c('p1' = 0.1, 'p2' = 1.5, 'p3' = 0.75))

ppomp <- panelPomp(
  object = list(mod1, mod2, mod3),
  shared = c("p1" = 0.1), # Shared-value parameters
  specific = c("p2", "p3") # Unit-Specific parameters
)


## --------------------------------------------------------------------------------------------------------------------
measSIR <- panelMeasles(
  units = c("Consett", "London", "Hastings"),
  first_year = 1948,
  last_year = 1948
  )


## ----plotMeas, eval=FALSE, echo=TRUE---------------------------------------------------------------------------------
plot(measSIR)


## ----plotMeasEval, fig.height=3, echo=FALSE, fig.width=7-------------------------------------------------------------
df <- sapply(measSIR@unit_objects, function(x) x@data) |> data.frame()
df$time <- measSIR@unit_objects$Consett@times |> lubridate::date_decimal()
df_long <- tidyr::pivot_longer(
  data = df, cols = -time, names_to = 'units', values_to = 'cases'
)
ggplot(df_long, aes(x = time, y = cases)) +
  geom_line() + 
  facet_wrap(~units, nrow = 1, scales = 'free') + 
  theme_classic() + 
  scale_x_datetime(date_labels = '%b %Y') + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0), axis.title.x = element_blank())


## --------------------------------------------------------------------------------------------------------------------
coef(measSIR)


## --------------------------------------------------------------------------------------------------------------------
shared(measSIR)
shared(measSIR) <- c(shared(measSIR), 'alpha' = 1)
shared(measSIR)


## ----ppompPfilter, fig.height=3--------------------------------------------------------------------------------------
pfilter(measSIR, Np = 1000) |> plot(unit = 'Hastings')


## ----ppompMif2, eval=TRUE--------------------------------------------------------------------------------------------
mif2Out <- mif2(
  measSIR, rw.sd = rw_sd(rho = 0.02, R_0 = 0.02), 
  Nmif = 10, Np = 200, cooling.fraction.50 = 0.5, 
  block = TRUE  # block = TRUE does MPIF, usually best + faster.
  )


## --------------------------------------------------------------------------------------------------------------------
unitLogLik(mif2Out)
sum(unitLogLik(mif2Out))
logLik(mif2Out)


## ----calcPfilter, cache=TRUE-----------------------------------------------------------------------------------------
library(foreach)
library(doParallel)
library(doRNG)
registerDoParallel(detectCores()-1)
foreach(rep=1:10, .combine = c) %dopar% {
  pfilter(measSIR, Np = 5000)
} -> pfOut  # Object of length 10


## ----logmeanexp------------------------------------------------------------------------------------------------------
panel_logmeanexp(
  sapply(pfOut, unitLogLik), MARGIN = 1, se = TRUE
)


## ----haitiDepartments, fig.height=3.5, echo=FALSE--------------------------------------------------------------------
library(haitipkg)
library(tidyverse)

dep_plot_df <- haitiCholera %>%
  select(-report) %>%
  pivot_longer(
    data = .,
    cols = -date_saturday,
    values_to = 'cases',
    names_to = 'dep'
  ) %>%
  mutate(
    date = as.Date(date_saturday),
    dep = gsub("\\.", "_", dep)
  ) %>%
  mutate(
    dep = case_when(dep == "Grand_Anse" ~ "Grande_Anse", TRUE ~ dep)
  )

dep_labeller <- as_labeller(
  c(
    'Artibonite' = 'Artibonite',
    'Sud_Est' = 'Sud-Est',
    'Sud.Est' = 'Sud-Est',
    'Nippes' = 'Nippes',
    'Nord_Est' = 'Nord-Est',
    'Nord.Est' = 'Nord-Est',
    'Ouest' = 'Ouest',
    'Centre' = 'Centre',
    'Nord' = 'Nord',
    'Sud' = 'Sud',
    'Nord_Ouest' = 'Nord-Ouest',
    'Nord.Ouest' = 'Nord-Ouest',
    'Grande_Anse' = 'Grand\'Anse',
    'Grand.Anse' = 'Grand\'Anse'
  )
)

# Vector used to arrange figures in the same order that was used in Lee et al (2020)
plot_order <- c(
  'Artibonite',
  'Sud_Est',
  'Nippes',
  'Nord_Est',
  'Ouest',
  'Centre',
  'Nord',
  'Sud',
  'Nord_Ouest',
  'Grande_Anse'
)

# Plot reported cholera cases, by department. Creates Figure 1.
ggplot(dep_plot_df, aes(x = date, y = cases + 1)) +
  facet_wrap(~dep, nrow = 2, labeller = dep_labeller) +
  geom_line() +
  theme(axis.title.x = element_blank()) +
  ylab('Reported Cases') +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_date(
    date_labels = "'%y",
    breaks = seq.Date(
      from = as.Date("2011-01-01"),
      to = as.Date("2019-01-01"),
      by = '2 years'
    )
  ) + 
  theme_bw() + 
  theme(axis.title.x = element_blank())


## --------------------------------------------------------------------------------------------------------------------
devtools::install_github('jeswheel/haitipkg')
sp_mod <- haitipkg::haiti3_spatPomp()

## ----evalBPF, cache = TRUE, echo=FALSE-------------------------------------------------------------------------------
tic <- Sys.time()
bpf_out <- bpfilter(sp_mod, Np = 1000, block_size = 1)
toc <- Sys.time()
bpf_time <- toc - tic


## ----IBPFexample, eval=FALSE-----------------------------------------------------------------------------------------
# ## NOT RUN
# ibpf(
#   h3_spat,
#   Nbpf = 100,  # Equivalent to Nmif
#   Np = 2000,
#   sharedParNames = <character vector of shared parameter names>,
#   unitParNames = <character vector of unit-specific parameter names>,
#   spat_regression = 0.1,
#   rw.sd = ...,
#   cooling.fraction.50 = ...,
#   block_size = 1
# )


## ----install, eval=FALSE---------------------------------------------------------------------------------------------
library(devtools)
install_github("kingaa/phylopomp@0.14.8.0")


## ----install-check---------------------------------------------------------------------------------------------------
packageVersion("phylopomp")
library(phylopomp)


## ----sim-sir-2, fig.width=6, fig.height=4, dpi=300, echo=FALSE-------------------------------------------------------
library(phylopomp)
set.seed(1234)
phylopomp::simulate(
  "SIR",time=2,
  Beta=2,gamma=1,psi=2,
  S0=1000,I0=5
)|>
  phylopomp::simulate(
    Beta=5,gamma=2,time=5,psi=3
  ) -> model.sir # update params
plot(model.sir)




## ----sim-seir1-2, fig.width=6, fig.height=4, dpi=300, echo=FALSE-----------------------------------------------------
set.seed(1234)
phylopomp::simulate(
  "SEIR",time=2,
  Beta=2,sigma=2,
  gamma=1,psi=2,
  S0=1000,E0=0,I0=5
) |>
  phylopomp::simulate(
    time=5,
    Beta=5,gamma=2,
    psi=3
  ) -> model.seir
plot(model.seir)




## ----sim-seir2-2, fig.width=6, fig.height=4, dpi=300,echo=FALSE------------------------------------------------------
plot(model.seir, obscure=FALSE)       # show E and I transitions




## ----sim-seir3-2, fig.width=6, fig.height=4, dpi=300, echo=FALSE-----------------------------------------------------
model.seir |>
  lineages( 
    obscure=FALSE
  ) |> 
  plot()


## ----infer-sir-1, eval=FALSE-----------------------------------------------------------------------------------------
set.seed(1234)
# simulate a geneal object from an SIR model
phylopomp::simulate("SIR",time=10,
  Beta=3,gamma=1,psi=2,omega=1,S0=100,I0=5) -> x
x |>        # build a pomp object
  sir_pomp(
    Beta=3,gamma=1,psi=2,omega=1,
    S0=100,I0=5,R0=0
  ) -> po
po |> pfilter(Np=5000) -> pf
pf |> logLik()




## ----infer-seir-1, eval=FALSE----------------------------------------------------------------------------------------
set.seed(1234)
phylopomp::simulate("SEIRS",      # simulate a genealogy from an SEIR model
  Beta=4,sigma=1,gamma=1,psi=1,omega=1,
  S0=100,E0=3,I0=5,R0=100, time=5
) -> G
G |>
  seirs_pomp(
    Beta=4,sigma=1,gamma=1,psi=1,omega=1,
    S0=100,E0=3,I0=5,R0=100
  ) |> pfilter(Np=1000) |>
  replicate(n=20) |> concat() -> pf
pf |> logLik() |> logmeanexp(se=TRUE)

