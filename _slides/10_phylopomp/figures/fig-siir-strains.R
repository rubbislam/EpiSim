## fig-siir-strains.R ------------------------------------------------------
## 10_phylopomp (Lecture C) -- a two-strain (SIIR) genealogy. With
## obscure=FALSE, phylopomp colours each lineage by the strain (deme) it
## belongs to, so you can see the two strains interleave in one tree.
##
## Run from the project root (where 10_phylopomp.Rproj lives):
##   Rscript figures/fig-siir-strains.R
## Writes: figures/fig-siir-strains.png
##
## phylopomp >= 0.19.4 required: simulator is runSIIR(); its initial states
## S_0, I1_0, I2_0, R_0 are FRACTIONS of `pop`.
## ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(phylopomp)
  library(ggplot2)
})

## Strain 2 spreads faster (Beta2 > Beta1); both are sampled equally.
set.seed(7)
GI <- runSIIR(
  time = 3,
  Beta1 = 4, Beta2 = 8,    # strain 2 is the faster spreader
  gamma = 1.5,             # recovery rate (shared)
  psi1 = 0.6, psi2 = 0.6,  # per-strain sampling rates
  pop = 300,
  S_0 = 0.95, I1_0 = 0.03, I2_0 = 0.02, R_0 = 0
)

nsample <- getInfo(GI, nsample = TRUE)$nsample

## obscure=FALSE keeps the strain labels, so plot() can colour each lineage
## by its strain (deme 1 vs deme 2). No points: the strain colours are the
## whole message here, and tip dots would only clutter them.
p <- plot(GI, obscure = FALSE) +
  labs(x = "time")

ggsave("figures/fig-siir-strains.png", p,
       width = 6.5, height = 4.0, dpi = 300, bg = "white")

cat("\n---- fig-siir-strains: verified numbers ----\n")
cat(sprintf("phylopomp version : %s\n", packageVersion("phylopomp")))
cat(sprintf("sampled tips      : %d\n", nsample))
cat("-------------------------------------------------\n")
