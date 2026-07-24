## fig-sir-genealogy.R -----------------------------------------------------
## 10_phylopomp (Lecture C) -- one simulated SIR genealogy, shown two ways:
##   1. the tree itself         -> fig-sir-genealogy.png
##   2. its lineages-through-time -> fig-sir-ltt.png
##
## Both come from the SAME simulated genealogy (one seed), so they always
## agree. Run from the project root (where 10_phylopomp.Rproj lives):
##   Rscript figures/fig-sir-genealogy.R
## Writes: figures/fig-sir-genealogy.png, figures/fig-sir-ltt.png
##
## phylopomp >= 0.19.4 required: the simulator is runSIR() (the old
## simulate("SIR", ...) string form was removed), and initial states are
## FRACTIONS of `pop`, not counts.
## ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(phylopomp)
  library(ggplot2)
})

## One modest, legible outbreak: ~150 sampled tips.
set.seed(2024)
G <- runSIR(
  time = 4,        # observation window, t = 0 .. 4
  Beta = 3,        # transmission rate
  gamma = 1,       # recovery rate
  psi = 0.4,       # per-lineage sampling rate (sets how many tips)
  pop = 500,       # total host population
  S0 = 0.97,       # initial fractions (sum to 1), NOT counts
  I0 = 0.03,
  R0 = 0
)

nsample <- getInfo(G, nsample = TRUE)$nsample

## 1. the genealogy -------------------------------------------------------
## points=TRUE draws the sampling events (tips) and coalescences (nodes).
p_tree <- plot(G, points = TRUE) +
  labs(x = "time")

ggsave("figures/fig-sir-genealogy.png", p_tree,
       width = 6.5, height = 4.0, dpi = 300, bg = "white")

## 2. lineages through time ----------------------------------------------
## How many lineages are "alive" in the tree at each moment: it rises as the
## epidemic grows and falls as it burns out -- the tree's shape carries the
## dynamics, which is exactly what phylodynamics exploits.
p_ltt <- plot(lineages(G)) +
  labs(x = "time", y = "lineages")

ggsave("figures/fig-sir-ltt.png", p_ltt,
       width = 6.5, height = 3.0, dpi = 300, bg = "white")

## the numbers this figure claims ----------------------------------------
cat("\n---- fig-sir-genealogy: verified numbers ----\n")
cat(sprintf("phylopomp version : %s\n", packageVersion("phylopomp")))
cat(sprintf("sampled tips      : %d\n", nsample))
cat("-------------------------------------------------\n")
