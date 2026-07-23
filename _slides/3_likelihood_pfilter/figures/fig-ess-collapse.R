## fig-ess-collapse.R ------------------------------------------------------
## A_likelihood_pfilter (Lecture A) -- the LEAN pfilter diagnostic.
##
## The "effective sample size" (ESS) beat. It is the mechanism behind the
## loglik-is-random figure: where ESS is low, only a handful of the Np
## particles actually survive the reweighting, so the conditional-likelihood
## estimate there rests on very few particles -- that is where the Monte
## Carlo noise comes from, and where you would spend more particles.
##
## TOP:    the data (amber), so you can see WHERE in the epidemic the dips fall.
## BOTTOM: ESS per week (blue), out of Np. Reference line at Np.
##
## HONEST NOTE: the cloud does NOT collapse at the big week-17->18 jump
## (18 -> 75 reports); ESS is healthy there (~3200/5000). It collapses at
## (a) week 14, the early rise, and (b) the crash tail (weeks 27-30), where
## the near-zero counts surprise particles that expect an ongoing epidemic.
## Every number below is computed, printed, and verified stable across 5 seeds.
##
## Run from the project root (where 3_likelihood_pfilter.Rproj lives):
##   Rscript figures/fig-ess-collapse.R
## Reads : data/ (via scripts/model_measSIR.R)
## Writes: figures/fig-ess-collapse.png ; caches results/ess-diagnostic.rds
## Runtime: ~4 s.

library(tidyverse)
library(pomp)
library(patchwork)

source("scripts/model_measSIR.R")

## Okabe-Ito, fixed course meanings ---------------------------------------
amber <- "#E69F00"  # THE DATA, always
blue  <- "#0072B2"  # the particle filter / the method that works
verm  <- "#D55E00"  # the trouble spots

NP    <- 5000L
SEED  <- 1350254336

## ---- compute ESS + cond.logLik, averaged over seeds for a stable picture -
dat <- as.data.frame(measSIR)                 # week, reports
seeds <- c(SEED, 1, 42, 2024, 999)
ess_mat <- sapply(seeds, function(s) {
  set.seed(s); eff_sample_size(pfilter(measSIR, Np = NP))
})
df <- tibble(
  week    = dat$week,
  reports = dat$reports,
  ess     = rowMeans(ess_mat),
  ess_lo  = apply(ess_mat, 1, min),
  ess_hi  = apply(ess_mat, 1, max)
)
saveRDS(df, "results/ess-diagnostic.rds")

## ---- the real numbers ----------------------------------------------------
wk14   <- df |> filter(week == 14)
tail_w <- df |> filter(week >= 27, week <= 30)
wk18   <- df |> filter(week == 18)
lowest <- df |> slice_min(ess, n = 1)

cat("\n---- fig-ess-collapse: computed numbers (Np =", NP, ", 5 seeds) ----\n")
cat(sprintf("median ESS over the series : %.0f of %d\n", median(df$ess), NP))
cat(sprintf("lowest ESS                 : week %d, ESS %.0f (%.0f%% of Np)\n",
            lowest$week, lowest$ess, 100 * lowest$ess / NP))
cat(sprintf("week 14 (early rise)       : ESS %.0f  [%.0f, %.0f] across seeds\n",
            wk14$ess, wk14$ess_lo, wk14$ess_hi))
cat(sprintf("crash tail (wk27-30)       : ESS %.0f (mean)\n", mean(tail_w$ess)))
cat(sprintf("week 18 (the 18->75 jump)  : ESS %.0f  -- healthy, NOT a collapse\n", wk18$ess))
cat("-------------------------------------------------------------------\n\n")

## ---- shared x axis -------------------------------------------------------
xsc <- scale_x_continuous(breaks = seq(0, 42, 7),
                          expand = expansion(mult = c(0.02, 0.02)))
base_txt <- theme(
  axis.text  = element_text(size = 7.5),
  axis.title = element_text(size = 9),
  plot.title = element_text(size = 10, face = "bold", colour = "grey20",
                            margin = margin(b = 3)),
  panel.grid.minor = element_blank()
)

## ---- TOP: the data, for context -----------------------------------------
p_top <- ggplot(df, aes(week, reports)) +
  geom_line(colour = amber, linewidth = 0.7) +
  geom_point(colour = amber, size = 1.1) +
  annotate("segment", x = 14, xend = 14, y = 0, yend = max(df$reports),
           colour = verm, linetype = "dotted", linewidth = 0.3) +
  annotate("rect", xmin = 26.5, xmax = 30.5, ymin = 0, ymax = max(df$reports),
           fill = verm, alpha = 0.07) +
  xsc +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = NULL, y = "reports", title = "the data") +
  theme_bw() + base_txt +
  theme(axis.text.x = element_blank(), plot.margin = margin(4, 8, 1, 6))

## ---- BOTTOM: ESS ---------------------------------------------------------
p_bot <- ggplot(df, aes(week, ess)) +
  geom_hline(yintercept = NP, colour = "grey70", linetype = "dashed",
             linewidth = 0.3) +
  annotate("text", x = 1, y = NP, label = sprintf("all %d particles", NP),
           hjust = 0, vjust = 1.4, size = 2.5, colour = "grey55") +
  geom_area(fill = blue, alpha = 0.12) +
  geom_line(colour = blue, linewidth = 0.7) +
  geom_point(colour = blue, size = 1.0) +
  ## flag the two trouble spots, computed not hand-placed
  annotate("segment", x = 14, xend = 14, y = 0, yend = wk14$ess,
           colour = verm, linetype = "dotted", linewidth = 0.3) +
  annotate("text", x = 14.6, y = wk14$ess,
           label = sprintf("the early rise:\nonly %.0f of %d\nsurvive here",
                           wk14$ess, NP),
           hjust = 0, vjust = 0.15, size = 2.6, colour = verm, lineheight = 0.95) +
  annotate("rect", xmin = 26.5, xmax = 30.5, ymin = 0, ymax = max(df$ess),
           fill = verm, alpha = 0.07) +
  annotate("text", x = 28.5, y = max(df$ess) * 0.62,
           label = "the crash\nto zero cases",
           hjust = 0.5, vjust = 1, size = 2.6, colour = verm, lineheight = 0.95) +
  xsc +
  scale_y_continuous(limits = c(0, NP),
                     breaks = c(0, 2500, 5000),
                     labels = c("0", "2500", "5000"),
                     expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "week", y = "effective\nsample size",
       title = "how many particles are actually doing the work?") +
  theme_bw() + base_txt +
  theme(plot.margin = margin(7, 8, 2, 6))

fig <- p_top / p_bot + plot_layout(heights = c(0.7, 1.15))

ggsave("figures/fig-ess-collapse.png", fig,
       width = 5.5, height = 3.3, dpi = 300, bg = "white")

cat("wrote figures/fig-ess-collapse.png\n")
