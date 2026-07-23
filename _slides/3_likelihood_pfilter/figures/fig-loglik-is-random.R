## fig-loglik-is-random.R -----------------------------------------------------
## SISMID / EpiSim -- Lecture A_likelihood_pfilter ("Motivation for iterated filtering").
##
## THE FIGURE THAT REPLACES JENSEN'S INEQUALITY.
## Replaces the "unbiased estimate of the likelihood implies a biased estimate
## of the log-likelihood" bullet (3_likelihood/main.qmd:542-559) and the Monte
## Carlo variability prose (4_pfilter_in_pomp/main.qmd:125-151).
##
## Story, told without algebra: 1000 independent pfilter runs at each of
## Np = 1e3, 1e4, 1e5, at the SAME parameter value, on the SAME data. As Np
## grows the cloud of answers gets TIGHTER *and* MARCHES UP toward the truth.
## Too few particles is not just noise -- it is a systematic under-estimate.
##
## Run from the project root (where 3_likelihood_pfilter.Rproj lives):
##   Rscript figures/fig-loglik-is-random.R
## Reads : results/loglikest-pfilter.rds
## Writes: figures/fig-loglik-is-random.png
##
## -----------------------------------------------------------------------------
## CAVEAT ON THE GREEN REFERENCE LINE -- read before presenting.
##
## The green line is logmeanexp(Np=1e5) = -128.18. Treat it as a good but
## slightly LOW estimate of the truth. Do NOT call it "the true value", and do
## NOT claim its confidence interval contains the truth. Why:
##
##  * logmeanexp averages on the LIKELIHOOD scale, so a handful of the largest
##    draws carry nearly all the weight. Effective sample size of that average,
##    out of 1000 runs:  Np=1e3 -> 5.9,  Np=1e4 -> 1.6,  Np=1e5 -> 8.6.
##    The green line is effectively an average over ~9 runs, not 1000.
##  * That is why lme(1e4) = -127.86 sits ABOVE lme(1e5) = -128.18, which looks
##    like the trend reversing. It is not a real reversal: at Np=1e4 the single
##    largest draw carries 79% of the weight, so that number is essentially one
##    lucky observation. This anomaly is deliberately NOT plotted.
##  * The se's are not trustworthy here: the lme(1e3) interval [-131.06,-129.17]
##    is DISJOINT from the lme(1e5) interval [-128.98,-127.38], yet both estimate
##    the same true log-likelihood. A bootstrap resamples the 1000 draws we have
##    and structurally cannot discover the unsampled upper tail, so its interval
##    is not a bracket on the truth either.
##  * Subsampling confirms non-convergence -- logmeanexp(1e5) still climbs
##    monotonically with the number of reps:
##      50 -> -128.58, 100 -> -128.46, 250 -> -128.32, 500 -> -128.29,
##      1000 -> -128.18.  It has not levelled off.
##
## HONEST STATEMENT: the reference is biased LOW by roughly 0.3-0.5 log units;
## the truth is plausibly at or above the upper end of any interval you compute.
## This is immaterial to the 5.4-log-unit teaching point and is conservative in
## direction -- a higher truth means a BIGGER gap and MORE runs falling below the
## line -- which is exactly why the figure survives the caveat.
##
## Corroboration that -128.18 is in the right place: results/pfilter-grid1.rds,
## nearest cell to (Beta=15, Gamma=0.5), spans -132.0 to -124.6. The reference
## sits comfortably inside. The Np=1e3 median of -135.95 also matches the
## documented pfilter(measSIR, Np=1000) ~= -136, confirming this cache is at the
## model's default parameters.
## -----------------------------------------------------------------------------

set.seed(20250715)   # nothing here is random, but the deck must be byte-reproducible

suppressPackageStartupMessages({
  library(tidyverse)
  library(pomp)
})

theme_set(theme_bw())

## Okabe-Ito, fixed course meanings -------------------------------------------
col_pf    <- "#0072B2"  # blue       : the particle filter / the method that works
col_naive <- "#D55E00"  # vermillion : the naive approach / what doesn't work
col_truth <- "#009E73"  # green      : truth / reference
col_grey  <- "grey60"   # de-emphasised context

## ---- data -------------------------------------------------------------------
ll <- readRDS("results/loglikest-pfilter.rds")
stopifnot(all(c("rep", "NP", "loglik") %in% names(ll)),
          nrow(ll) == 3000, all(is.finite(ll$loglik)))
ll <- ll |> mutate(NP = as.numeric(NP))

## ---- the numbers that go on the figure --------------------------------------
smry <- ll |>
  group_by(NP) |>
  summarise(
    n      = n(),
    median = median(loglik),
    mean   = mean(loglik),
    sd     = sd(loglik),
    lo     = min(loglik),
    hi     = max(loglik),
    lme    = logmeanexp(loglik),   # unbiased in the LIKELIHOOD -> best estimate of truth
    lme_se = logmeanexp(loglik, se = TRUE)[2],
    .groups = "drop"
  )

## best available estimate of the true log-likelihood: the Np = 1e5 logmeanexp
ref <- smry$lme[smry$NP == 1e5]

med <- setNames(smry$median, format(smry$NP, scientific = TRUE))
gap_median <- smry$median[smry$NP == 1e3] - smry$median[smry$NP == 1e5]
bias_1e3   <- smry$median[smry$NP == 1e3] - ref

cat("\n---- log-likelihood estimates by particle count (1000 reps each) ----\n")
print(as.data.frame(smry), digits = 6)
cat("\nreference (Np=1e5 logmeanexp) = ", round(ref, 3),
    "  (se ", round(smry$lme_se[smry$NP == 1e5], 3), ")\n", sep = "")
cat("median gap, Np=1e3 vs Np=1e5  = ", round(gap_median, 2), " log units\n", sep = "")
cat("bias of Np=1e3 median vs ref  = ", round(bias_1e3, 2), " log units\n", sep = "")
cat("spread (sd) shrinks ", round(smry$sd[smry$NP == 1e3], 2), " -> ",
    round(smry$sd[smry$NP == 1e5], 2), " log units\n", sep = "")
for (np in c(1e3, 1e4, 1e5))
  cat("Np=", format(np, scientific = TRUE), " fraction of runs below reference: ",
      round(mean(ll$loglik[ll$NP == np] < ref), 3), "\n", sep = "")
cat("\n")

## ---- layout -----------------------------------------------------------------
y_lo <- -146.5
y_hi <- -117.0

x_gap <- 10^3.5      # dead band between the 1e3 and 1e4 violins
y_cap <- -117.6      # internal caption ("what am I looking at")
y_sd  <- -145.0      # spread labels, in the free band below every violin

fig <- ggplot(ll, aes(x = NP, y = loglik)) +

  ## reference: the truth ----------------------------------------------------
  geom_hline(yintercept = ref, colour = col_truth, linewidth = 0.8,
             linetype = "22") +

  ## guides anchoring the two ends of the gap arrow to the medians they -------
  ## actually measure (1e3 and 1e5) -- without these the arrow head floats in
  ## blank space and reads as "the gap between 1e3 and 1e4", which is 3.4.
  annotate("segment", x = 1e3, xend = x_gap,
           y = med["1e+03"], yend = med["1e+03"],
           colour = "grey70", linewidth = 0.3) +
  annotate("segment", x = x_gap, xend = 1e5,
           y = med["1e+05"], yend = med["1e+05"],
           colour = "grey70", linewidth = 0.3) +

  ## the clouds of answers ---------------------------------------------------
  geom_violin(aes(group = NP, fill = factor(NP), colour = factor(NP)),
              width = 0.66, alpha = 0.45, linewidth = 0.35,
              trim = TRUE, scale = "width") +

  ## the march of the medians ------------------------------------------------
  geom_line(data = smry, aes(x = NP, y = median),
            colour = "grey25", linewidth = 0.45) +
  geom_point(data = smry, aes(x = NP, y = median, colour = factor(NP)),
             size = 2.2, show.legend = FALSE) +

  ## the gap: median at 1e3 -> median at 1e5 ---------------------------------
  annotate("segment", x = x_gap, xend = x_gap,
           y = med["1e+03"], yend = med["1e+05"],
           colour = "grey15", linewidth = 0.45,
           arrow = arrow(ends = "both", type = "closed",
                         length = unit(0.05, "in"))) +
  annotate("text", x = x_gap, y = med["1e+03"] - 0.9,
           label = paste0(sprintf("%.1f", abs(gap_median)), "\nlog units\nlow"),
           hjust = 0.5, vjust = 1, size = 2.8, colour = "grey15",
           lineheight = 0.92, fontface = "bold") +

  ## truth label -------------------------------------------------------------
  annotate("text", x = 10^2.56, y = ref + 0.9,
           label = "best estimate of the truth",
           hjust = 0, vjust = 0, size = 2.95, colour = col_truth,
           fontface = "bold") +

  ## setup caption: these are repeated runs at ONE parameter value ------------
  ## "one fixed parameter value" is load-bearing: it stops a student comparing
  ## the -128.18 here against the -107.46 MLE on the ridge slide and concluding
  ## the deck contradicts itself. These runs are NOT at the MLE.
  annotate("text", x = 10^4, y = y_cap,
           ## "particle-filter" is dropped, not lost: the x-axis already reads
           ## "number of particles". Spelling it out here pushes the line to
           ## 1439px inside a 1478px panel, i.e. hard against the border.
           label = "1,000 independent runs at each setting - same model, same data, one fixed parameter value",
           hjust = 0.5, vjust = 1, size = 2.6, colour = "grey35") +

  ## spread labels -----------------------------------------------------------
  geom_text(data = smry, aes(x = NP, y = y_sd,
                             label = paste0("spread ", sprintf("%.1f", sd))),
            inherit.aes = FALSE, size = 2.8, colour = "grey25", vjust = 0.5) +

  ## scales ------------------------------------------------------------------
  scale_x_log10(
    breaks = c(1e3, 1e4, 1e5),
    labels = c("1,000", "10,000", "100,000"),
    expand = expansion(mult = c(0.07, 0.07))
  ) +
  scale_y_continuous(breaks = seq(-145, -120, by = 5),
                     expand = expansion(mult = c(0, 0))) +
  coord_cartesian(ylim = c(y_lo, y_hi)) +
  scale_fill_manual(values   = c("1000" = col_naive, "10000" = col_grey,
                                 "1e+05" = col_pf)) +
  scale_colour_manual(values = c("1000" = col_naive, "10000" = col_grey,
                                 "1e+05" = col_pf)) +
  labs(x = "number of particles", y = "log-likelihood estimate") +
  theme(
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text          = element_text(size = 8),
    axis.title         = element_text(size = 9.5),
    plot.margin        = margin(4, 5, 3, 3)
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/fig-loglik-is-random.png", fig,
       width = 5.5, height = 3.0, dpi = 300, bg = "white")

cat("wrote figures/fig-loglik-is-random.png\n")
