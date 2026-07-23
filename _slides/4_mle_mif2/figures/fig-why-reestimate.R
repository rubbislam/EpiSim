## fig-why-reestimate.R ---------------------------------------------------
## B_mle_mif2, Lecture B. Answers: "Why do we need to estimate the likelihood in a
## different manner?"  i.e. why the number mif2 prints is NOT your likelihood.
##
## During mif2 the parameters are perturbed at every observation. The value
## mif2 reports is therefore the likelihood of a model that has extra noise
## in it -- not the likelihood of the actual model at the final parameters.
## The honest number comes from re-running a plain pfilter (several times,
## then logmeanexp) at the final parameter vector.
##
## Precisely (verified against the cache):
##   vermillion = logLik(mif) == traces(mif)[50,"loglik"], the iteration-50
##                filtering pass over a PERTURBED parameter swarm;
##   blue       = pfilter at coef(mif) == traces(mif)[51,], the final vector.
## Trace row 51 has loglik = NA in all 20 reps: mif2 never filtered at the
## parameters it hands you. The iter50 -> iter51 move is 1.2%-6.2% relative
## (per-rep max over Beta/Rho/Eta), so these are NOT the same parameters.
##
## Run from the project root (where 4_mle_mif2.Rproj lives):
##   Rscript figures/fig-why-reestimate.R
## Reads : results/local_search.rds, results/lik_local.rds
## Writes: figures/fig-why-reestimate.png
## Nothing here is random; the figure is byte-reproducible. (No seed needed.)

suppressPackageStartupMessages({
  library(pomp)
  library(tidyverse)
})

theme_set(theme_bw())

## Okabe-Ito, course-fixed meanings ---------------------------------------
col_pf   <- "#0072B2"  # blue       = the particle filter / the method that works
col_mif  <- "#D55E00"  # vermillion = the misleading number / what doesn't work
col_grey <- "grey60"

## ------------------------------------------------------------------------
## 1. The two numbers, per replicate
## ------------------------------------------------------------------------
local_search <- readRDS("results/local_search.rds")   # mif2List, 20 reps, Np=2000, Nmif=50
lik_local    <- readRDS("results/lik_local.rds")      # honest pfilter re-estimate at final params

stopifnot(length(local_search) == nrow(lik_local))
## confirm row i of lik_local really is replicate i of local_search
stopifnot(all.equal(unname(sapply(local_search, coef, "Beta")), lik_local$Beta))

dat <- tibble(
  rep = seq_along(local_search),
  mif = sapply(local_search, logLik),   # what mif2 reported at its last iteration
  pf  = lik_local$loglik,               # logmeanexp of plain pfilter reps at final params
  se  = lik_local$loglik.se
) |>
  mutate(gap = pf - mif)

## ------------------------------------------------------------------------
## 2. Real numbers (printed to console; used verbatim in the annotations)
## ------------------------------------------------------------------------
n_up      <- sum(dat$gap > 0)
n_rep     <- nrow(dat)
gap_med   <- median(dat$gap)
gap_max   <- max(dat$gap)
gap_min   <- min(dat$gap)
se_mean   <- mean(dat$se)
best_mif  <- which.max(dat$mif)
best_pf   <- which.max(dat$pf)

cat("n replicates             :", n_rep, "\n")
cat("pfilter > mif2 in        :", n_up, "/", n_rep, "replicates\n")
cat("gap (pf - mif) median    :", round(gap_med, 3), "\n")
cat("gap range                :", round(gap_min, 3), "to", round(gap_max, 3), "\n")
cat("mean Monte Carlo se of pf:", round(se_mean, 3), "\n")
cat("best replicate by mif2   :", best_mif, "  by pfilter:", best_pf, "\n")

## ------------------------------------------------------------------------
## 3. Order rows by the honest number, best at the top
## ------------------------------------------------------------------------
dat <- dat |>
  arrange(pf) |>
  mutate(row = row_number())          # row 1 = worst, row 20 = best (top of plot)

## Anchors are FRACTIONS of the data span, never absolute log-unit offsets.
## An absolute offset silently pushes text outside the panel the moment the
## cache changes -- ggplot then DROPS the annotation with only a "Removed 1
## row" warning, and the punchline vanishes from the slide unnoticed.
x_lo  <- min(dat$mif) - 0.35
x_hi  <- max(dat$pf)  + 0.35
x_sp  <- x_hi - x_lo
## The key goes ABOVE the rows rather than into a "corner": which corner is
## empty depends on the data, but the strip above the top row never has any.
y_key <- n_rep + 2.4

## ------------------------------------------------------------------------
## 4. The figure
## ------------------------------------------------------------------------
p <- ggplot(dat, aes(y = row)) +
  ## the move from the printed number to the honest number
  geom_segment(
    aes(x = mif, xend = pf, yend = row),
    colour = col_grey, linewidth = 0.45,
    arrow = arrow(length = unit(0.05, "in"), type = "closed"),
    lineend = "butt", linejoin = "mitre"
  ) +
  geom_point(aes(x = mif), colour = col_mif, size = 1.7) +
  geom_point(aes(x = pf),  colour = col_pf,  size = 1.7) +
  ## key, drawn in the clear strip ABOVE the rows -- no legend needed
  annotate("point", x = x_lo + 0.03 * x_sp, y = y_key, colour = col_mif, size = 1.7) +
  annotate("text",  x = x_lo + 0.07 * x_sp, y = y_key, hjust = 0, size = 2.75,
           colour = col_mif, fontface = "bold",
           label = "the number mif2 printed") +
  annotate("point", x = x_lo + 0.52 * x_sp, y = y_key, colour = col_pf, size = 1.7) +
  annotate("text",  x = x_lo + 0.56 * x_sp, y = y_key, hjust = 0, size = 2.75,
           colour = col_pf, fontface = "bold",
           label = "a plain pfilter at the final parameters") +
  ## punchline, bottom-right, anchored to the right edge so it cannot fall off
  annotate(
    "text", x = x_hi - 0.02 * x_sp, y = 8, hjust = 1, vjust = 1,
    size = 2.75, colour = "grey20", lineheight = 1.15,
    label = paste0(n_up, " of ", n_rep, "\nruns move the same way.\n",
                   "mif2 reads low by\n", sprintf("%.2f", gap_med),
                   " log units\n(median; up to ", sprintf("%.1f", gap_max), ").")
  ) +
  scale_x_continuous(
    limits = c(x_lo, x_hi),
    breaks = scales::breaks_width(0.5)
  ) +
  ## clip="off" is NOT used: everything must live inside the panel. The y range
  ## is opened at the top to make the key's strip, not to overflow.
  coord_cartesian(ylim = c(0.4, y_key + 1.0)) +
  scale_y_continuous(breaks = NULL) +
  labs(
    x = "log likelihood     (each row = one of 20 mif2 runs)",
    y = NULL
  ) +
  theme(
    axis.text.x  = element_text(size = 7.5),
    axis.text.y  = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin  = margin(4, 6, 3, 3)
  )

ggsave(
  "figures/fig-why-reestimate.png", p,
  width = 5.5, height = 3.0, dpi = 300, bg = "white"
)

cat("wrote figures/fig-why-reestimate.png\n")
