## fig-direct-sim-fails.R -----------------------------------------------------
## A_likelihood_pfilter, Lecture A.  "ATTEMPT 1: direct simulation" -- the figure that
## kills guess-and-check.  Replaces the direct-simulation Monte Carlo algebra
## (L(theta) = E[prod_n f(y*_n|X_n)]) in 3_likelihood/main.qmd:345-388.
## Complements ../graphics/directsim_graph.pdf, which shows the IDEA; this shows
## the REALITY on the Consett 1948 data.
##
## Run from the project root (where 3_likelihood_pfilter.Rproj lives):
##   Rscript figures/fig-direct-sim-fails.R
## Reads : data/Measles_Consett_1948.csv (via scripts/model_measSIR.R)
## Writes: figures/fig-direct-sim-fails.png
##
## -- Parameter choice ---------------------------------------------------------
## We simulate at the best parameter we have for this model (Beta=22.31,
## Gamma=1.049; the logmeanexp MLE of results/pfilter-grid1.rds, particle-filter
## loglik ~ -107.5), NOT at the shipped default Beta=15, Gamma=0.5.  Deliberate,
## and it is the fair test: if unconditional simulation cannot estimate the
## likelihood even where the model fits best, it cannot do it anywhere.  (At the
## default, direct simulation lands near the truth by luck -- ESS 1.5 -- which
## would teach students exactly the wrong lesson.)
##
## -- What the figure actually argues ------------------------------------------
## The naive intuition to kill is "plenty of my simulations look like the data,
## so I am sampling the right region".  Some genuinely do: 150 of the 500 produce
## a real epidemic (peak >= 5) and 53 of them peak between 50 and 110 reports
## somewhere in weeks 15-24, i.e. right on top of the observed peak.  But the
## likelihood needs ONE trajectory to clear ALL 42 weeks at once, and only 27 of
## the 500 manage it -- one of those carries 92% of the weight.  The funnel
## 150 -> 27 -> 1 is the curse of dimensionality in plain sight, and it is
## exactly what resampling (the particle filter) fixes.  Do not retitle the left
## panel "nothing goes near the data": that is false here, and the render shows
## it is false (at week 20, ten simulations sit at or above the observed 47).
##
## -- Why there is no "middle 95%" band statistic here --------------------------
## An earlier version annotated "the data sits inside the middle 95% of the
## simulations at 37 of 42 weeks".  That number is real but VACUOUS: 350 of the
## 500 never take off, so the 2.5% quantile is 0 at every one of the 42 weeks and
## the lower bound never binds (sum(dat < lo) == 0).  "Inside the middle 95%" was
## therefore arithmetically identical to "at or below the 97.5th percentile", and
## it invited the false reading "the simulations bracket the data".  Restricting
## to the 150 that take off does NOT repair it -- the 2.5% quantile is still 0 at
## all 42 weeks there too (verified), so the "39 of 42" version of the statistic
## carries exactly the same defect.  Any band statistic on this ensemble is
## vacuous by construction.  We count take-offs instead, which is a real claim.
## -----------------------------------------------------------------------------

suppressMessages({
  library(tidyverse)
  library(pomp)
  library(patchwork)
})

theme_set(theme_bw())

## Okabe-Ito, course-fixed meanings
col_data <- "#E69F00"  # amber      : THE DATA (Consett reports), always
col_bad  <- "#D55E00"  # vermillion : the naive approach / what doesn't work
col_grey <- "grey60"   # grey       : de-emphasised context / latent

source("scripts/model_measSIR.R")

## ---- 1. the experiment ------------------------------------------------------
NSIM <- 500
SEED <- 1948

theta <- coef(measSIR)
theta["Beta"]  <- 22.31
theta["Gamma"] <- 1.049

set.seed(SEED)
sims <- simulate(measSIR, params = theta, nsim = NSIM, format = "arrays")

## log weight of each simulation: sum_n log f(y*_n | X_n)
ld <- dmeasure(
  measSIR,
  y      = obs(measSIR),
  x      = sims$states,
  times  = time(measSIR),
  params = theta,
  log    = TRUE
)                                   # NSIM x N
logw <- rowSums(ld)

## normalised weights (stable); simulations with logw = -Inf get weight 0
w  <- exp(logw - max(logw))
w  <- w / sum(w)
ws <- sort(w, decreasing = TRUE)

## ---- 2. the real numbers ----------------------------------------------------
ess       <- 1 / sum(w^2)
top1      <- ws[1]
n_ok      <- sum(is.finite(logw))     # sims with nonzero prob at EVERY week
best      <- which.max(logw)
ll_direct <- logmeanexp(logw)

wk   <- time(measSIR)
N_WK <- length(wk)
rep  <- sims$obs["reports", , ]                   # NSIM x N
dat  <- as.numeric(obs(measSIR)[1, ])

## Did the simulation produce an epidemic at all, or fizzle out?  A take-off is
## a run that ever exceeds TAKEOFF reports.  No band/quantile statistic is used:
## with 350 fizzles the 2.5% quantile is 0 at every week, so any "middle 95%"
## claim here is vacuous (see header).  These two counts partition the 500.
TAKEOFF <- 5
peak    <- apply(rep, 1, max)
is_take <- peak >= TAKEOFF
n_take  <- sum(is_take)                  # produced a real epidemic
n_fizz  <- sum(!is_take)                 # never took off

## backs the words "several track the data closely": peak of roughly the right
## size, at roughly the right time (data peaks at 75 in week 18).
pk_week <- wk[apply(rep, 1, which.max)]
n_close <- sum(peak >= 50 & peak <= 110 & pk_week >= 15 & pk_week <= 24)

stopifnot(n_take + n_fizz == NSIM, n_ok <= n_take)

message(sprintf("ESS                 = %.2f of %d", ess, NSIM))
message(sprintf("top-1 weight        = %.1f%%", 100 * top1))
message(sprintf("took off (peak>=%d)  = %d of %d", TAKEOFF, n_take, NSIM))
message(sprintf("fizzled             = %d of %d", n_fizz, NSIM))
message(sprintf("peak-matching sims  = %d", n_close))
message(sprintf("finite-lik sims     = %d of %d", n_ok, NSIM))
message(sprintf("the funnel          = %d -> %d -> 1", n_take, n_ok))
message(sprintf("direct-sim loglik   = %.2f", ll_direct))
message(sprintf("winner sim          = %d", best))

## ---- 3. left panel: the spaghetti -------------------------------------------
spag <- as.data.frame(as.table(rep)) |>
  setNames(c("sim", "tix", "reports")) |>
  mutate(week = wk[as.integer(tix)], sim = as.integer(sim),
         took_off = is_take[sim])

d_dat <- data.frame(week = wk, reports = dat)
d_win <- data.frame(week = wk, reports = rep[best, ])

ymax <- max(rep, dat)

p_left <- ggplot() +
  ## The fizzles are drawn thicker, NOT darker.  All 350 lie on top of each other
  ## along y = 0, so at any alpha they composite to a solid line; the reason the
  ## single commonest outcome was invisible is that it was a HAIRLINE, not that it
  ## was faint.  Width is what makes 70% of the ensemble read as dense.
  geom_line(
    data = filter(spag, !took_off), aes(week, reports, group = sim),
    colour = col_grey, alpha = 0.10, linewidth = 0.45
  ) +
  geom_line(
    data = filter(spag, took_off), aes(week, reports, group = sim),
    colour = col_grey, alpha = 0.15, linewidth = 0.18
  ) +
  ## white halos so the two signal lines lift out of the cloud
  geom_line(data = d_win, aes(week, reports), colour = "white", linewidth = 1.5) +
  geom_line(data = d_win, aes(week, reports), colour = col_bad, linewidth = 0.65) +
  geom_line(data = d_dat, aes(week, reports), colour = "white", linewidth = 1.9) +
  geom_line(data = d_dat, aes(week, reports), colour = col_data, linewidth = 1.05) +
  ## The setup, and it is a true claim, not a band statistic: 150 runs really do
  ## produce an epidemic and 53 of them peak at the right size at the right time.
  ## This is the belief the right panel then demolishes.  Boxed in translucent
  ## white so the spaghetti behind it cannot muddy the text.
  annotate("label", x = 0, y = ymax * 1.03,
           label = sprintf(
             "%d of the %d produce a real epidemic;\nseveral track the data closely",
             n_take, NSIM),
           hjust = 0, vjust = 1, size = 2.6, colour = "grey30",
           lineheight = 1.0, fill = alpha("white", 0.82), border.colour = NA,
           label.padding = unit(0.05, "lines")) +
  ## Name the zero line: it is the single commonest outcome and it is otherwise
  ## just "the x-axis" to a viewer.  Parked in weeks 33-42 above the flat amber
  ## tail (data <= 4 there, 0.03*ymax); only 12 of the 500 ever enter this box,
  ## and the translucent fill covers those.
  annotate("label", x = 42, y = ymax * 0.24,
           label = sprintf("%d of the %d\nnever take off", n_fizz, NSIM),
           hjust = 1, vjust = 1, size = 2.7, colour = "grey30",
           lineheight = 1.0, fill = alpha("white", 0.82), border.colour = NA,
           label.padding = unit(0.05, "lines")) +
  annotate("segment", x = 38.5, xend = 38.5,
           y = ymax * 0.115, yend = ymax * 0.025,
           colour = "grey55", linewidth = 0.3,
           arrow = arrow(length = unit(0.045, "in"), type = "closed")) +
  ## Colour key (no legend: colours are named here and in the talk).  Parked in
  ## the far-right column, which is verifiably empty: no simulation exceeds 49
  ## reports after week 31 (0.36*ymax), the setup note spans only weeks 0-27, and
  ## the "never take off" label sits far below at 0.24*ymax, so nothing shares an
  ## x-range or a baseline with these.  Plain text, no white boxes -- boxes here
  ## punched holes in the vermillion line at its peak.
  annotate("text", x = 41.8, y = ymax * 0.80,
           label = "the data", hjust = 1, vjust = 1,
           size = 2.9, colour = col_data, fontface = "bold") +
  annotate("text", x = 41.8, y = ymax * 0.70,
           label = "the one\nthat counted", hjust = 1, vjust = 1,
           size = 2.9, colour = col_bad, fontface = "bold", lineheight = 1.0) +
  scale_x_continuous(breaks = seq(0, 40, 10), expand = expansion(c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(c(0.01, 0.07))) +
  labs(x = "week", y = "weekly reports",
       title = "500 simulations. Looks promising.") +
  theme(
    plot.title  = element_text(size = 9.5, face = "bold", hjust = 0),
    axis.title  = element_text(size = 9),
    axis.text   = element_text(size = 7.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(2, 4, 2, 2)
  )

## ---- 4. right panel: where the weight actually went -------------------------
cum   <- data.frame(n = 0:NSIM, share = c(0, cumsum(ws)))
ideal <- data.frame(n = c(0, NSIM), share = c(0, 1))

p_right <- ggplot() +
  geom_line(data = ideal, aes(n, share),
            colour = col_grey, linetype = "22", linewidth = 0.45) +
  geom_step(data = cum, aes(n, share),
            colour = col_bad, linewidth = 0.85) +
  annotate("point", x = 1, y = top1, colour = col_bad, size = 1.5) +
  annotate("text", x = 34, y = 0.90,
           label = sprintf("1 simulation carries %.0f%%\nof the total weight",
                           100 * top1),
           hjust = 0, vjust = 1, size = 2.85, colour = col_bad,
           lineheight = 1.0, fontface = "bold") +
  annotate("text", x = 300, y = 0.415,
           label = "if all 500\ncounted equally",
           hjust = 0.5, vjust = 1, size = 2.7, colour = "grey40",
           lineheight = 1.0) +
  annotate("segment", x = 300, xend = 300, y = 0.44, yend = 0.575,
           colour = "grey55", linewidth = 0.3,
           arrow = arrow(length = unit(0.045, "in"), type = "closed")) +
  annotate("text", x = 496, y = 0.175,
           label = sprintf("effective sample size: %.1f of %d", ess, NSIM),
           hjust = 1, vjust = 0.5, size = 2.9, colour = "grey10",
           fontface = "bold") +
  annotate("text", x = 496, y = 0.015,
           label = sprintf("only %d of the %d clear\nall %d weeks at once",
                           n_ok, NSIM, N_WK),
           hjust = 1, vjust = 0, size = 2.7, colour = "grey35",
           lineheight = 1.0) +
  scale_x_continuous(breaks = seq(0, 500, 100), expand = expansion(c(0.015, 0.025))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(c(0.01, 0.05))) +
  labs(x = "number of simulations (heaviest first)",
       y = "share of total weight",
       title = "Essentially one of them counted.") +
  theme(
    plot.title  = element_text(size = 9.5, face = "bold", hjust = 0),
    axis.title  = element_text(size = 9),
    axis.text   = element_text(size = 7.5),
    panel.grid.minor = element_blank(),
    plot.margin = margin(2, 2, 2, 4)
  )

## ---- 5. assemble ------------------------------------------------------------
fig <- p_left | p_right

ggsave(
  "figures/fig-direct-sim-fails.png", fig,
  width = 5.5, height = 3.0, dpi = 300, bg = "white"
)

message("wrote figures/fig-direct-sim-fails.png")
