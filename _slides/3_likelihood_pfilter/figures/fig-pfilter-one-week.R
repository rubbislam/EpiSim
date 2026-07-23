## fig-pfilter-one-week.R -------------------------------------------------
## ONE observation interval of the Consett 1948 measles particle filter,
## shown on real data: PREDICT -> WEIGHT -> RESAMPLE.
##
## Replaces the prediction / filtering / two-step-recursion display equations
## (old 3_likelihood/main.qmd:434-468).
##
## Run from the project root (where 3_likelihood_pfilter.Rproj lives):
##   Rscript figures/fig-pfilter-one-week.R
## Reads : data/Measles_Consett_1948.csv (via scripts/model_measSIR.R)
##         results/pfilter-grid1.rds     (for the MLE parameter vector)
## Writes: figures/fig-pfilter-one-week.png

source("scripts/model_measSIR.R")
library(ggplot2)
library(dplyr)

theme_set(theme_bw())

## Okabe-Ito, fixed course meanings
col_pf   <- "#0072B2"   # the particle filter
col_bad  <- "#D55E00"   # what doesn't work: particles the data kill
col_data <- "#E69F00"   # THE DATA

## ---- settings ----------------------------------------------------------
## Week 16 -> 17 is the interval where the weighting visibly bites: the
## model's cloud is still growing (median ~38 expected reports) while the
## data say 18, so the filter has to drag the cloud down onto the data.
## Weeks 19-21 are useless for teaching (ESS ~95% of J: nothing happens).
## The seed is fixed so the step is legible. This draw is sharper than the
## median for the week (12 killed / 4 copies / ESS 51%, against a week-17
## median of 10 / 3 / 61%), but nothing shown is unusual for the week.
n_part <- 30     # few particles so individuals stay visible
n_week <- 17     # the observation interval shown: week 16 -> week 17
seed   <- 238

## ---- parameters: the MLE from the cached Beta-Gamma grid ---------------
## DELIBERATE: the grid MLE, not coef(measSIR). Do not "fix" this to the
## defaults. At the defaults the week-17 datum falls outside the particle
## cloud entirely, so every particle is bad and the figure teaches "the
## model is wrong" -- true, and Lecture A's LATER punchline, but not the
## point of THIS slide, which is the mechanics of one filtering step.
## At the MLE the datum sits inside the cloud, so kill-and-clone is visible.
readRDS("results/pfilter-grid1.rds") |>
  group_by(Beta, Gamma, Rho, k, Eta, N) |>
  summarise(loglik = logmeanexp(loglik), .groups = "drop") |>
  slice_max(loglik, n = 1) -> best
theta <- c(Beta = best$Beta, Gamma = best$Gamma, Rho = best$Rho,
           k = best$k, Eta = best$Eta, N = best$N)

## ---- 0. filtering particles at week n-1 --------------------------------
set.seed(seed)
pf <- pfilter(measSIR, params = theta, Np = n_part, save.states = "filter")
x_prev <- saved_states(pf)[[n_week - 1]]           # 4 x J filtering particles

## emeasure wants a 3-d array: nvar x nrep x ntime
x_prev_a <- array(x_prev, dim = c(nrow(x_prev), ncol(x_prev), 1L),
                  dimnames = list(name = rownames(x_prev), NULL, NULL))
e_prev <- as.numeric(emeasure(measSIR, x = x_prev_a,
                              times = n_week - 1, params = theta)[1, , 1])

## ---- 1. PREDICT: simulate each particle forward one week ---------------
set.seed(seed * 1000 + n_week)
x_pred <- rprocess(measSIR, x0 = x_prev, t0 = n_week - 1,
                   times = n_week, params = theta)
e_pred <- as.numeric(emeasure(measSIR, x = x_pred,
                              times = n_week, params = theta)[1, , 1])

## ---- 2. WEIGHT: the datum arrives --------------------------------------
y_n   <- obs(measSIR)[, n_week, drop = FALSE]
ystar <- as.numeric(y_n)
w  <- as.numeric(dmeasure(measSIR, y = y_n, x = x_pred,
                          times = n_week, params = theta, log = FALSE))
wn  <- w / sum(w)
ess <- 1 / sum(wn^2)

## ---- 3. RESAMPLE: systematic resampling (what pomp does internally) ----
systematic <- function (wn, n, u) findInterval((u + seq_len(n) - 1) / n,
                                               cumsum(wn)) + 1L
set.seed(seed * 1000 + n_week + 7)
idx  <- systematic(wn, n_part, runif(1))
mult <- tabulate(idx, nbins = n_part)          # copies made of each parent

n_kept   <- sum(mult > 0)
n_killed <- sum(mult == 0)
max_mult <- max(mult)

## ---- report every number that reaches the figure -----------------------
cat("\n=== fig-pfilter-one-week: computed numbers ===\n")
cat(sprintf("params (grid MLE): Beta=%.5f Gamma=%.6f Rho=%.1f k=%d Eta=%.2f N=%d\n",
            theta["Beta"], theta["Gamma"], theta["Rho"], theta["k"],
            theta["Eta"], theta["N"]))
cat(sprintf("  cached grid loglik at this cell = %.2f\n", best$loglik))
cat(sprintf("interval: week %d -> %d ; y* = %d reports ; J = %d ; seed = %d\n",
            n_week - 1, n_week, ystar, n_part, seed))
cat(sprintf("week %d filter cloud (expected reports): %.1f to %.1f (median %.1f)\n",
            n_week - 1, min(e_prev), max(e_prev), median(e_prev)))
cat(sprintf("week %d PREDICT cloud: %.1f to %.1f (median %.1f)\n",
            n_week, min(e_pred), max(e_pred), median(e_pred)))
cat(sprintf("weights: min=%.3g max=%.3g ; max normalised=%.1f%% (uniform %.1f%%)\n",
            min(w), max(w), 100 * max(wn), 100 / n_part))
cat(sprintf("ESS = %.1f of %d\n", ess, n_part))
cat(sprintf("resample: %d kept, %d killed, most copies of one particle = %d\n",
            n_kept, n_killed, max_mult))
cat(sprintf("RESAMPLED cloud: %.1f to %.1f (median %.1f)\n",
            min(e_pred[idx]), max(e_pred[idx]), median(e_pred[idx])))
cat(sprintf("mean |E - y*|: predict %.1f -> resampled %.1f  (%.0f%% closer)\n",
            mean(abs(e_pred - ystar)), mean(abs(e_pred[idx] - ystar)),
            100 * (1 - mean(abs(e_pred[idx] - ystar)) / mean(abs(e_pred - ystar)))))
cat(sprintf("sd of cloud: predict %.1f -> resampled %.1f\n",
            sd(e_pred), sd(e_pred[idx])))
cat("=============================================\n\n")

## ---- assemble plotting frames ------------------------------------------
x0 <- 0; x1 <- 1; x2 <- 2; x3 <- 3          # the four columns

seg1 <- data.frame(x = x0, xend = x1, y = e_prev, yend = e_pred)
d0   <- data.frame(x = x0, y = e_prev)
d1   <- data.frame(x = x1, y = e_pred)
d2   <- data.frame(x = x2, y = e_pred, w = wn)

## survivors fan out horizontally: one dot per copy, so clones are countable.
## The killed get their own narrow sub-column so a cross can never land on
## top of a survivor that happens to share its value.
x_surv <- x3 - 0.13
x_dead <- x3 + 0.24
surv <- bind_rows(lapply(which(mult > 0), function (j) {
  m <- mult[j]
  data.frame(x = x_surv + (seq_len(m) - (m + 1) / 2) * 0.105, y = e_pred[j])
}))
dead <- data.frame(x = x_dead, y = e_pred[mult == 0])

## the particle cloned the most: label it
j_star <- which.max(mult)
y_star_part <- e_pred[j_star]

ytop <- max(e_pred) + 7
ybot <- min(e_pred) - 3

p <- ggplot() +
  ## the datum: amber reference, starting at PREDICT. It must NOT reach back
  ## under the week-16 column: y* = 18 is week 17's datum, week 16's is 22.
  annotate("segment", x = x1 - 0.15, xend = x_surv + 0.19, y = ystar, yend = ystar,
           colour = col_data, linewidth = 0.7, linetype = "22") +
  ## stage 1: each particle simulates forward -> a spreading cloud
  geom_segment(data = seg1, aes(x = x, xend = xend, y = y, yend = yend),
               colour = col_pf, alpha = 0.3, linewidth = 0.3) +
  geom_point(data = d0, aes(x, y), colour = col_pf, size = 1.5, alpha = 0.85) +
  geom_point(data = d1, aes(x, y), colour = col_pf, size = 1.5, alpha = 0.85) +
  ## stage 2: grey ring = the particle; blue area = its weight
  geom_point(data = d2, aes(x, y), shape = 21, colour = "grey65",
             fill = NA, size = 1.6, stroke = 0.3) +
  geom_point(data = d2, aes(x, y, size = w), colour = col_pf, alpha = 0.9) +
  scale_size_area(max_size = 3.4, guide = "none") +
  annotate("point", x = x2, y = ystar, colour = col_data,
           size = 3.4, shape = 18) +
  ## stage 3: the new cloud, and the particles the data killed
  geom_point(data = surv, aes(x, y), colour = col_pf, size = 1.5, alpha = 0.9) +
  geom_point(data = dead, aes(x, y), colour = col_bad, shape = 4,
             size = 1.8, stroke = 0.65) +
  ## annotations
  annotate("text", x = x1 + 0.08, y = ystar - 3.9,
           label = paste0("the data: ", ystar, " reports"),
           colour = col_data, size = 2.8, fontface = "bold", hjust = 0) +
  annotate("text", x = x_dead + 0.16, y = max(e_pred) + 1,
           label = paste0(n_killed, " killed:\ntoo far from\nthe data"),
           colour = col_bad, size = 2.7, fontface = "bold",
           hjust = 0, vjust = 1, lineheight = 0.95) +
  ## the best particle sits on the data, so label its clone row directly
  annotate("text", x = x_surv + 0.25, y = y_star_part,
           label = paste0("cloned ", max_mult, "x"),
           colour = col_pf, size = 2.7, fontface = "bold",
           hjust = 0, vjust = 0.5) +
  scale_x_continuous(
    breaks = c(x0, x1, x2, x3),
    labels = c("particles at\nweek 16",
               "1. PREDICT\nsimulate to week 17",
               "2. WEIGHT\ndot area = fit to data",
               "3. RESAMPLE\nkill and clone"),
    limits = c(x0 - 0.3, x3 + 1.05)
  ) +
  scale_y_continuous(limits = c(ybot, ytop)) +
  labs(x = NULL, y = "expected weekly reports") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 7.6, colour = "grey20"),
    axis.text.y = element_text(size = 7.5),
    axis.title.y = element_text(size = 9),
    plot.margin = margin(4, 4, 2, 4)
  )

ggsave("figures/fig-pfilter-one-week.png", p,
       width = 5.5, height = 3.0, dpi = 300, bg = "white")

cat("wrote figures/fig-pfilter-one-week.png\n")
