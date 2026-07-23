## fig-data-vs-hidden.R ---------------------------------------------------
## A_likelihood_pfilter (Lecture A) -- the figure that REPLACES the joint-density
## factorization and the L(theta) = integral ... dx_{0:N} formula.
##
## TOP:    what we HAVE  -- the 42 weekly Consett measles reports (amber).
## BOTTOM: what we DON'T -- the full latent SIR state, 50 draws each, all three
##         compartments (S/I/R, coloured, alpha 0.5) on a sqrt y-scale so R
##         (~36k) does not crush S (~2k) and I (~150). We never see any of it.
##
## The talk: "The likelihood asks -- across every possible hidden epidemic,
##  how plausible are MY 42 numbers? That's the integral. We are not going to
##  write it down. We are going to compute it."
##
## Run from the project root (where 3_likelihood_pfilter.Rproj lives):
##   Rscript figures/fig-data-vs-hidden.R
## Reads : data/ (via scripts/model_measSIR.R), results/pfilter-grid1.rds
## Writes: figures/fig-data-vs-hidden.png
## Runtime: ~5 s.

library(tidyverse)
library(pomp)
library(patchwork)

source("scripts/model_measSIR.R")

## Okabe-Ito, fixed course meanings ---------------------------------------
amber <- "#E69F00"  # THE DATA, always
grey_latent <- "grey45"  # latent / unobserved states

## The latent paths are drawn at the MLE, taken straight from the cached
## 40x40 Beta-Gamma grid (Np=5000, 9 reps/cell -> 14400 rows) rather than
## hard-coded, so the figure tracks the cache. Rho, k, Eta and N are held
## fixed across that grid, so they carry through unchanged. This yields
## Beta=22.31, Gamma=1.049, loglik=-107.46.
grid <- readRDS("results/pfilter-grid1.rds")
stopifnot(sapply(grid[, c("Rho", "k", "Eta", "N")],
                 function(z) length(unique(z))) == 1)

mle_cell <- grid |>
  group_by(Beta, Gamma) |>
  summarise(loglik = logmeanexp(loglik), .groups = "drop") |>
  slice_max(loglik, n = 1)

mle <- c(Beta = mle_cell$Beta, Gamma = mle_cell$Gamma,
         Rho = grid$Rho[1], k = grid$k[1],
         Eta = grid$Eta[1], N = grid$N[1])

cat("MLE from grid: Beta =", signif(mle[["Beta"]], 6),
    " Gamma =", signif(mle[["Gamma"]], 6),
    " loglik =", round(mle_cell$loglik, 2), "\n")

n_traj <- 50
np <- 2000

## The data: 42 weekly reports ---------------------------------------------
dat <- as.data.frame(measSIR)   # week, reports

## The hidden epidemics ----------------------------------------------------
## Each pfilter run with filter.traj=TRUE returns ONE draw from the
## smoothing distribution -- i.e. one latent path that is consistent with the
## observed 42 numbers. 50 independent runs => 50 independent hidden
## epidemics, any of which could have produced this data. Unconditional
## simulate() is NOT used here on purpose: from I=1 most epidemics fizzle,
## which is the *next* figure's story (attempt 1 fails), not this one's.
set.seed(1350254336)

## filter_traj() leaves the time dimension unnamed, so build it from the
## pomp object: t0 (week 0) followed by the 42 observation times.
wk <- c(timezero(measSIR), time(measSIR))

traj <- map_dfr(seq_len(n_traj), function(i) {
  pf <- pfilter(measSIR, params = mle, Np = np, filter.traj = TRUE)
  ft <- filter_traj(pf)                       # [var, rep, time]
  stopifnot(dim(ft)[3] == length(wk))
  tibble(
    .id  = i,
    week = wk,
    S    = as.numeric(ft["S", 1, ]),
    I    = as.numeric(ft["I", 1, ]),
    R    = as.numeric(ft["R", 1, ])
  )
})

## ---- REAL NUMBERS (printed, and used in the report) ---------------------
peak_I    <- traj |> group_by(.id) |> summarise(pk = max(I), wk = week[which.max(I)])
n_obs     <- nrow(dat)
data_peak <- max(dat$reports)
data_peak_wk <- dat$week[which.max(dat$reports)]

cat("\n---- verified numbers ----\n")
cat("observations                :", n_obs, "\n")
cat("data peak                   :", data_peak, "reports at week", data_peak_wk, "\n")
cat("data total reports          :", sum(dat$reports), "\n")
cat("latent trajectories         :", n_traj, "(Np =", np, "each)\n")
cat("peak I, median              :", median(peak_I$pk), "\n")
cat("peak I, range               :", paste(range(peak_I$pk), collapse = " - "), "\n")
cat("peak I week, median         :", median(peak_I$wk), "\n")
cat("peak I week, range          :", paste(range(peak_I$wk), collapse = " - "), "\n")
cat("max I across all traj       :", max(traj$I), "\n")
cat("--------------------------\n\n")

## ---- panels -------------------------------------------------------------
base_txt <- theme(
  axis.text  = element_text(size = 8),
  axis.title = element_text(size = 9.5),
  plot.title = element_text(size = 10.5, face = "bold",
                            margin = margin(b = 3)),
  plot.margin = margin(3, 6, 2, 6)
)

xsc <- scale_x_continuous(limits = c(0, 42.5), breaks = seq(0, 40, 10),
                          expand = expansion(mult = c(0.01, 0.02)))

## The headline for each panel goes in the title, not inside the panel: at
## this size any in-panel banner collides with the data (top) or with the
## peak of the latent band (bottom).

## TOP: what we have. Solid, certain, real.
p_top <- ggplot(dat, aes(week, reports)) +
  geom_line(colour = amber, linewidth = 0.6) +
  geom_point(colour = amber, size = 1.5) +
  annotate("text", x = 0.5, y = 76, hjust = 0, vjust = 1,
           label = "this is all we ever get",
           colour = "grey35", size = 2.9) +
  xsc +
  scale_y_continuous(limits = c(0, 82), breaks = seq(0, 75, 25),
                     expand = expansion(mult = c(0.02, 0.03))) +
  labs(y = "weekly reports", title = "observed: 42 numbers") +
  theme_bw() + base_txt +
  ## Title is grey25 (10.4:1 on white), NOT amber: amber text is 2.25:1 and is
  ## unreadable from the back of the room. The amber data marks immediately
  ## below already carry the amber=data association; thick graphical marks and
  ## text are judged differently. Matches the bottom panel title.
  theme(plot.title = element_text(size = 10.5, face = "bold", colour = "grey25",
                                  margin = margin(b = 3)),
        axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

## BOTTOM: what we don't have -- the full latent SIR state. TWO y-axes:
## S and I share the LEFT axis (0..~2.3k); R lives near 36k and is read off
## the RIGHT axis. R is linearly rescaled into left-axis coordinates for
## drawing, and sec_axis restores its true scale on the right.
comp_cols <- c(S = "#0072B2", I = "#D55E00", R = "#009E73")
green <- comp_cols[["R"]]

y1_max <- 2400                     # left axis: S and I
r_lo   <- 34500; r_hi <- 37500     # right axis: R
resc_R <- function(x) (x - r_lo) / (r_hi - r_lo) * y1_max

traj_long <- traj |>
  pivot_longer(c(S, I, R), names_to = "compartment", values_to = "count") |>
  mutate(compartment = factor(compartment, levels = c("S", "I", "R")),
         yval = if_else(compartment == "R", resc_R(count), count))
stopifnot(max(traj_long$yval) <= y1_max, min(traj_long$yval) >= 0)

p_bot <- ggplot(traj_long, aes(week, yval,
                               group = interaction(.id, compartment),
                               colour = compartment)) +
  geom_line(alpha = 0.5, linewidth = 0.4) +
  scale_colour_manual(values = comp_cols, name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1.2))) +
  annotate("text", x = 29, y = 750, hjust = 0, vjust = 1,
           label = "50 of infinitely many:\nwe never see S, I or R",
           colour = "grey35", size = 2.9, lineheight = 1.05) +
  xsc +
  scale_y_continuous(
    name = "S, I",
    limits = c(0, y1_max),
    breaks = c(0, 500, 1000, 1500, 2000),
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    sec.axis = sec_axis(~ . / y1_max * (r_hi - r_lo) + r_lo, name = "R",
                        breaks = c(35000, 36000, 37000),
                        labels = scales::label_number(
                          scale_cut = scales::cut_short_scale())),
    expand = expansion(mult = c(0.02, 0.04))) +
  labs(x = "week", title = "unobserved: the entire epidemic") +
  theme_bw() + base_txt +
  theme(plot.title = element_text(size = 10.5, face = "bold", colour = "grey25",
                                  margin = margin(b = 3)),
        axis.title.y.right = element_text(colour = green),
        axis.text.y.right  = element_text(colour = green),
        legend.position = "right",
        legend.key.size = unit(0.85, "lines"),
        legend.text = element_text(size = 8.5),
        legend.margin = margin(0, 0, 0, 0),
        ## keep this title off the baseline of the panel above it
        plot.margin = margin(9, 6, 2, 6))

fig <- p_top / p_bot + plot_layout(heights = c(1, 1.1))

ggsave("figures/fig-data-vs-hidden.png", fig,
       width = 5.5, height = 3.4, dpi = 300, bg = "white")

cat("wrote figures/fig-data-vs-hidden.png\n")
