## fig-likelihood-ridge.R -------------------------------------------------
## B_mle_mif2 (Lecture B) -- "This isn't a hill. It's a knife-edge ridge,
## and it runs along constant R0."
##
## Replaces the prose at 4_pfilter_in_pomp/main.qmd:298-312 (peaks/valleys/
## ridges described in words) and :476-495 ("the surface is wedge-shaped, so
## its curvature varies considerably...").
##
## Run from the project root (where 4_mle_mif2.Rproj lives):
##   Rscript figures/fig-likelihood-ridge.R
## Reads : results/pfilter-grid1.rds
## Writes: figures/fig-likelihood-ridge.png
##
## Deterministic: the cache is read, averaged with logmeanexp, and
## interpolated. Nothing here is random; the seed is set for house-style
## conformity only.
## ------------------------------------------------------------------------

set.seed(1948)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pomp)
  library(patchwork)
})

theme_set(theme_bw())

## Okabe-Ito, fixed course meanings ---------------------------------------
blue  <- "#0072B2"  # the particle filter / the method that works
verm  <- "#D55E00"  # the naive approach / across the ridge / what doesn't work
green <- "#009E73"  # truth / the R0 ridge / reference

## ---- 1. the surface -----------------------------------------------------
## 40 x 40 Beta-Gamma grid, 9 particle-filter replicates per cell (Np=5000).
## logmeanexp averages the replicates on the *likelihood* scale, which is the
## right thing to do: the pfilter is unbiased for the likelihood, not for the
## log-likelihood. Note this is not itself an unbiased log-likelihood -- it is
## the log of an unbiased likelihood estimate, so Jensen still biases it low.
## Immaterial here (median within-cell replicate sd on the crest is 0.24, so
## the bias is ~0.01 log units against drops of 1-5 units), but this deck's
## whole point is precision on exactly this trap, so do not call it unbiased.

readRDS("results/pfilter-grid1.rds") |>
  group_by(Beta, Gamma) |>
  summarise(loglik = logmeanexp(loglik), .groups = "drop") -> surf

mle    <- surf |> slice_max(loglik, n = 1)
loglik_max <- mle$loglik

## the crest: best Gamma for each Beta. Beta >= 15 avoids the low-Beta corner
## where the ridge has not yet formed and the grid is coarse relative to it.
surf |>
  group_by(Beta) |>
  slice_max(loglik, n = 1) |>
  ungroup() |>
  mutate(R0 = Beta / Gamma) -> crest

R0_ridge <- mean(crest$R0[crest$Beta >= 15])
R0_sd    <- sd(crest$R0[crest$Beta >= 15])

## ---- 2. bilinear interpolator on the averaged surface --------------------
bs <- sort(unique(surf$Beta))
gs <- sort(unique(surf$Gamma))
M  <- matrix(NA_real_, length(bs), length(gs))
M[cbind(match(surf$Beta, bs), match(surf$Gamma, gs))] <- surf$loglik

bilin <- function(x, y) {
  i  <- pmin(findInterval(x, bs), length(bs) - 1)
  j  <- pmin(findInterval(y, gs), length(gs) - 1)
  tx <- (x - bs[i]) / (bs[i + 1] - bs[i])
  ty <- (y - gs[j]) / (gs[j + 1] - gs[j])
  M[cbind(i, j)]     * (1 - tx) * (1 - ty) +
  M[cbind(i + 1, j)] * tx       * (1 - ty) +
  M[cbind(i, j + 1)] * (1 - tx) * ty       +
  M[cbind(i + 1, j + 1)] * tx   * ty
}

## ---- 3. along the crest vs across it ------------------------------------
## Fair, scale-free comparison, both anchored at the MLE:
##   ALONG  : move Beta AND Gamma by the same % -> R0 unchanged.
##   ACROSS : move Gamma alone by that % -> R0 changes by that %.
## +/-30% keeps both walks strictly inside the grid, so no extrapolation.

pct <- seq(-30, 30, length.out = 241)

walk <- bind_rows(
  tibble(pct = pct, dir = "along",
         loglik = bilin(mle$Beta * (1 + pct/100), mle$Gamma * (1 + pct/100))),
  tibble(pct = pct, dir = "across",
         loglik = bilin(rep(mle$Beta, length(pct)), mle$Gamma * (1 + pct/100)))
)

## How far can you walk before losing 2 log units? Searched out to the grid
## edge (further than the +/-30% plotted), so it is never truncated by the plot
## window; no extrapolation beyond the cached grid.
##
## BOTH directions, always. The surface is asymmetric -- reporting only the
## increasing side makes the ridge look longer and flatter than it is (+30.6%
## vs -21.9% along), and a one-sided ratio quietly overstates the punchline.
## sgn = +1 increases the parameter(s), sgn = -1 decreases them.
reach <- function(d, sgn = 1) {
  lim <- if (d == "along") {
    if (sgn > 0) (max(bs)/mle$Beta - 1) else (1 - min(bs)/mle$Beta)
  } else {
    if (sgn > 0) (max(gs)/mle$Gamma - 1) else (1 - min(gs)/mle$Gamma)
  }
  e <- seq(0, lim, length.out = 4000)
  v <- if (d == "along") bilin(mle$Beta * (1 + sgn*e), mle$Gamma * (1 + sgn*e))
       else              bilin(rep(mle$Beta, length(e)), mle$Gamma * (1 + sgn*e))
  i <- which(loglik_max - v >= 2)[1]
  if (is.na(i)) NA_real_ else sgn * e[i] * 100
}

## ---- 4. left panel: the surface -----------------------------------------
## Clip 25 log units below the max -- the deck's existing idiom.
surf |>
  mutate(loglik = ifelse(loglik > loglik_max - 25, loglik, NA_real_)) -> surf_clip

## the constant-R0 line, drawn only where it is inside the panel
R0_line <- tibble(Beta = seq(min(bs), max(bs), length.out = 400)) |>
  mutate(Gamma = Beta / R0_ridge) |>
  filter(Gamma >= min(gs), Gamma <= max(gs))

## the "across" line: hold Beta at the MLE and vary Gamma. This is the vertical
## slice the orange curve in the right panel walks along -- stepping off the
## crest at fixed Beta, so R0 = Beta/Gamma changes.
across_line <- tibble(Beta = mle$Beta, Gamma = c(min(gs), max(gs)))

## The colour ramp is stretched towards the top: nearly all of the story
## lives in the few log units below the max, so a linear ramp over 25 units
## saturates the crest into a featureless blob and hides the knife edge.
ramp_at <- loglik_max - c(25, 15, 9, 5, 2, 0)

p_surf <- ggplot(surf_clip, aes(Beta, Gamma)) +
  geom_raster(aes(fill = loglik), interpolate = TRUE) +
  ## Contours read the UNCLIPPED surface: the clip is a fill decision, not a
  ## data decision. Feeding the clipped surface here NA'd 515 cells, which made
  ## the contourer work on a ragged grid and emitted a "removed 515 rows
  ## containing non-finite values" warning on every run. Every break is above
  ## max-25, so no line is drawn in the pale region either way -- this only
  ## keeps the deepest contour honest at the clip boundary, and keeps the run
  ## silent so a real warning would be noticed.
  geom_contour(data = surf, aes(z = loglik),
               breaks = loglik_max - c(1, 3, 6, 10, 15, 20),
               colour = "white", alpha = 0.55, linewidth = 0.3) +
  ## green R0 line, with a white halo so it separates from the blue field
  geom_line(data = R0_line, colour = "white", linewidth = 2.0, alpha = 0.95) +
  geom_line(data = R0_line, colour = green, linewidth = 0.95) +
  ## orange "across" line (fixed Beta = MLE), same white-halo style so it pairs
  ## one-to-one with the orange curve on the right panel
  geom_line(data = across_line, colour = "white", linewidth = 2.0, alpha = 0.95) +
  geom_line(data = across_line, colour = verm, linewidth = 0.95) +
  ## the MLE
  geom_point(data = mle, colour = "black", fill = "white",
             shape = 21, size = 2.0, stroke = 0.6) +
  ## white label box: the crest is dark, so green-on-blue would not read.
  ## angle 38 matches the line's on-page slope at 5.5 x 3.0 in.
  annotate("label", x = 14.6, y = 14.6/R0_ridge + 0.055,
           label = "constant~R[0] == 21.1", parse = TRUE, colour = green,
           size = 2.7, fontface = "bold", angle = 38,
           fill = "white", alpha = 0.88, border.colour = NA,
           label.padding = unit(1.1, "pt")) +
  ## orange label near the top of the vertical line (a pale region), matching
  ## the green label's white-box style
  annotate("label", x = mle$Beta, y = max(gs) - 0.10,
           label = "vary~gamma~alone", parse = TRUE, colour = verm,
           size = 2.7, fontface = "bold",
           fill = "white", alpha = 0.88, border.colour = NA,
           label.padding = unit(1.1, "pt")) +
  annotate("text", x = mle$Beta + 0.6, y = mle$Gamma - 0.075, label = "MLE",
           colour = "white", size = 2.7, hjust = 0, fontface = "bold") +
  ## this annotation replaces both the colourbar and the deck's prose
  ## ("all points ... less than 25 units below the maximum are shown in grey")
  annotate("text", x = 29.4, y = 0.455, label = "pale: > 25 log units\nbelow the maximum",
           colour = "grey45", size = 2.6, hjust = 1, vjust = 0, lineheight = 0.9) +
  scale_fill_gradientn(colours = c("grey95", "#cfe0ee", "#8fbcdb",
                                   "#3f8fc0", blue, "#004c78"),
                       values = scales::rescale(ramp_at),
                       na.value = "grey97", guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = range(bs), ylim = range(gs)) +
  labs(x = expression(beta), y = expression(gamma),
       subtitle = "Darker = better fit. The good region is a ridge.") +
  theme(
    plot.subtitle   = element_text(size = 8.2, face = "bold"),
    axis.text       = element_text(size = 7.2),
    axis.title      = element_text(size = 9.5),
    panel.grid      = element_blank(),
    plot.margin     = margin(2, 3, 2, 2)
  )

## ---- 5. right panel: flat along, steep across ---------------------------
p_walk <- ggplot(walk, aes(pct, loglik, colour = dir)) +
  geom_hline(yintercept = loglik_max - 2, linetype = "dashed",
             colour = "grey55", linewidth = 0.3) +
  geom_line(linewidth = 1) +
  geom_point(data = tibble(pct = 0, loglik = loglik_max, dir = "along"),
             colour = "black", size = 1.3) +
  ## labels sit in the empty corners: green above its own curve at top-right,
  ## vermillion in the dead space at bottom-centre.
  annotate("text", x = 30, y = -105.9, label = "along the crest",
           colour = green, size = 2.7, hjust = 1, vjust = 1, fontface = "bold") +
  annotate("text", x = 30, y = -107.0, label = "(R[0]~fixed)", parse = TRUE,
           colour = green, size = 2.7, hjust = 1, vjust = 1) +
  annotate("text", x = 0, y = -117.9, label = "across the crest",
           colour = verm, size = 2.7, hjust = 0.5, vjust = 0, fontface = "bold") +
  annotate("text", x = 0, y = -119.1, label = "(gamma~alone)", parse = TRUE,
           colour = verm, size = 2.7, hjust = 0.5, vjust = 0) +
  annotate("text", x = 30, y = -110.2, label = "2 log units down",
           colour = "grey40", size = 2.6, hjust = 1, vjust = 1) +
  scale_colour_manual(values = c(along = green, across = verm), guide = "none") +
  scale_x_continuous(breaks = seq(-30, 30, by = 15),
                     labels = function(x) paste0(x, "%")) +
  coord_cartesian(xlim = c(-30, 30), ylim = c(-119.2, -105.8)) +
  labs(x = "% change from the MLE",
       y = "log-likelihood",
       subtitle = "Flat along it, steep across it") +
  theme(
    plot.subtitle = element_text(size = 8.2, face = "bold"),
    axis.text     = element_text(size = 7.2),
    axis.title    = element_text(size = 9.5),
    panel.grid.minor = element_blank(),
    plot.margin   = margin(2, 2, 2, 4)
  )

## ---- 6. assemble --------------------------------------------------------
p <- p_surf + p_walk + plot_layout(widths = c(1.2, 1))

ggsave("figures/fig-likelihood-ridge.png", p,
       width = 5.5, height = 3.0, dpi = 300, bg = "white")

## ---- 7. the numbers this figure claims ----------------------------------
cat("\n---- fig-likelihood-ridge: verified numbers ----\n")
cat(sprintf("MLE            : Beta=%.2f Gamma=%.3f loglik=%.2f  (R0=%.2f)\n",
            mle$Beta, mle$Gamma, loglik_max, mle$Beta/mle$Gamma))
cat(sprintf("ridge R0       : %.2f +/- %.2f  (Beta >= 15)\n", R0_ridge, R0_sd))
cat(sprintf("crest Beta 15->30 at R0=%.2f: loglik falls %.2f units total\n",
            R0_ridge, diff(range(bilin(seq(15, 30, length.out = 200),
                                       seq(15, 30, length.out = 200)/R0_ridge)))))
for (d in c("along", "across")) {
  w <- walk |> filter(dir == d)
  cat(sprintf("%-7s: drop at -20%%=%.2f  +20%%=%.2f | 2 log units lost at %+.1f%% / %+.1f%%\n",
              d,
              w$loglik[which.min(abs(w$pct + 20))] - loglik_max,
              w$loglik[which.min(abs(w$pct - 20))] - loglik_max,
              reach(d, -1), reach(d, +1)))
}

## The asymmetry is the whole reason these are printed two-sided. Say the
## defensible range out loud so nobody lifts a single flattering number.
r_dn <- reach("along", -1) / reach("across", -1)
r_up <- reach("along", +1) / reach("across", +1)
cat(sprintf("2-log reach is %.1fx-%.1fx further along than across (down / up)\n",
            min(r_dn, r_up), max(r_dn, r_up)))
w_al <- walk |> filter(dir == "along"); w_ac <- walk |> filter(dir == "across")
ratio20 <- c(
  (w_ac$loglik[which.min(abs(w_ac$pct + 20))] - loglik_max) /
  (w_al$loglik[which.min(abs(w_al$pct + 20))] - loglik_max),
  (w_ac$loglik[which.min(abs(w_ac$pct - 20))] - loglik_max) /
  (w_al$loglik[which.min(abs(w_al$pct - 20))] - loglik_max))
cat(sprintf("across is %.1fx-%.1fx steeper than along at a +/-20%% move\n",
            min(ratio20), max(ratio20)))
cat("SPOKEN LINE: a 20% move along the crest costs UNDER TWO log units\n")
cat("             (-1.69 down, -0.92 up). Do NOT say 'under one' -- false at -20%.\n")
cat("-------------------------------------------------\n")
