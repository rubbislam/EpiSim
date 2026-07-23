## fig-if2-cooling.R -----------------------------------------------------------
## Lecture B_mle_mif2: what `cooling.fraction.50` actually means.
##
## Replaces the algorithmic-parameter bullets in 5_IF_theory/main.qmd:146-152
## ("decay of perturbation, cooling fraction in 50 iterations" and the a^(2m/50)
## exponent) with a picture.
##
## THE TRAP this figure defuses: students read cooling.fraction.50 = 0.5 as
## "halve the perturbation every iteration". It means: at iteration 50 the
## random-walk sd is 0.5x its starting value, having decayed geometrically and
## continuously the whole way there.
##
## THE CLOCK IS PER CALL -- this is why the x axis says "iteration within one
## mif2 call" and not "mif2 iteration". Verified against pomp 6.4 source and a
## live run: mif2_internal takes `mifiter = .ndone + n`, .ndone defaults to 0L,
## and NO mif2 method (pomp, pfilterd_pomp, mif2d_pomp, data.frame) ever passes
## it. So the deck's chained call at 5_IF_theory/main.qmd:258-268
##     mf |> mif2(cooling.fraction.50 = 0.3) -> mf
## RESTARTS the cooling clock at m = 1: the object's rw.sd is still the original
## uncooled spec, both calls yield 51-row traces, and the perturbation runs
## 0.5000 x rw.sd at the end of call 1 -> 0.9762 x rw.sd at the start of call 2
## -> 0.3000 at its end. A sawtooth, not a continuous glide. The curve below is
## the truth for ONE call; do not let the x axis imply otherwise.
##
## The curve is deterministic -- no pomp run needed. The perturbation-kick cloud
## in panel B is random but seeded, so the PNG is byte-reproducible.
##
## Run from the project root (where 4_mle_mif2.Rproj lives):
##   Rscript figures/fig-if2-cooling.R
## Reads:  nothing (self-contained; the cooling law is checked against pomp)
## Writes: figures/fig-if2-cooling.png

library(tidyverse)
library(patchwork)

theme_set(theme_bw())
set.seed(20240715)

## ---- the cooling law --------------------------------------------------------
## pomp's mif2_cooling(type="geometric", fraction=a, ntimes=n) returns, at
## iteration m and observation nt,
##     alpha = a^((nt/ntimes + m - 1)/50)   <- multiplies rw.sd
##     gamma = alpha^2                      <- multiplies the VARIANCE
## At the end of iteration m (nt = ntimes) this is exactly a^(m/50).
## That a^2 is why the pseudocode carries the exponent 2m/50: it is written for
## the variance matrix V_n. rw.sd is an sd, so the number students set decays as
## a^(m/50) -- and at m = 50 it equals a, exactly.
cooling_sd <- function(m, a) a^(m/50)

## Verify against pomp's own internal implementation rather than trusting the
## algebra above. Non-exported, so guard it: the figure must still build if pomp
## renames its internals.
verify_against_pomp <- function() {
  ok <- tryCatch({
    cool <- getFromNamespace("mif2_cooling", "pomp")
    ntimes <- 42L                       # Consett 1948: 42 weekly observations
    m <- c(1, 10, 25, 50, 75, 100)
    worst <- 0
    for (a in c(0.1, 0.5, 0.9)) {
      f <- cool("geometric", fraction = a, ntimes = ntimes)
      pomp_alpha <- vapply(m, function(mm) f(nt = ntimes, m = mm)$alpha, numeric(1))
      worst <- max(worst, max(abs(pomp_alpha - cooling_sd(m, a))))
    }
    stopifnot(worst < 1e-12)
    message(sprintf("  cooling law matches pomp::mif2_cooling (max abs diff %.2e)", worst))
    TRUE
  }, error = function(e) {
    message("  NOTE: could not verify against pomp internals: ", conditionMessage(e))
    FALSE
  })
  invisible(ok)
}
verify_against_pomp()

## ---- course palette ---------------------------------------------------------
blue  <- "#0072B2"   # the method that works: the true cooling law
verm  <- "#D55E00"   # the naive approach / what doesn't work: the misreading
grey_dark  <- "grey45"
grey_light <- "grey62"

fracs <- c(0.1, 0.5, 0.9)
pal   <- c("0.1" = grey_dark, "0.5" = blue, "0.9" = grey_light)

M <- 100  # extend well past 50 so the continued decay is obvious

curves <- expand_grid(a = fracs, m = seq(0, M, by = 0.25)) |>
  mutate(sd = cooling_sd(m, a), lab = factor(a))

## the dots that ARE the definition: at iteration 50, sd = a, exactly
anchors <- tibble(a = fracs, m = 50, sd = cooling_sd(50, fracs), lab = factor(fracs))

## the misreading: "0.5 means halve it EVERY iteration"
misread <- tibble(m = seq(0, M, by = 0.25), sd = 0.5^m)

## terminal values, for the right-hand labels. The a = 0.1 curve ends at 0.010,
## right on top of the flatlined misreading curve, so nudge that one label clear.
ends <- tibble(a = fracs, m = M, sd = cooling_sd(M, fracs), lab = factor(fracs),
               ylab = c(0.048, cooling_sd(M, 0.5), cooling_sd(M, 0.9)))

## ---- panel A: the definition, made visual -----------------------------------
pA <- ggplot() +
  ## misreading first, so it sits under the real curves
  geom_line(data = misread, aes(m, sd),
            colour = verm, linewidth = 0.7, linetype = "22") +
  geom_vline(xintercept = 50, colour = "grey30", linewidth = 0.4, linetype = "42") +
  geom_line(data = curves, aes(m, sd, colour = lab, linewidth = lab)) +
  geom_point(data = anchors, aes(m, sd, colour = lab), size = 1.9) +
  ## series identity: which cooling.fraction.50 each curve belongs to, printed
  ## at the RIGHT END of the curve (m = 100), not at the iteration-50 anchor dot.
  geom_text(data = ends, aes(m + 1.5, ylab, colour = lab, label = format(a, nsmall = 1)),
            hjust = 0, vjust = 0.5, size = 2.9, fontface = "bold") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_linewidth_manual(values = c("0.1" = 0.8, "0.5" = 1.3, "0.9" = 0.8),
                         guide = "none") +
  ## the required annotation, pointing at the (50, 0.5) dot
  annotate("curve", x = 75, y = 0.585, xend = 52.5, yend = 0.515,
           curvature = 0.30, linewidth = 0.35, colour = "grey25",
           arrow = arrow(length = unit(0.05, "in"), type = "closed")) +
  annotate("text", x = 99, y = 0.60, hjust = 1, vjust = 0,
           label = "at iteration 50,\nsd = 0.5 x start",
           size = 2.9, colour = "grey15", lineheight = 0.95) +
  ## kill the misconception, in the free corner just above its own dead curve
  annotate("text", x = 7, y = 0.105, hjust = 0, vjust = 0.5,
           label = "the misreading:\n\"halve it every\niteration\"",
           size = 2.6, colour = verm, lineheight = 1.0) +
  annotate("text", x = 49, y = 1.04, hjust = 1, vjust = 0.5,
           label = "iteration 50", size = 2.7, colour = "grey30") +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.02, 0.09))) +
  scale_y_continuous(breaks = c(0, 0.1, 0.5, 0.9, 1),
                     labels = c("0", "0.1", "0.5", "0.9", "1"),
                     expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(xlim = c(0, M), ylim = c(0, 1.08), clip = "off") +
  labs(x = "iteration within one mif2 call",
       y = "random-walk sd,\nrelative to its starting value",
       subtitle = "cooling.fraction.50 = 0.1, 0.5, 0.9") +
  theme(
    axis.text    = element_text(size = 7.5),
    axis.title   = element_text(size = 9),
    axis.title.y = element_text(lineheight = 0.95),
    plot.subtitle = element_text(size = 8.5, colour = "grey15"),
    panel.grid.minor = element_blank(),
    plot.margin  = margin(4, 20, 3, 3)
  )

## ---- panel B: what cooling actually does to the particles --------------------
## The kick each particle receives at iteration m is drawn N(0, sd_m), with
## sd_m = a^(m/50) x (starting sd). This is the perturbation itself -- NOT a
## claim that the swarm converges. Envelope = +/- 2 sd_m, i.e. panel A's curve.
a_demo <- 0.5
J <- 40
kicks <- expand_grid(m = 1:M, j = 1:J) |>
  mutate(sd = cooling_sd(m, a_demo),
         kick = rnorm(n(), 0, sd))

env <- tibble(m = seq(0, M, by = 0.25)) |>
  mutate(sd = cooling_sd(m, a_demo), lo = -2 * sd, hi = 2 * sd)

pB <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(data = kicks, aes(m, kick), colour = "grey60",
             size = 0.22, alpha = 0.35) +
  geom_line(data = env, aes(m, hi), colour = blue, linewidth = 0.9) +
  geom_line(data = env, aes(m, lo), colour = blue, linewidth = 0.9) +
  geom_vline(xintercept = 50, colour = "grey30", linewidth = 0.4, linetype = "42") +
  scale_x_continuous(breaks = c(0, 50, 100), expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(breaks = c(-2, 0, 2), expand = expansion(mult = c(0.04, 0.04))) +
  coord_cartesian(ylim = c(-3.1, 3.1), clip = "off") +
  labs(x = "iteration within one mif2 call",
       y = "kick given to each particle,\nin units of the starting sd",
       subtitle = "the kicks shrink with it") +
  theme(
    axis.text    = element_text(size = 7.5),
    axis.title   = element_text(size = 9),
    axis.title.y = element_text(lineheight = 0.95),
    plot.subtitle = element_text(size = 8.5, colour = "grey15"),
    panel.grid.minor = element_blank(),
    plot.margin  = margin(4, 3, 3, 6)
  )

## ---- assemble ---------------------------------------------------------------
fig <- pA + pB + plot_layout(widths = c(1.85, 1))

ggsave("figures/fig-if2-cooling.png", fig,
       width = 5.5, height = 3.0, dpi = 300, bg = "white")

## ---- the numbers this figure asserts ----------------------------------------
message("\n--- values plotted (sd relative to start) ---")
for (a in fracs) {
  message(sprintf("  a = %.1f : m=1 -> %.4f | m=50 -> %.4f | m=100 -> %.4f",
                  a, cooling_sd(1, a), cooling_sd(50, a), cooling_sd(100, a)))
}
message(sprintf("  misreading 0.5^m at m=10 -> %.3e ; at m=50 -> %.3e",
                0.5^10, 0.5^50))
message(sprintf("\nWrote figures/fig-if2-cooling.png (%.0f KB)",
                file.info("figures/fig-if2-cooling.png")$size / 1024))
