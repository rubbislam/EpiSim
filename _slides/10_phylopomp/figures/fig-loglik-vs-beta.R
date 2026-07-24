## fig-loglik-vs-beta.R ----------------------------------------------------
## 10_phylopomp (Lecture C) -- inference: the likelihood of an SIR model
## GIVEN a genealogy, evaluated across a grid of Beta. The tree was simulated
## at Beta = 3, and the likelihood peaks right there. This is the SAME
## particle filter, and the SAME logmeanexp, from Lectures A and B -- only
## the "data" has changed, from a case-count series to a genealogy.
##
## Run from the project root (where 10_phylopomp.Rproj lives):
##   Rscript figures/fig-loglik-vs-beta.R
## Reads/writes cache: results/loglik-vs-beta.rds
## Writes figure     : figures/fig-loglik-vs-beta.png
##
## The pfilter sweep takes ~1 minute; the result is cached, so the figure
## redraws instantly on a re-run. Delete the .rds to recompute.
## ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(phylopomp)
  library(pomp)
  library(tidyverse)
})

blue  <- "#0072B2"   # the method / the estimate
green <- "#009E73"   # the truth (Beta = 3)

## the same simulated genealogy as fig-sir-genealogy.R -------------------
set.seed(2024)
G <- runSIR(time = 4, Beta = 3, gamma = 1, psi = 0.4,
            pop = 500, S0 = 0.97, I0 = 0.03, R0 = 0)

beta_true <- 3
beta_grid <- seq(1.5, 5, by = 0.5)

cache <- "results/loglik-vs-beta.rds"
if (file.exists(cache)) {
  ll <- readRDS(cache)
} else {
  ll <- map_dfr(beta_grid, function(b) {
    set.seed(2024)                       # same filter randomness at each Beta
    G |>
      sir_pomp(Beta = b, gamma = 1, psi = 0.4,
               S0 = 0.97, I0 = 0.03, R0 = 0, pop = 500) |>
      pfilter(Np = 1000) |>
      replicate(n = 8) |>
      concat() -> pf
    est <- logmeanexp(sapply(pf, logLik), se = TRUE)
    tibble(Beta = b, loglik = est[1], se = est[2])
  })
  saveRDS(ll, cache)
}

mle <- ll |> slice_max(loglik, n = 1)

## plot -------------------------------------------------------------------
p <- ggplot(ll, aes(Beta, loglik)) +
  geom_vline(xintercept = beta_true, linetype = "dashed",
             colour = green, linewidth = 0.5) +
  annotate("text", x = beta_true, y = min(ll$loglik - ll$se), vjust = 0,
           hjust = -0.1, label = "truth", colour = green,
           size = 3.2, fontface = "bold") +
  geom_errorbar(aes(ymin = loglik - se, ymax = loglik + se),
                width = 0.08, colour = blue, alpha = 0.7) +
  geom_line(colour = blue, linewidth = 0.7) +
  geom_point(colour = blue, size = 2) +
  geom_point(data = mle, colour = "black", fill = "white",
             shape = 21, size = 2.6, stroke = 0.7) +
  labs(x = expression(beta),
       y = "log-likelihood of the genealogy",
       subtitle = "Same particle filter as Lectures A and B; the peak sits at the truth") +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 9, face = "bold"),
        panel.grid.minor = element_blank())

ggsave("figures/fig-loglik-vs-beta.png", p,
       width = 6.5, height = 3.6, dpi = 300, bg = "white")

## the numbers this figure claims ----------------------------------------
cat("\n---- fig-loglik-vs-beta: verified numbers ----\n")
cat(sprintf("phylopomp version : %s\n", packageVersion("phylopomp")))
cat(sprintf("MLE on the grid   : Beta=%.1f  loglik=%.2f (se %.2f)\n",
            mle$Beta, mle$loglik, mle$se))
cat(sprintf("at the truth Beta=3: loglik=%.2f (se %.2f)\n",
            ll$loglik[ll$Beta == 3], ll$se[ll$Beta == 3]))
cat("-------------------------------------------------\n")
