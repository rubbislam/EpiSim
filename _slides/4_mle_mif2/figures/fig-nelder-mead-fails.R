## fig-nelder-mead-fails.R -----------------------------------------------------
## Lecture A: ATTEMPT 2 -- throw the particle-filter likelihood at
## optim(method="Nelder-Mead").  It fails.  This figure SHOWS the failure and
## shows WHY.
##
## THE MECHANISM (read this before editing the annotations -- an earlier version
## of this figure told the story BACKWARDS):
##   The optimizer does NOT wander along the flat crest.  It never reaches the
##   crest.  It dies on the RISING FLANK, ~300 evaluations in, still climbing.
##   10 of the 12 runs stop below Beta=19; the true maximum is at Beta=22.1.
##   The 2 runs that DO reach the crest are the 2 SUCCESSES (honest -107.5,
##   -107.6, within 0.11 of the top).  Reaching the crest = winning.
##   The failure is dying on the climb, and the reason is the noise-to-slope
##   ratio DOWN THERE: at Beta=15 one Np=1000 call has a 95% spread of ~5.7 log
##   units, while the hill only rises ~1.2 log units per unit of Beta.  A single
##   objective call at the start is worth +/- 5 units of Beta of pure noise.
##   That is why the simplex cannot find the uphill direction.  Noise at the
##   crest (~1.4) is ~4x SMALLER -- measuring it there understates the case.
##
##   "The optimizer isn't broken.  The function you handed it is."
##
## Run from the project root (where 4_mle_mif2.Rproj lives):
##     Rscript figures/fig-nelder-mead-fails.R
##
## Reads : data/ (via scripts/model_measSIR.R), results/pfilter-grid1.rds
## Writes: results/nelder-mead-fails.rds   (cache -- the deck never recomputes)
##         figures/fig-nelder-mead-fails.png
## Set RECOMPUTE <- TRUE to force the ~70 s of particle filtering to re-run.
## Seeds are set inside every worker, so results do not depend on worker count.
## -----------------------------------------------------------------------------

library(tidyverse)
library(pomp)
library(patchwork)
library(foreach)
library(doParallel)

RECOMPUTE <- FALSE
CACHE     <- "results/nelder-mead-fails.rds"
NWORKERS  <- 12L      # results are worker-count independent; lower this freely

## ---- house style ------------------------------------------------------------
## Okabe-Ito, course-fixed meanings.  #E69F00 amber = the data: it does not
## appear here -- this figure is about the naive optimizer meeting the truth.
## Vermillion is used for BOTH the Np=1000 objective and the optimizer's
## answers, which is the same thread: "the naive approach / what doesn't work".
## (fig-loglik-is-random.R already maps Np=1000 -> vermillion.)
theme_set(theme_bw())
verm  <- "#D55E00"   # the naive approach / what doesn't work / Np = 1000
green <- "#009E73"   # truth / the R0 ridge / reference
grey  <- "grey60"    # de-emphasised context

## ---- the ridge, defined by the cached 40x40 grid -----------------------------
## pfilter-grid1.rds: 40 Beta x 40 Gamma, 9 pfilter reps per cell at Np=5000.
grid <- readRDS("results/pfilter-grid1.rds")
attr(grid, "rng") <- NULL; attr(grid, "out.attrs") <- NULL
cells <- as.data.frame(grid) |>
  group_by(Beta, Gamma) |>
  summarise(loglik = logmeanexp(loglik), .groups = "drop")
grid_mle <- cells |> slice_max(loglik, n = 1)
R0_ridge <- grid_mle$Beta / grid_mle$Gamma          # 21.27
cat(sprintf("cached grid MLE: Beta=%.3f Gamma=%.4f loglik=%.2f  -> R0 ridge=%.3f\n",
            grid_mle$Beta, grid_mle$Gamma, grid_mle$loglik, R0_ridge))

## ---- the experiment ---------------------------------------------------------
if (RECOMPUTE || !file.exists(CACHE)) {

  source("scripts/model_measSIR.R")
  registerDoParallel(NWORKERS)
  wall0 <- Sys.time()

  fixed_pars <- c(Rho = 0.5, k = 10, Eta = 0.06, N = 38000)
  ## unname(): b and g often arrive already named (e.g. exp(par[1]) inside the
  ## optimizer), and c(Beta = c(Beta = 15)) would silently build "Beta.Beta".
  pars   <- function(b, g) c(Beta = unname(b), Gamma = unname(g), fixed_pars)
  betas  <- seq(14, 26, length.out = 41)

  ## (1) TRUTH along the crest of the ridge: Gamma follows Beta so that
  ##     R0 = Beta/Gamma stays at the ridge value.  Np=10000 x 36 reps -- heavy
  ##     on purpose, so that the green curve really is smooth and any wiggle the
  ##     audience sees in the vermillion curve is the objective, not my estimate.
  truth <- foreach(i = seq_along(betas), .combine = rbind, .packages = "pomp") %dopar% {
    set.seed(20260715L + i)
    b  <- betas[i]
    ll <- replicate(36, logLik(pfilter(measSIR, params = pars(b, b / R0_ridge), Np = 10000)))
    data.frame(Beta = b, loglik = logmeanexp(ll), se = logmeanexp(ll, se = TRUE)[2])
  }

  ## (2) WHAT THE OPTIMIZER SEES: the same transect, but ONE particle filter
  ##     run at Np=1000 with a fresh seed at every Beta -- one honest call of
  ##     the objective function, exactly as optim() would make it.
  jagged <- foreach(i = seq_along(betas), .combine = rbind, .packages = "pomp") %dopar% {
    set.seed(31337L + i)
    b <- betas[i]
    data.frame(Beta = b, loglik = logLik(pfilter(measSIR, params = pars(b, b / R0_ridge), Np = 1000)))
  }

  ## (3) ATTEMPT 2: Nelder-Mead on that objective.  Same start, same code,
  ##     same settings -- only the random seed differs between runs.
  ##     Beta and Gamma are estimated on the log scale; Eta is held at 0.06
  ##     so that the cached grid is the right reference.
  objective <- function(par) {
    ll <- tryCatch(logLik(pfilter(measSIR, params = pars(exp(par[1]), exp(par[2])), Np = 1000)),
                   error = function(e) NA_real_)
    if (!is.finite(ll)) 1e10 else -ll
  }
  start <- c(Beta = 15, Gamma = 0.5)
  nm <- foreach(s = 1:12, .combine = rbind, .packages = "pomp") %dopar% {
    set.seed(1000L + s)
    f <- optim(log(start), objective, method = "Nelder-Mead", control = list(maxit = 1000))
    ## Honest value of the likelihood AT the answer the optimizer handed back
    ## (Np=10000 x 6 reps) -- the optimizer never gets to see this.
    set.seed(99000L + s)
    hh <- replicate(6, logLik(pfilter(measSIR,
                                      params = pars(exp(f$par[1]), exp(f$par[2])), Np = 10000)))
    data.frame(run = s, Beta = exp(f$par[1]), Gamma = exp(f$par[2]),
               reported = -f$value, honest = logmeanexp(hh),
               nfev = unname(f$counts[1]), convergence = f$convergence)
  }

  ## (4) Monte Carlo noise of ONE objective call -- measured at THREE places,
  ##     crucially INCLUDING where the runs actually die.  Measuring only at the
  ##     crest (the old version did) understates the noise ~4x and describes a
  ##     region no failing run ever visited.
  ##       Beta = 15   : the start, and where most runs stop
  ##       Beta = 18   : the top of the pile of failures
  ##       Beta = 22.4 : the crest, near the true maximum
  noise_betas <- c(15, 18, 22.4)
  ngrid <- expand.grid(rep = 1:240, bi = seq_along(noise_betas))
  noise <- foreach(r = seq_len(nrow(ngrid)), .combine = rbind, .packages = "pomp") %dopar% {
    j <- ngrid$bi[r]; i <- ngrid$rep[r]
    set.seed(5000L + 1000L * j + i)
    b <- noise_betas[j]
    data.frame(Beta = b, loglik = logLik(pfilter(measSIR, params = pars(b, b / R0_ridge), Np = 1000)))
  }

  wall <- as.numeric(difftime(Sys.time(), wall0, units = "secs"))
  cat(sprintf("compute wall-clock: %.1f s on %d workers\n", wall, NWORKERS))

  saveRDS(list(truth = truth, jagged = jagged, nm = nm, noise = noise,
               R0_ridge = R0_ridge, grid_mle = grid_mle, start = start,
               noise_betas = noise_betas, wall_secs = wall),
          CACHE)
}

res <- readRDS(CACHE)
truth <- res$truth; jagged <- res$jagged; nm <- res$nm

## ---- the numbers this figure is allowed to state ----------------------------
best <- max(truth$loglik)
best_beta <- truth$Beta[which.max(truth$loglik)]

## Noise of ONE Np=1000 call, by Beta.  This is the headline evidence.
nsum <- res$noise |>
  group_by(Beta) |>
  summarise(sd     = sd(loglik),
            lo     = unname(quantile(loglik, 0.025)),
            hi     = unname(quantile(loglik, 0.975)),
            spread = hi - lo, .groups = "drop")

## SIGNAL the optimizer is trying to climb: the slope of the true curve on the
## flank where the runs actually die.  Compare this with the noise above.
flank      <- truth |> filter(Beta >= 14, Beta <= 19)
flank_slope <- unname(coef(lm(loglik ~ Beta, data = flank))[2])

## Secondary (true, but NOT the headline): the crest really is flat.  This is
## why the 2 successful runs land at 22.6/23.1 rather than exactly at 22.1.
crest    <- truth |> filter(Beta >= 21, Beta <= 25.4)
crest_dl <- max(crest$loglik) - min(crest$loglik)

died_low  <- nm |> filter(Beta < 19)
reached   <- nm |> filter(Beta >= 19)
n_low     <- nrow(died_low)

cat(sprintf("\n-- true maximum: loglik %.2f at Beta %.1f\n", best, best_beta))
cat(sprintf("-- SIGNAL: true slope on the flank (Beta 14-19) = %.2f loglik per unit Beta\n",
            flank_slope))
cat("-- NOISE: one Np=1000 call (240 reps each):\n")
for (i in seq_len(nrow(nsum)))
  cat(sprintf("     Beta %5.1f : sd %.3f, 95%% spread %.2f log units  (= %.1f units of Beta)\n",
              nsum$Beta[i], nsum$sd[i], nsum$spread[i], nsum$spread[i] / flank_slope))
cat(sprintf("   -> noise at the start is %.1fx the noise at the crest\n",
            nsum$spread[nsum$Beta == 15] / nsum$spread[nsum$Beta == 22.4]))
cat(sprintf("-- WHERE THEY STOPPED: %d of %d below Beta=19 (%.2f to %.2f); only %d reached the crest\n",
            n_low, nrow(nm), min(died_low$Beta), max(died_low$Beta), nrow(reached)))
cat(sprintf("   the %d that died low are %.2f to %.2f log units short (median %.2f)\n",
            n_low, min(best - died_low$honest), max(best - died_low$honest),
            median(best - died_low$honest)))
cat(sprintf("   the %d that reached the crest are %.2f and %.2f short -- they WON\n",
            nrow(reached), (best - reached$honest)[1], (best - reached$honest)[2]))
cat(sprintf("-- all %d overstate their own answer, by %.2f to %.2f (mean %.2f);\n",
            sum(nm$reported > nm$honest), min(nm$reported - nm$honest),
            max(nm$reported - nm$honest), mean(nm$reported - nm$honest)))
cat(sprintf("   selection bias tracks the noise: %.2f at Beta=%.1f -> %.2f at Beta=%.1f\n",
            (nm$reported - nm$honest)[which.min(nm$Beta)], min(nm$Beta),
            (nm$reported - nm$honest)[which.max(nm$Beta)], max(nm$Beta)))
cat(sprintf("-- secondary: crest (Beta 21.0-25.4) is flat, true loglik moves only %.2f units\n",
            crest_dl))
cat(sprintf("-- R0 = Beta/Gamma is the robust part: %.2f to %.2f (ridge %.2f)\n",
            min(nm$Beta / nm$Gamma), max(nm$Beta / nm$Gamma), R0_ridge))
cat(sprintf("-- convergence codes: %s (10 = simplex degenerated, NOT converged); evals %d-%d\n",
            paste(sort(unique(nm$convergence)), collapse = ","), min(nm$nfev), max(nm$nfev)))

## ---- shared scales ----------------------------------------------------------
## Same x and same y on both panels, so B reads against the curve from A.
## NB the dots in B need NOT lie ON that curve: the green curve is the transect
## along the R0 ridge (Gamma = Beta/R0_ridge), whereas each dot is the honest
## loglik at that run's OWN (Beta, Gamma).  A dot sitting below the curve is a
## run that also drifted off the ridge -- e.g. run 4 (Beta 15.0, R0 22.9).
## ylim reaches -118.6 to fit the Beta=15 noise bar; the crest is flat anyway.
ylim <- c(-118.6, -106.3)
xlim <- c(13.7, 26.3)
ax <- theme(
  plot.title  = element_text(size = 9.5, face = "bold", hjust = 0, margin = margin(b = 3)),
  axis.title  = element_text(size = 9),
  axis.text   = element_text(size = 7.5),
  plot.margin = margin(3, 6, 2, 3)
)

## ---- PANEL A: the function you hand the optimizer ---------------------------
## The bars are the headline: one Np=1000 call scatters by ~5.7 log units at the
## start and only ~1.4 at the crest, against a hill that rises ~1.2 per unit Beta.
bars <- nsum |> filter(Beta %in% c(15, 22.4))
pA <- ggplot() +
  geom_line(data = truth, aes(Beta, loglik), colour = green, linewidth = 1.2) +
  geom_line(data = jagged, aes(Beta, loglik), colour = verm, linewidth = 0.4, alpha = 0.55) +
  geom_point(data = jagged, aes(Beta, loglik), colour = verm, size = 0.6, alpha = 0.55) +
  geom_errorbar(data = bars, aes(x = Beta, ymin = lo, ymax = hi),
                colour = verm, width = 0.55, linewidth = 0.85) +
  annotate("text", x = 14.1, y = -106.9, label = "the truth", colour = green,
           size = 2.9, hjust = 0, fontface = "bold") +
  annotate("text", x = 16.0, y = -117.4, hjust = 0, size = 2.7, colour = verm,
           lineheight = 1.05,
           label = sprintf("one call here scatters\nby %.1f log units",
                           nsum$spread[nsum$Beta == 15])) +
  annotate("text", x = 25.9, y = -109.8, hjust = 1, size = 2.7, colour = verm,
           lineheight = 1.05,
           label = sprintf("up here, only %.1f", nsum$spread[nsum$Beta == 22.4])) +
  ## The slope is a property of the green curve, so it is labelled green: the
  ## contrast the audience must see is verm noise (5.5) vs green signal (1.2).
  annotate("text", x = 20.4, y = -114.2, hjust = 0, size = 2.6, colour = green,
           lineheight = 1.05,
           label = sprintf("the hill only rises\n%.1f per unit of beta", flank_slope)) +
  coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "The climb is buried in noise",
       x = expression(beta), y = "log likelihood") +
  ax

## ---- PANEL B: 12 runs, and they die on the way up ---------------------------
pB <- ggplot() +
  geom_hline(yintercept = best, colour = green, linetype = "dashed", linewidth = 0.4) +
  geom_line(data = truth, aes(Beta, loglik), colour = green, linewidth = 1.2) +
  geom_vline(xintercept = unname(res$start["Beta"]), colour = grey,
             linetype = "dashed", linewidth = 0.4) +
  ## shape 21 + white stroke: the dots at Beta 17.9/18.2/18.6 sit on top of one
  ## another, and the audience must be able to count 12 of them.
  geom_point(data = nm, aes(Beta, honest), shape = 21, fill = verm,
             colour = "white", stroke = 0.3, size = 2.2) +
  annotate("text", x = 14.1, y = -106.9, label = "the maximum", colour = green,
           size = 2.75, hjust = 0, fontface = "bold") +
  annotate("text", x = 15.3, y = -117.9, label = "all 12 runs start here",
           colour = grey, size = 2.6, hjust = 0) +
  annotate("text", x = 25.9, y = -112.8, hjust = 1, size = 2.7, colour = verm,
           lineheight = 1.05,
           label = sprintf("%d of %d never got\npast beta 19 --\nthey quit on the climb",
                           n_low, nrow(nm))) +
  annotate("text", x = 21.4, y = -108.6, hjust = 0.5, size = 2.6, colour = green,
           lineheight = 1.05, label = "the 2 that got\nhere were right") +
  coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "Where 12 runs stopped",
       x = expression(beta), y = NULL) +
  ax + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p <- pA + pB + plot_layout(widths = c(1, 1))

ggsave("figures/fig-nelder-mead-fails.png", p,
       width = 5.5, height = 3.0, dpi = 300, bg = "white")
cat("\nwrote figures/fig-nelder-mead-fails.png\n")
