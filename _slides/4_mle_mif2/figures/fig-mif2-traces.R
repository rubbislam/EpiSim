## fig-mif2-traces.R ------------------------------------------------------
## B_mle_mif2 (Lecture B) headline evidence figure: IF2 convergence traces.
##
## Run from the project root (where 4_mle_mif2.Rproj lives):
##     Rscript figures/fig-mif2-traces.R
##
## Reads : results/local_search.rds  (mif2List, 20 replicated mif2 runs)
##         results/lik_local.rds     (honest pfilter re-evaluation of the 20 endpoints)
## Writes: figures/fig-mif2-traces.png
##
## THE POINT: the loglik panel converges (all 20 runs agree the fit is good);
## the individual parameter panels do NOT (Beta lands anywhere from 16 to 23).
## Both are true at once. The last panel says why: Beta and Eta trade off along
## a ridge, and the product Beta*Eta is what the data actually pins down.
##
## Every number quoted above is recomputed from the cache at run time and
## printed to the console -- do not hand-edit them. The cache MUST come from
## the same model as scripts/model_measSIR.R (Gamma=0.5). The legacy
## mif/results/local_search.rds is a DIFFERENT model (Gamma=2) and will
## silently produce a figure that contradicts the exercises.
## -------------------------------------------------------------------------

suppressMessages({
  library(pomp)
  library(tidyverse)
  library(patchwork)
})

set.seed(1350254336)   # nothing here is random; fixed for byte-reproducibility

## ---- Okabe-Ito course palette (fixed meanings) --------------------------
BLUE  <- "#0072B2"  # the method that works -- IF2 traces
VERM  <- "#D55E00"  # what doesn't work -- the params that never agree
GREEN <- "#009E73"  # truth / reference / the ridge-invariant quantity

## ---- Load ---------------------------------------------------------------
mifs <- readRDS("results/local_search.rds")
liks <- readRDS("results/lik_local.rds")

stopifnot(is(mifs, "mif2List"))

n_rep  <- length(mifs)
n_mif  <- mifs[[1]]@Nmif
n_part <- unique(as.integer(mifs[[1]]@Np))[1]
cool   <- mifs[[1]]@cooling.fraction.50

## ---- Traces -> long ------------------------------------------------------
## pomp stores rows 0..Nmif. Row i holds the parameters entering iteration i
## and the loglik of that iteration's (perturbed) filtering pass; the final
## row therefore has NA loglik -- which is exactly why lik_local.rds exists.
tr <- pomp::melt(pomp::traces(mifs)) |>
  as_tibble() |>
  rename(rep = .L1) |>
  mutate(iteration = iteration - 1L)   # -> 0..Nmif

wide <- tr |>
  pivot_wider(names_from = name, values_from = value) |>
  mutate(BetaEta = Beta * Eta)

## ---- Real numbers, computed (never typed from memory) --------------------
fin  <- wide |> filter(iteration == n_mif)          # final parameter values
ll0  <- wide |> filter(iteration == 0)              # starting filter loglik

starts <- wide |> filter(iteration == 0) |> summarise(
  Beta = unique(round(Beta, 6)) |> paste(collapse = "/"),
  Rho  = unique(round(Rho, 6))  |> paste(collapse = "/"),
  Eta  = unique(round(Eta, 6))  |> paste(collapse = "/")
)

rng <- function(x) c(min(x), max(x))

## Gamma is NOT estimated here -- mif2's rw.sd covers Beta, Rho, Eta only, so
## Gamma is pinned at whatever the pomp object carried in. Read it off the cache
## rather than typing it: if local_search.rds is ever regenerated at a different
## Gamma, the card below re-labels itself instead of lying.
gamma_fix <- unique(wide$Gamma)
stopifnot(length(gamma_fix) == 1L)

ll_spread_0   <- diff(rng(ll0$loglik))
ll_fin_rng    <- rng(liks$loglik)          # honest pfilter logliks
ll_fin_spread <- diff(ll_fin_rng)
beta_rng      <- rng(fin$Beta)
rho_rng       <- rng(fin$Rho)
eta_rng       <- rng(fin$Eta)
be_rng        <- rng(fin$BetaEta)
be_cv         <- sd(fin$BetaEta) / mean(fin$BetaEta)
beta_cv       <- sd(fin$Beta) / mean(fin$Beta)
r_beta_eta    <- cor(fin$Beta, fin$Eta)

cat("\n---- fig-mif2-traces: computed numbers ----\n")
cat(sprintf("replicates            : %d\n", n_rep))
cat(sprintf("mif2 iterations       : %d (Np=%d, cooling.fraction.50=%.2f)\n",
            n_mif, n_part, cool))
cat(sprintf("starting point (all)  : Beta=%s Rho=%s Eta=%s\n",
            starts$Beta, starts$Rho, starts$Eta))
cat(sprintf("Gamma (FIXED, not est): %g\n", gamma_fix))
cat(sprintf("loglik spread @iter 0 : %.1f units\n", ll_spread_0))
cat(sprintf("loglik final (pfilter): %.2f to %.2f  (spread %.2f)\n",
            ll_fin_rng[1], ll_fin_rng[2], ll_fin_spread))
cat(sprintf("Beta  final           : %.1f to %.1f   (CV %.1f%%)\n",
            beta_rng[1], beta_rng[2], 100 * beta_cv))
cat(sprintf("Rho   final           : %.3f to %.3f\n", rho_rng[1], rho_rng[2]))
cat(sprintf("Eta   final           : %.4f to %.4f\n", eta_rng[1], eta_rng[2]))
cat(sprintf("Beta*Eta final        : %.2f to %.2f  (CV %.1f%%)\n",
            be_rng[1], be_rng[2], 100 * be_cv))
cat(sprintf("cor(Beta, Eta)        : %.3f\n", r_beta_eta))
cat("-------------------------------------------\n\n")

## ---- Panel builder -------------------------------------------------------
base_thm <- theme_bw(base_size = 8) +
  theme(
    plot.title       = element_text(size = 8.5, face = "bold", hjust = 0,
                                    margin = margin(b = 1.5)),
    axis.text        = element_text(size = 7),
    axis.title       = element_text(size = 9),
    axis.title.y     = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    plot.margin      = margin(2, 4, 2, 2)
  )

trace_panel <- function(var, title, ylim = NULL, xlab = NULL, breaks = waiver()) {
  d <- wide |> filter(!is.na(.data[[var]]))
  p <- ggplot(d, aes(iteration, .data[[var]], group = rep)) +
    geom_line(colour = BLUE, alpha = 0.30, linewidth = 0.3) +
    scale_x_continuous(breaks = c(0, 25, 50), expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = breaks) +
    labs(title = title, x = xlab) +
    base_thm
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  if (is.null(xlab)) p <- p + theme(axis.title.x = element_blank())
  p
}

lab <- function(p, x, y, text, colour, hjust = 0, vjust = 1, size = 2.6) {
  p + annotate("text", x = x, y = y, label = text, colour = colour,
               hjust = hjust, vjust = vjust, size = size, lineheight = 0.95)
}

## ---- axis limits, derived from the data ---------------------------------
## These were once hardcoded, which silently CLIPPED the traces the day the
## underlying cache changed. Never hardcode them: derive from the data and
## reserve an explicit slot for the label so text can never sit on a trace.
##   lab_at="bottom"/"top" reserves lab_frac of the data span for the label.
auto_ylim <- function(var, lab_at = "none", lab_frac = 0.20, pad = 0.06) {
  r <- range(wide[[var]], na.rm = TRUE); s <- diff(r)
  lo <- r[1] - pad * s; hi <- r[2] + pad * s
  if (lab_at == "bottom") lo <- lo - lab_frac * s
  if (lab_at == "top")    hi <- hi + lab_frac * s
  c(lo, hi)
}
lab_y <- function(ylim, at = "bottom") {
  s <- diff(ylim)
  if (at == "bottom") ylim[1] + 0.02 * s else ylim[2] - 0.02 * s
}

## ---- Panel 1: loglik -- CONVERGES ---------------------------------------
## Show the FULL climb, not just the endgame. The story a first-time viewer
## needs is "it rockets up, then all twenty agree" -- and at full scale the
## converged band is ~1% of panel height, which reads as a single line and
## makes the agreement obvious on sight. The exact residual spread is small
## enough to be worth stating in words rather than pixels, so the green label
## carries it.
##
## DO NOT use `fin` here. `fin` is the final trace row and its loglik is NA by
## construction (mif2 never filters the final parameter vector -- that is the
## point of this whole lecture). min(fin$loglik) is NA, which silently turns
## coord_cartesian into a no-op and drops every annotation. Use the traces.
## Label wraps to TWO lines deliberately: as one line at size 2.6 it is wider
## than the panel and gets clipped at the left edge. Keep the wrap.
ll_ylim <- auto_ylim("loglik", lab_at = "bottom", lab_frac = 0.26)
p_ll <- trace_panel("loglik", "log likelihood", ylim = ll_ylim) |>
  lab(x = 50, y = lab_y(ll_ylim),
      text = sprintf("they agree:\nfinal spread %.1f units", ll_fin_spread),
      colour = GREEN, hjust = 1, vjust = 0, size = 2.6)

## ---- inset: a zoom of the converged tail --------------------------------
## At the full -150..-106 scale the last iterations collapse onto a single
## line, hiding the small run-to-run differences. This inset re-plots the tail
## on its own narrow y-scale, so the 20 runs read as 20 distinct lines. Derive
## the window from the data (never hardcode); it floats in the empty centre of
## the panel, below the converged band and above the green label.
tail_from <- 25L
d_tail    <- wide |> filter(!is.na(loglik), iteration >= tail_from)
tail_ylim <- range(d_tail$loglik) + c(-1, 1) * 0.10 * diff(range(d_tail$loglik))

inset_ll <- ggplot(d_tail, aes(iteration, loglik, group = rep)) +
  geom_line(colour = BLUE, alpha = 0.55, linewidth = 0.3) +
  coord_cartesian(xlim = c(tail_from, n_mif), ylim = tail_ylim, expand = FALSE) +
  scale_x_continuous(breaks = c(tail_from, n_mif)) +
  scale_y_continuous(n.breaks = 3) +
  theme_bw(base_size = 6) +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_text(size = 5.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.15),
    plot.background  = element_rect(fill = "white", colour = "grey45",
                                    linewidth = 0.4),
    plot.margin      = margin(1, 1, 1, 1)
  )

p_ll <- p_ll +
  patchwork::inset_element(inset_ll, left = 0.33, bottom = 0.24,
                           right = 0.985, top = 0.68, align_to = "panel")

## ---- Panels 2-4: the estimated parameters -- DO NOT CONVERGE ------------
## Beta and Beta*Eta share a COMMON RELATIVE window: mean * (1 +/- REL_FRAC).
## This is the honest way to make "Beta is loose, Beta*Eta is tight" visible.
## Auto-scaling each panel to its own data range would make both fill the
## panel equally and destroy the comparison; a hand-picked range for each
## would be an axis trick. A shared fractional window is neither: the
## fraction of panel height each trace occupies is then proportional to its
## coefficient of variation, which is exactly the quantity being compared
## (Beta CV 9.6% vs Beta*Eta CV 2.1% -- so Beta should look ~4.6x wider).
REL_FRAC <- 0.32
rel_ylim <- function(var) {
  m <- mean(wide[[var]], na.rm = TRUE)
  c(m * (1 - REL_FRAC), m * (1 + REL_FRAC))
}

beta_ylim <- rel_ylim("Beta")
p_beta <- trace_panel("Beta", expression(bold(beta)), ylim = beta_ylim) |>
  lab(x = 50, y = lab_y(beta_ylim),
      text = sprintf("ends anywhere from %.0f to %.0f", beta_rng[1], beta_rng[2]),
      colour = VERM, hjust = 1, vjust = 0)

rho_ylim <- auto_ylim("Rho", lab_at = "top")
p_rho <- trace_panel("Rho", expression(bold(rho)), ylim = rho_ylim) |>
  lab(x = 50, y = lab_y(rho_ylim, "top"),
      text = sprintf("ends %.2f to %.2f", rho_rng[1], rho_rng[2]),
      colour = VERM, hjust = 1)

eta_ylim <- auto_ylim("Eta", lab_at = "top")
p_eta <- trace_panel("Eta", expression(bold(eta)), ylim = eta_ylim,
                     xlab = "mif2 iteration") |>
  lab(x = 2, y = lab_y(eta_ylim, "top"),
      text = sprintf("ends %.3f to %.3f", eta_rng[1], eta_rng[2]),
      colour = VERM)

## ---- Panel 5: the derived combination -- CONVERGES ----------------------
## SAME relative window as Beta (see REL_FRAC above). The trace should read as
## a near-flat line next to Beta's visibly wandering band -- that contrast is
## the whole figure, and it is earned, not manufactured.
be_ylim <- rel_ylim("BetaEta")
p_be <- trace_panel("BetaEta", expression(bold(beta %*% eta)),
                    ylim = be_ylim, xlab = "mif2 iteration") +
  geom_hline(yintercept = mean(fin$BetaEta), colour = GREEN,
             linetype = "dashed", linewidth = 0.35)
p_be <- lab(p_be, x = 50, y = lab_y(be_ylim),
            text = sprintf("this one DOES converge:\n%.2f to %.2f", be_rng[1], be_rng[2]),
            colour = GREEN, hjust = 1, vjust = 0)

## ---- Panel 6: the caption card ------------------------------------------
## Generous y-range + clip="off" so nothing gets cut at the panel edge.
## ggplot text is sized in POINTS, not data units, so stacking multi-line
## blocks by eye reliably overlaps. Instead: one row per LINE, on a fixed
## pitch. The panel is ~1.5in (~108pt); pitch 0.088 of a 1.05-unit range is
## ~9.0pt against ~6.6pt text, so lines cannot collide by construction.
## Pitch is 0.088 and not 0.095 because 12 rows at 0.095 put the last line at
## y=-0.045, and scale_y_continuous(limits=c(0,1.05)) DROPS out-of-range rows
## silently rather than clipping them. Keep 11*pitch <= 1.0 if you add a line.
pitch <- 0.088
card_txt <- tribble(
  ~text,                                ~col,    ~face,
  sprintf("%d runs · %d iterations", n_rep, n_mif), "grey30", "plain",
  "Same starting parameters;",           "grey30", "plain",
  "only the random seed differs.",       "grey30", "plain",
  ## Load-bearing: gamma is held fixed. If gamma were estimated too, this would
  ## be this lecture's earlier Beta/Gamma R0 ridge, not the Beta*Eta one here.
  sprintf("γ fixed at %g; β, ρ, η estimated.", gamma_fix), "grey30", "plain",
  "",                                    "grey30", "plain",
  "The likelihood agrees.",              GREEN,    "bold",
  "β, ρ, η do not.",                     VERM,     "bold",
  "β×η does.",                           GREEN,    "bold",
  "",                                    "grey30", "plain",
  "Another ridge, like the",             "grey30", "plain",
  "R0 ridge: the data pins a",           "grey30", "plain",
  "combination, not the pieces.",        "grey30", "plain"
) |>
  mutate(y = 1.0 - (row_number() - 1) * pitch)

card <- ggplot(card_txt, aes(x = 0, y = y, label = text)) +
  geom_text(aes(colour = col, fontface = face), hjust = 0, vjust = 1,
            size = 2.6, show.legend = FALSE) +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1.05)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(2, 2, 2, 6))

## ---- Assemble ------------------------------------------------------------
fig <- (p_ll | p_beta | p_rho) / (p_eta | p_be | card)

ggsave("figures/fig-mif2-traces.png", fig,
       width = 5.5, height = 3.4, dpi = 300, bg = "white")

cat("wrote figures/fig-mif2-traces.png\n")
