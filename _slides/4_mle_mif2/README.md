# Lecture B — Finding the MLE with iterated filtering

*SISMID: Simulation-based Inference for Epidemiological Dynamics.*
Slot: July 24, 11:00–12:30pm. Instructor: Qianying (Ruby) Lin.

**Finding the parameters that *maximise* the likelihood (the MLE).**
Built from the former `5_IF_theory`, plus the likelihood-surface material
from `4_pfilter_in_pomp`. Follows directly from `A_likelihood_pfilter`.

Lecture A evaluated `L(θ)` at a fixed θ. This lecture finds
`argmax_θ L(θ)`: why a standard optimizer can't (the noisy, ridged surface)
and the algorithm that can (iterated filtering, `mif2`).

## For students

1. Open `4_mle_mif2.Rproj` **first**. Every path assumes the project root is
   the working directory.
2. Work through `scripts/EXERCISES.R` (`mif2` basics).
3. Answers are in `scripts/SOLUTIONS.R`.

One script works on macOS, Linux and Windows via `makePSOCKcluster()`.

## Layout

| Path | What it is |
|:--|:--|
| `main.qmd` → `main.pdf` | the deck |
| `data/` | `Measles_Consett_1948.csv` |
| `scripts/model_measSIR.R` | **the** model. Byte-identical to `A_likelihood_pfilter`'s copy. |
| `scripts/model_measSEIR.R` | the SIR + a latent `E` compartment (`Sigma=1`, a 1-week latent period). Used only by the lunch-break exercise. |
| `scripts/EXERCISES.R` | the student worksheet (`mif2` basics; Ex 6 fits the SEIR) |
| `scripts/SOLUTIONS.R` | worked answers, with the real numbers |
| `scripts/slide-code.R` | the key code shown on the slides, as one runnable script |
| `figures/fig-*.R` | one script per figure |
| `results/local_search.rds` | 20 replicated `mif2` runs, Np=2000, Nmif=50 |
| `results/lik_local.rds` | honest `pfilter` re-estimates at those endpoints |
| `results/pfilter-grid1.rds` | 40×40 β–γ likelihood grid (the surface / ridge figures) |
| `results/nelder-mead-fails.rds` | the failed-optimizer demo |

The surface figures — `fig-likelihood-ridge.{R,png}` and
`fig-nelder-mead-fails.{R,png}` — moved here from Lecture A: they are about
*finding* the maximum, so they belong with the MLE. Both read
`results/pfilter-grid1.rds` (a copy also lives in Lecture A, which still
needs it for its own figures).

## What this lecture covers

The MLE (`argmax L`) → the likelihood **surface is a knife-edge ridge**
(constant R₀) → a standard optimizer **fails** on that noisy ridge →
**iterated filtering** (`mif2`): the four knobs, cooling, reading the
convergence traces, and why the number `mif2` prints is not your likelihood.

The slides stay on the SIR throughout. **Exercise 6 is a lunch-break run**:
students fit the SEIR (the SIR + latent `E` compartment they built in
Lecture 2) with the same `mif2` machinery and compare it to the SIR fit. The
verified result: SEIR reaches ≈ −104.3 vs the SIR's ≈ −105.4 — about 1.1 log
units for one extra parameter (`Sigma`), an AIC tie — and `Sigma` is barely
identified (it ranges ~1.2–3.4 across searches). The latent period does not
clearly earn its keep for a single outbreak.

## The one trap that will bite you

Build the parameter transformation **once, outside** any `foreach` loop:

```r
measSIR |> pomp(
  partrans=parameter_trans(log="Beta", logit=c("Rho","Eta")),
  paramnames=c("Beta","Rho","Eta")
) -> po
```

Inside a parallel loop it forces a recompile on every worker. This is the
issue flagged in the 2026 planning document. Also: `foreach` finds
`measSIR`/`po` on its own — passing `.export=` prints a warning.

## Caches must come from *our* model

`results/local_search.rds` **must** be generated from
`scripts/model_measSIR.R` (`Gamma=0.5`).

> The legacy `mif/results/local_search.rds` was built with `Gamma=2` — a
> genuinely different model. Dropping it in produces a figure that silently
> contradicts the exercises (loglik ≈ −115 with β ending 26–36, versus
> ≈ −106 with β ending 16–23). `figures/fig-mif2-traces.R` reads `Gamma` off
> the cache and re-labels itself, but the numbers will still disagree.

## Conventions

Same as `A_likelihood_pfilter` — see `../A_likelihood_pfilter/README.md`.
In short: sketch figures carry the idea, generated figures carry the
evidence; amber is always the data; never hardcode axis limits; never write
a number you did not compute.
