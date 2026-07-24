# Lecture C — Genealogies and phylodynamics with phylopomp

*SISMID: Simulation-based Inference for Epidemiological Dynamics.*
Instructor: Qianying (Ruby) Lin. Built from the phylopomp material in the
former `9_extensions`, updated to the current package.

**The data are now a genealogy (a pathogen family tree), not a case-count
series — but the inference is the same POMP pipeline.** Lecture A evaluated
`L(θ)` and Lecture B maximised it, both for a case-count model. This lecture
shows that once a genealogy is turned into a `pomp` object, `pfilter`,
`logmeanexp`, and `mif2` all transfer unchanged.

## Requires the up-to-date phylopomp

```r
# install.packages("pak")
pak::pak("kingaa/phylopomp")   # 0.19.4.2 or newer
```

**The API changed from older versions.** Simulate with the per-model
`run...()` functions (`runSIR`, `runSIIR`, ...); the old
`simulate("SIR", ...)` string form was removed. Initial states (`S0`, `I0`,
...) are now **fractions** of `pop`, not counts.

## For students

1. Open `10_phylopomp.Rproj` **first**. Every path assumes the project root
   is the working directory.
2. Work through `scripts/EXERCISES.R`.
3. Answers, with the real numbers, are in `scripts/SOLUTIONS.R`.

## Layout

| Path | What it is |
|:--|:--|
| `main.qmd` → `main.pdf` | the deck |
| `scripts/EXERCISES.R` | the student worksheet (simulate, two strains, infer) |
| `scripts/SOLUTIONS.R` | worked answers, with verified numbers |
| `scripts/slide-code.R` | the key code shown on the slides, as one runnable script |
| `figures/fig-*.R` | one script per figure; each writes its `fig-*.png` |
| `results/loglik-vs-beta.rds` | cached β-sweep for the inference figure |

There is **no `data/` folder and no `model_*.R`**: the "data" is a genealogy
that `phylopomp` simulates, and the models are built in to the package
(`runSIR`, `sir_pomp`, ...).

## Figures

All four are generated from the package (run from the project root):

```sh
Rscript figures/fig-sir-genealogy.R   # fig-sir-genealogy.png + fig-sir-ltt.png
Rscript figures/fig-siir-strains.R    # fig-siir-strains.png (two strains, coloured)
Rscript figures/fig-loglik-vs-beta.R  # fig-loglik-vs-beta.png (+ caches the sweep)
```

Each script sets a seed, so the trees and numbers reproduce exactly. The
concept figures on the motivation slides (`nextstrain_geneal.png`,
`phylogen_phylodyn.pdf`, `MGP_events.pdf`) are shared, in `../graphics/`.

## What this lecture covers

A genealogy is a **partially observed** epidemic → **simulate** trees from
compartmental models (`runSIR`, and `runSIIR` for two strains) → read the
tree and its **lineages-through-time** → **evaluate the likelihood** of a
model given the tree with the *same* particle filter as Lectures A and B →
the likelihood **peaks at the truth**, and `mif2` would maximise it.

## The one thing that will trip you up

The simulator interface is the modern one: **per-model `run...()`
functions**, and **fractional** initial states. `runSIR(..., S0 = 0.97,
I0 = 0.03)` means 97%/3% of `pop`, not 0.97 and 0.03 individuals.

## Conventions

Same as `A_likelihood_pfilter` and `B_mle_mif2` — see
`../A_likelihood_pfilter/README.md`. In short: generated figures carry the
evidence and recompute their own annotations; never write a number you did
not compute. The inference tools in `phylopomp` are newer than `pfilter`
itself, so keep trees modest and particle counts generous while exploring.
