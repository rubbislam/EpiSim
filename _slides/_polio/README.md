## Polio lesson

The polio lesson exemplifies a more sophisticated workflow for a real scientific problem.

----------------------------

### Run-levels

The workflow illustrated in this lesson makes use of "run-levels", as controlled by an **R** variable `run_level`.
The choice of run-level affects various algorithmic parameters in such a way that higher values of the run-level result in more extensive (and expensive), but more thorough calculations.
At run-level 1, the codes are quickly completed: this is useful for debugging.
At run-level 2, more computational resources are mobilized: this is useful for checking that the computations will scale.
Run-level 3 is intended to produce meaningful results.

In the codes herein, one changes the run-level by setting the environment variable `RUNLEVEL` to 1, 2, or 3.
Alternatively, one can edit the codes to set the run-level explicitly.

The slides and notes are generated using run-level 3: the directory `results` contains the archived computational results.
The results of running the codes in `main.R` at run-levels 1 and 2 are contained in the directories `runlevel1`, `runlevel2`, respectively.


----------------------------

### Manifest

- index.md: main lesson page
- main.Rnw, main.R: content for notes and slides
- notes.Rnw, slides.Rnw: shells for notes and slides
- CLUSTER.R: local cluster setup
- polio_wisconsin.csv: data
- params.csv: database of parameter estimates
- polio_fig1A.png: model diagram
- results/: archive files
- run3.sbat: Slurm script
- algorithmic-parameters-exercise.Rmd: worked solution
- convergence-exercise.Rmd: worked solution
- demography-exercise.Rmd: worked solution
- subdirectories:
  - results: archive files for `main.Rnw` (RUNLEVEL=3)
  - runlevel1: results of running `main.R` with RUNLEVEL=1
  - runlevel2: results of running `main.R` with RUNLEVEL=2
  - initial-values-exercise: worked solution
  - starting-values-exercise: worked solution
