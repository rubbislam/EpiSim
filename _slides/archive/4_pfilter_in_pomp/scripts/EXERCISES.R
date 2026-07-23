library(tidyverse)
library(pomp)
library(iterators)
library(doFuture)


# Setup -------------------------------------------------------------------

# Building the starting model
source("scripts/model_measSIR.R")

# Exercise 1 --------------------------------------------------------------

# -   How much computer processing time does a particle filter take?
# -   How does this scale with the number of particles?

## HINT: system.time( <code> )[3] can be used to time something. 
## HINT: You don't have to, but you can create a foreach loop, and time 
## how long it takes. Something like:
## 
## my_vec <- c(1, 1, 1)
## foreach(i=my_vec, .combine = c) %dopar% {system.time(<code>)[3]} -> results


# Exercise 2 --------------------------------------------------------------

# Compute several 1-dim likelihood slices in the $\eta$ direction.
# 
# Hint: Use the slice_design function like in the slides. If you need help, 
#       try looking at the help page: ?slice_design.
# Hint: Try using parallel computing to speed this up. 

# Exercise 3 --------------------------------------------------------------

# Compute a two-dimensional slice of the likelihood in the $\beta$-$\eta$ plane.
# 
# Hint: Use the expand.grid function like in the slides. If you need help, 
#       try looking at the help page: ?expand.grid. 
# Hint: Try using parallel computing to speed this up. 

