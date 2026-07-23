# ==============================================================================
# IN-CLASS EXERCISES: The Zombie Apocalypse (Independent Modeling)
# ==============================================================================
# In this exercise, you will translate your conceptual diagram of the zombie
# outbreak into a working pomp object. You will need to rely on your previous
# experience building the SEIR and SIR models!
# ==============================================================================

library(pomp)
library(tidyverse)

# 1. Load the simulated dataset
# zombie_df <- read.csv('raw-data/zombieCases.csv')
# head(zombie_df)

# ==============================================================================
# EXERCISE 1: Conceptualization & Variables (No Code Required)
# ==============================================================================
# Discuss with your group and write down:
# 1. What are your State Variables? (e.g., S, E, Z, R)
# 2. What is your Incidence Tracker? (You need this to match the weekly 'cases' data)
# 3. What are your Parameters? (Transmission rates, reporting fraction, etc.)
#
# (Keep this list handy; you will need it for the statenames and paramnames!)


# ==============================================================================
# EXERCISE 2: Writing the C-Snippets
# ==============================================================================

# ------------------------------------------------------------------------------
# A. The Process Model (rprocess)
# ------------------------------------------------------------------------------
# BUILD YOUR TRANSITIONS HERE.
#
# HINT 1: The Exposed class splits into two pathways (Zombies or Immune).
#         You should use the reulermultinom() function for this specific transition.
# HINT 2: (Optional) If you want the mutation to "turn on" suddenly at t=0,
#         C-code uses the following syntax for conditional logic:
#         double my_rate = (t < 0) ? 0.0 : normal_rate;
# HINT 3: Don't forget to add new zombies to your incidence tracking variable!

zombie_rprocess <- Csnippet("
  // --- YOUR CODE HERE ---



")

# ------------------------------------------------------------------------------
# B. The Measurement Model (rmeasure & dmeasure)
# ------------------------------------------------------------------------------
# BUILD YOUR MEASUREMENT CODE HERE.
#
# HINT: The scenario states that reporting infrastructure decays as the total
# Zombie horde grows. You will need to calculate a dynamic reporting rate
# (e.g., using an exponential decay function) before evaluating the likelihood.

zombie_dmeas <- Csnippet("
  // --- YOUR CODE HERE ---


")

zombie_rmeas <- Csnippet("
  // --- YOUR CODE HERE ---


")

# ------------------------------------------------------------------------------
# C. The Initial State (rinit)
# ------------------------------------------------------------------------------
# Initialize your compartments here.

zombie_rinit <- Csnippet("
  // --- YOUR CODE HERE ---


")


# ==============================================================================
# EXERCISE 3: Compiling the pomp Object
# ==============================================================================
# Bring it all together. Fill in the missing arguments based on your snippets.

# zombie_pomp <- zombie_df |>
#   pomp(
#     times = "time",
#     t0 = -60,                  # Hint: Patient zero infected 60 days before data starts
#     rprocess = euler(zombie_rprocess, delta.t = 1/7),
#     rinit = zombie_rinit,
#     rmeasure = zombie_rmeas,
#     dmeasure = zombie_dmeas,
#     accumvars = c(__________), # DON'T FORGET THIS!
#     statenames = c(__________),
#     paramnames = c(__________)
#   )


# ==============================================================================
# EXERCISE 4: Simulating and Testing
# ==============================================================================
# Before trying to fit the model with mif2, pick a reasonable set of starting
# parameters and run a simulation. Does the model compile? Do the simulated
# trajectories look somewhat like an epidemic curve?

# test_params <- c(
#   N = 100000,
#   # --- ADD YOUR OTHER PARAMETERS AND GUESSES HERE ---
# )
#
# sims <- simulate(zombie_pomp, params = test_params, nsim = 20,
#                  format = "data.frame", include.data = TRUE)
#
# ggplot(sims, aes(x = time, y = cases, group = .id, color = (.id == "data"))) +
#   geom_line(alpha = 0.5) +
#   theme_bw() +
#   labs(title = "Initial Simulation Check")
