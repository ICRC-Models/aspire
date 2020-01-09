
setwd("/gscratch/csde/kpeebles")

# Load list of particles
load(file = paste0(getwd(), "/particles.RDATA"))

# Source functions
source(paste0(getwd(), "/set-up.R"))

# Run set of simulations
lapply(X   = particles[i_start:i_end],
       FUN = function(x) { run_sim(lambda    = as.numeric(x[, "lambda"]),
                                   cond_rr   = as.numeric(x[, "cond_rr"]),
                                   c         = as.numeric(x[, "c"]),
                                   s         = as.numeric(x[, "s"]),
                                   rr_ai     = as.numeric(x[, "rr_ai"]),
                                   i         = as.numeric(x[, "i"]),
                                   p_rate_rr = as.numeric(x[, "p_rate_rr"]),
                                   base_male_hiv_incidence = as.numeric(x[, "base_male_hiv_incidence"]), 
                                   prev_ai   = "prop_ai_obs",
                                   prop_adh  = "prop_adh_obs",
                                   rr_ring   = 0.25,
                                   calibrate = T) })