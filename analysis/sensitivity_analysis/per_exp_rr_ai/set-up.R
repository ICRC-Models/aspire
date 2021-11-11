## Purpose: Create scripts to run sensitivity analyses of the relative risk of RAI
## Author:  Kathryn Peebles
## Date:    12 November 2019

library(data.table)

# Load particles from final iteration of calibration process
load("~/Documents/code/aspire/abc/hyak/t13/particles_acc_13.RDATA")

particles_acc_13[, wt_rho := rho * w]
setorderv(particles_acc_13, cols = "wt_rho", order = -1)

# Identify the top-fitting set of 10 particles that represent the range of values in the prior distribution for rr_ai.
sensitivity_particles <- particles_acc_13[c(1, 2, 3, 5, 7, 9, 10, 12, 17, 51)]

setorderv(sensitivity_particles, cols = "rr_ai", order = 1)
sensitivity_particles[, rr_ai_round := round(rr_ai)]

setwd("~/Documents/code/aspire")

# Attach packages
library(data.table)
library(readr)

# Specify the combinations of proportion of acts that are RAI and ring efficacy to simulate
params <- as.data.table(expand.grid(prop_ai = sort(c(seq(0.05, 0.30, 0.05), 0.063)), rr_ring = 0.30, re_randomize = T))

params <- params[rep(1:.N, 10)]
sensitivity_particles <- sensitivity_particles[rep(1:.N, 7)]

setkey(params, prop_ai)

params <- as.data.table(cbind(params, sensitivity_particles))

# Read in template script
template <- read_lines(file = paste0(getwd(), "/analysis/sensitivity_analysis/per_exp_rr_ai/template-script.R"))

if(!dir.exists(paths = paste0(getwd(), "/analysis/sensitivity_analysis/per_exp_rr_ai/sim-scripts-ai-sens"))) {
  dir.create(path = paste0(getwd(), "/analysis/sensitivity_analysis/per_exp_rr_ai/sim-scripts-ai-sens"))
}

for(i in 1:nrow(params)) {
  
  param_sets_to_sim <- c(paste0("sim_prop_ai      <- ", params[i, prop_ai]),
                         paste0("sim_rr_ring      <- ", params[i, rr_ring]),
                         paste0("sim_re_randomize <- ", params[i, re_randomize]),
                         paste0("val_base_male_hiv_incidence <- ", params[i, base_male_hiv_incidence]),
                         paste0("val_c                       <- ", params[i, c]),
                         paste0("val_cond_rr                 <- ", params[i, cond_rr]),
                         paste0("val_lambda                  <- ", params[i, lambda]),
                         paste0("val_p_rate_rr               <- ", params[i, p_rate_rr]),
                         paste0("val_rr_ai                   <- ", params[i, rr_ai]),
                         paste0("val_s                       <- ", params[i, s]),
                         paste0("rr_ai_round                 <- ", params[i, rr_ai_round]))
  
  sim_script <- c(param_sets_to_sim, template)
  
  filename <- paste0(getwd(), "/analysis/sensitivity_analysis/per_exp_rr_ai/sim-scripts-ai-sens/sim_", i, ".R")
  
  write_lines(x = sim_script, path = filename)
}

# Create bash script
top_line <- "#!/bin/bash"

n_jobs <- nrow(params)

jobs_loop <- paste0("for i in {1..", n_jobs, "}; do\n\texport SIMNO=$i\n\tsbatch runsim_aspire_ai_sens_analysis.sh\ndone")

master <- c(top_line, jobs_loop)

write_lines(master, paste0(getwd(), "/analysis/sensitivity_analysis/per_exp_rr_ai/master_aspire_ai_sens_analysis.sh"))
