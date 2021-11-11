## Purpose: Run test simulations of analysis parameters
## Author:  Kathryn Peebles
## Date:    16 October 2019

setwd("~/Documents/code/aspire")

# Attach packages
library(data.table)
library(readr)

# Specify the combinations of prevalence of adherence and RAI to simulate
params <- as.data.table(expand.grid(prop_ai = sort(c(seq(0.05, 0.30, 0.05), 0.063)), rr_ring = c(0.25, 0.30, 0.35), re_randomize = c(T, F)))
addl_params <- as.data.table(expand.grid(prop_ai = sort(c(0.01, 0.027, 0.167)), rr_ring = c(0.25, 0.30, 0.35), re_randomize = T))

params <- rbindlist(l = list(params, addl_params))

# Read in template script
template <- read_lines(file = paste0(getwd(), "/analysis/per_exp_rr_ai/template-script.R"))

sim_dir <- "/analysis/per_exp_rr_ai/sim-scripts-ai"

if(!dir.exists(paths = paste0(getwd(), sim_dir))) {
  dir.create(path = paste0(getwd(), sim_dir))
}

for(i in 1:nrow(params)) {
  param_sets_to_sim <- c(paste0("sim_prop_ai      <- ", params[i, prop_ai]),
                         paste0("sim_rr_ring      <- ", params[i, rr_ring]),
                         paste0("sim_re_randomize <- ", params[i, re_randomize]))
  
  sim_script <- c(param_sets_to_sim, template)
  
  filename <- paste0(getwd(), sim_dir, "/sim_", i, ".R")
  
  write_lines(x = sim_script, path = filename)
}

# Create bash script
top_line <- "#!/bin/bash"

n_jobs <- nrow(params)

jobs_loop <- paste0("for i in {1..", n_jobs, "}; do\n\texport SIMNO=$i\n\tsbatch runsim_aspire_ai_analysis.sh\ndone")

master <- c(top_line, jobs_loop)

write_lines(master, paste0(getwd(), "/analysis/per_exp_rr_ai/master_aspire_ai_analysis.sh"))
