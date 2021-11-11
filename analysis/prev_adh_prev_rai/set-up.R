## Purpose: Run test simulations of analysis parameters
## Author:  Kathryn Peebles
## Date:    16 October 2019

setwd("~/Documents/code/aspire")

# Attach packages
library(data.table)
library(readr)

# Specify the combinations of prevalence of adherence and RAI to simulate
params <- as.data.table(expand.grid(prop_adh = seq(0.5, 1, 0.05), prev_ai = seq(0, 30, 5), rr_ring = c(0.20, 0.25, 0.30, 0.35, 0.40)))

# Read in template script
template <- read_lines(file = paste0(getwd(), "/analysis/prev_adh_prev_rai/template-script.R"))

sim_dir <- "/analysis/prev_adh_prev_rai/sim-scripts-prev-adh-ai"

if(!dir.exists(paths = paste0(getwd(), sim_dir))) {
  dir.create(path = paste0(getwd(), sim_dir))
}

for(i in 1:nrow(params)) {
  param_sets_to_sim <- c(paste0("sim_prev_ai  <- \"", params[i, prev_ai], "\""),
                         paste0("sim_prop_adh <- ", params[i, prop_adh]),
                         paste0("sim_rr_ring  <- ", params[i, rr_ring]))
  
  sim_script <- c(param_sets_to_sim, template)
  
  filename <- paste0(getwd(), sim_dir, "/sim_", i, ".R")
  
  write_lines(x = sim_script, path = filename)
}

# Create bash script
top_line <- "#!/bin/bash"

n_jobs <- nrow(params)

jobs_loop <- paste0("for i in {1..", n_jobs, "}; do\n\texport SIMNO=$i\n\tsbatch runsim_aspire_prev_adh_ai.sh\ndone")

master <- c(top_line, jobs_loop)

write_lines(master, paste0(getwd(), "/analysis/prev_adh_prev_rai/master_aspire_prev_adh_ai.sh"))
