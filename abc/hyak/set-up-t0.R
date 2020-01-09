## Purpose: Create first set of particles and create individual scripts to run simulations in parallel on Hyak
## Author:  Kathryn Peebles
## Date:    22 March 2018

setwd("~/Documents/code/aspire")

load(paste0(getwd(), "/abc/priors.RDATA"))

options(scipen = 999)

N <- 100000

# Sample particles from priors
particles <- data.table(lambda                  = rbeta(n      = N,
                                                        shape1 = priors[param == "lambda", alpha],
                                                        shape2 = priors[param == "lambda", beta]), 
                        cond_rr                 = runif(n   = N,
                                                        min = priors[param == "cond_rr", min],
                                                        max = priors[param == "cond_rr", max]),
                        c                       = runif(n   = N,
                                                        min = priors[param == "c", min],
                                                        max = priors[param == "c", max]),
                        s                       = runif(n   = N,
                                                        min = priors[param == "s", min],
                                                        max = priors[param == "s", max]),
                        rr_ai                   = runif(n   = N,
                                                        min = priors[param == "rr_ai", min],
                                                        max = priors[param == "rr_ai", max]),
                        p_rate_rr               = runif(n   = N,
                                                        min = priors[param == "p_rate_rr", min],
                                                        max = priors[param == "p_rate_rr", max]),
                        base_male_hiv_incidence = runif(n   = N,
                                                        min = priors[param == "base_male_hiv_incidence", min],
                                                        max = priors[param == "base_male_hiv_incidence", max]),
                        i                       = 1:N)

# Convert data table to list
particles <- lapply(1:nrow(particles), function(x) { particles[x, ] })

if(!dir.exists(paths = paste0(getwd(), "/abc/hyak/t0"))) {
  dir.create(path = paste0(getwd(), "/abc/hyak/t0"))
}

save(particles, file = paste0(getwd(), "/abc/hyak/t0/particles.RDATA"))

# Create vectors for start and end particles for each script
n_new_particles <- N
batch_size <- floor(N/14)

starts <- seq(1, n_new_particles - batch_size - 1, batch_size)
ends   <- seq(batch_size, n_new_particles, batch_size)
ends[length(ends)] <- n_new_particles

library(readr)

## Read in template script
template <- read_lines(file = paste0(getwd(), "/abc/hyak/template-script.R"))

if(!dir.exists(paths = paste0(getwd(), "/abc/hyak/t0/sim-scripts"))) {
  dir.create(path = paste0(getwd(), "/abc/hyak/t0/sim-scripts"))
}

for(i in 1:length(starts)) {
  particles_to_sim <- c(paste0("i_start = ", starts[i]),
                        paste0("i_end = ", ends[i]))
  
  sim_script <- c(particles_to_sim, template)
  
  filename <- paste0(getwd(), "/abc/hyak/t0/sim-scripts/sim_", i, ".R")
  
  write_lines(x = sim_script, path = filename)
}


# Create bash script
top_line <- "#!/bin/bash"

n_nodes <- 14

jobs_loop <- paste0("for i in {1..", n_nodes, "}; do\n\texport SIMNO=$i\n\tsbatch runsim_abc_aspire.sh\ndone")

master <- c(top_line, jobs_loop)

write_lines(master, paste0(getwd(), "/abc/hyak/t0/master_t0.sh"))
