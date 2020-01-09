#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   31 March 2018
# abc_samples_process: Function to process sample of particles in iteration t and create sample of particles to run in iteration t + 1
#
# input:  Values for t, N, alpha, the stopping criterion, and the number of nodes available to run the next sample of particles
# output: Either none or the set of particles and accompanying scripts to simulate the next sample of particles.
#
#######################################################################################

abc_samples_process <- function(t, N, alpha, p_acc_min, n_nodes) {
  library(mnormt)
  library(data.table)
  setwd("~/Documents/code/aspire/abc/hyak")
  source("~/Documents/code/aspire/abc/abc-fx.R")
  load("~/Documents/code/aspire/abc/priors.RDATA", envir = .GlobalEnv)
  
  t_prev <- t - 1
  t_curr <- t
  t_next <- t + 1
  
  ## Create empty data table in which to store particles
  n_particles <- length(list.files(path = paste0(getwd(), "/t", t_curr, "/particles")))
  particles <- data.table(sapply(c("lambda", "cond_rr", "c", "s", "rr_ai", "p_rate_rr", "base_male_hiv_incidence", "rho", "i"), function(x) { x = rep(NA_real_, n_particles) }))
  
  ## Load all particles into empty data table
  for(i in 1:n_particles) {
    load(paste0(getwd(), "/t", t_curr, "/particles/particle_", i, ".RDATA"))
    particles[i, names(particles) := particle]
  }
  
  ## Calculate weights
  if(t_curr == 0) {
    ## Set w to 1 for all particles in iteration 0
    particles[, w := 1]
  } else if(t_curr > 0) {
    ## Load particles from iteration t-1
    load(paste0(getwd(), "/t", t_prev, "/particles_acc_", t_prev, ".RDATA"))
    prev_particles <- eval(parse(text = paste0("particles_acc_", t_prev)))
    
    particles[, w := calc_weights(prev_iter_particles = prev_particles[, .(lambda, cond_rr, c, s, rr_ai, p_rate_rr, base_male_hiv_incidence)], 
                                  weights             = prev_particles[, w], 
                                  curr_iter_particles = particles[, .(lambda, cond_rr, c, s, rr_ai, p_rate_rr, base_male_hiv_incidence)])]
  }
  
  ## Specify iteration
  particles[, t := t_curr]
  
  if(t > 0) {
    ## Append particles in iteration t-1 to particles in iteration t.
    particles <- rbind(prev_particles, particles)
  }
  
  ## Rho is likelihood value of parameters given data. We want to maximize this value, so order in descending order, and keep top alpha-proportion
  setorderv(x = particles, cols = "rho", order = -1)
  
  particles <- particles[1:(nrow(particles) * alpha)]
  
  ## Normalize weights
  particles[, w := w/sum(w)]
  
  ## Save accepted particles from iteration t in appropriate directory
  particles_set_name <- paste0("particles_acc_", t_curr)
  assign(x = particles_set_name, value = particles)
  file_name <- paste0("particles_acc_", t, ".RDATA")
  save(list = particles_set_name, file = paste0(getwd(), "/t", t_curr, "/", file_name))
  
  ## Calculate the proportion of accepted particles and compare to the minimum acceptance criterion
  if(t == 0) {
    p_acc <- 1
  } else {
    p_acc <- nrow(particles[t == t_curr])/nrow(particles)
  }
  
  ## If the proportion of accepted particles is greater than the minimum acceptance criterion, increase value of t and select new batch of particles to simulate
  if(p_acc > p_acc_min) {
    prev_particles <- particles
    
    new_particles <- select_particles(dt = prev_particles, n_new_particles = N - N * alpha)
    
    new_particles <- perturb_particles(dt = new_particles, cov_mat = (0.5 * cov.wt(x = as.matrix(prev_particles[, .(lambda, cond_rr, c, s, rr_ai, p_rate_rr, base_male_hiv_incidence)]), wt = prev_particles[, w])$cov))
    
    new_particles <- new_particles[, .(lambda_new, cond_rr_new, c_new, s_new, rr_ai_new, p_rate_rr_new, base_male_hiv_incidence_new)]
    new_particles[, i := 1:.N]
    setnames(new_particles, old = c("lambda_new", "cond_rr_new", "c_new", "s_new", "rr_ai_new", "p_rate_rr_new", "base_male_hiv_incidence_new"), new = c("lambda", "cond_rr", "c", "s", "rr_ai", "p_rate_rr", "base_male_hiv_incidence"))
    
    ## Convert data table of new particles to simulate to list
    particles <- lapply(1:nrow(new_particles), function(x) { new_particles[x, ] })
    
    if(!dir.exists(paths = paste0(getwd(), "/t", t_next))) {
      dir.create(path = paste0(getwd(), "/t", t_next))
    }
    
    save(particles, file = paste0(getwd(), "/t", t_next, "/particles.RDATA"))
    
    ## Create vectors for start and end particles for each script
    n_new_particles <- length(particles)
    batch_size <- floor(n_new_particles/n_nodes)
    starts <- seq(1, n_new_particles - batch_size - 1, batch_size)
    ends   <- seq(batch_size, n_new_particles, batch_size)
    ends[length(ends)] <- n_new_particles
    
    library(readr)
    
    ## Read in template script
    template <- read_lines(file = paste0(getwd(), "/template-script.R"))
    
    if(!dir.exists(paths = paste0(getwd(), "/t", t_next, "/sim-scripts"))) {
      dir.create(path = paste0(getwd(), "/t", t_next, "/sim-scripts"))
    }
    
    for(i in 1:length(starts)) {
      particles_to_sim <- c(paste0("i_start = ", starts[i]),
                            paste0("i_end = ", ends[i]))
      
      sim_script <- c(particles_to_sim, template)
      
      filename <- paste0(getwd(), "/t", t_next, "/sim-scripts/sim_", i, ".R")
      
      write_lines(x = sim_script, path = filename)
    }
    
    # Create bash script
    top_line <- "#!/bin/bash"
    
    jobs_loop <- paste0("for i in {1..", n_nodes, "}; do\n\texport SIMNO=$i\n\tsbatch runsim_abc_aspire.sh\ndone")
    
    master <- c(top_line, jobs_loop)
    
    write_lines(master, paste0(getwd(), "/t", t_next, "/master_t", t_next, ".sh"))
    
    ## Transfer sim-scripts folder, master_t#.sh, and particles.RDATA to Hyak. bash master_t#.sh
  } else {
    print("Proportion accepted particles less than minimum acceptance criterion. Approximate convergence achieved - yay!")
  }
}
