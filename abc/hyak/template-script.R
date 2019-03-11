
setwd("/gscratch/csde/kpeebles")

## Load list of particles
load(file = paste0(getwd(), "/particles.RDATA"))

## Source functions
source(paste0(getwd(), "/set-up.R"))
source(paste0(getwd(), "/fx/run_sim_abc.R"))

library(parallel)

## Run set of simulations
lapply(X   = particles[i_start:i_end],
       FUN = function(x) { run_sim_abc(lambda  = as.numeric(x[, "lambda"]),
                                       cond_rr = as.numeric(x[, "cond_rr"]),
                                       c       = as.numeric(x[, "c"]),
                                       s       = as.numeric(x[, "s"]),
                                       rr_ai   = as.numeric(x[, "rr_ai"]),
                                       i       = as.numeric(x[, "i"])) })