## Purpose: Create data table of priors for use in calculating numerator of ABC iterations and in sampling particles in iteration 0.
## Author:  Kathryn Peebles
## Date:    27 January 2020

library(data.table)

param_names <- c("lambda", "cond_rr", "c", "s", "rr_ai", "p_rate_rr", "base_male_hiv_incidence")

priors <- data.table(param = param_names, 
                     dist  = c("beta", rep("unif", length(param_names) - 1)),
                     min   = c(NA, 1, 0, 0,  5,  1, 0.0097),
                     max   = c(NA, 3, 1, 1, 16, 10, 0.0988),
                     alpha = c(5, rep(NA, length(param_names) - 1)),
                     beta  = c(1664, rep(NA, length(param_names) - 1)))

save(priors, file = "~/Documents/code/aspire/abc/priors.RDATA")
     