## Purpose: Create data table of priors for use in calculating numerator of ABC iterations and in sampling particles in iteration 0.
## Author:  Kathryn Peebles
## Date:    22 March 2018

library(data.table)

param_names <- c("prop_ai", "a", "b")

priors <- data.table(param = param_names, 
                     dist = rep("unif", length(param_names)),
                     min = c(0.00, 0.01, 0.01),
                     max = c(1.00, 7.00, 7.00))

save(priors, file = "~/Documents/code/aspire/abc/priors.RDATA")
     