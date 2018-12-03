#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   19 March 2018
# assign_prob_condom: function to assign individual probability of condom use.
#
# input:  Subset of f_dt containing only "country" variable; data table of country-specific condom probability.
# output: Vector of values bounded in [0, 1] indicating individual-specific per-act probability of condom use
#
#######################################################################################

assign_prob_condom <- function(countries) {
  
  load(paste0(getwd(), "/data-public/cond_dt.RData"))
  
  prob_condom <- unlist(unname(sapply(unique(countries), function(x) {
    rbeta(n = sum(countries == x), shape1 = cond_dt[country == x, alpha], shape2 = cond_dt[country == x, beta])
  })))
  
  # ## Check that mean of condom probabilities by country approximates proportion reported in cond_dt
  # sapply(unique(countries), function(x) { mean(prob_condom[f_dt[, which(country == x)]]) })
  # sapply(unique(countries), function(x) { min(prob_condom[f_dt[, which(country == x)]]) })
  # sapply(unique(countries), function(x) { max(prob_condom[f_dt[, which(country == x)]]) })
  
  return(prob_condom)
}