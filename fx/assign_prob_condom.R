#######################################################################################
# 
# Kathryn Peebles
# 2017-11-19
# assign_prob_condom: function to assign individual probability of condom use. Called in creation of data table to store individual participant characteristics.
#
# input:  Subset of f_dt containing only "country" variable; data tables of country-specific coital frequency and condom use.
# output: Vector of values bounded in [0, 1] indicating individual-specific per-act probability of condom use
#
#######################################################################################

assign_prob_condom <- function(countries) {
  
  load(paste0(getwd(), "/data/cond_dt.RData"))
  load(paste0(getwd(), "/data/acts_dt.RData"))
  
  lambda <- cond_dt[, cond_prob] * acts_dt[, mean * 30/7] # Mean number of condom-protected acts. Acts reported for previous 7 days; convert to 30-day increments
  names(lambda) <- cond_dt$country
  prob_condom <- unname(sapply(countries, function(x) {
    rpois(n = 1, lambda = lambda[x])/acts_dt[country == x, mean * 30/7]
  }))
  prob_condom[prob_condom > 1] <- 1
  
  # ## Check that mean of condom probabilities by country approximates proportion reported in cond_dt
  # sapply(unique(countries), function(x) { mean(prob_condom[f_dt[, which(country == x)]]) })
  
  return(prob_condom)
}