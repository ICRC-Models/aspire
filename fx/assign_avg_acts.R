#######################################################################################
# 
# Kathryn Peebles
# 2017-11-19
# assign_avg_acts: function to assign mean number of monthly coital acts to each woman.
#
# input:  Subset of f_dt containing only "country" variable; data tables of country-specific weekly coital frequency and condom use.
# output: Vector of values indicating individual-specific average monthly coital frequency
#
#######################################################################################

assign_avg_acts <- function(countries) {
  
  load(paste0(getwd(), "/data/acts_dt.RData"))
  
  n_acts <- unname(sapply(countries, function(x) {
    rgamma(n = 1, shape = acts_dt[country == x, alpha], rate = acts_dt[country == x, beta])
  }))
  
  n_acts <- n_acts * 30/7 # Mean acts reported in acts_dt is for previous 7 days. Adjust to 30-day increments.
  
  # ## Check that mean number of acts by country approximates proportion reported in acts_dt
  # sapply(unique(countries), function(x) { mean(n_acts[f_dt[, which(country == x)]]) })

  return(n_acts)
}