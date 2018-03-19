#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   15 March 2018
# assign_male_hiv_status: Function to assign male partner HIV status as a binomial draw with probability equal to age- and country-specific prevalence.
#
# input:  Subset of m_dt that includes country and age category of male partners to whom HIV status will be assigned.
# output: Vector of 0/1 values indicating HIV status of male partners.
# 
#######################################################################################

assign_male_hiv_status <- function(dt) {

  hiv <- unlist(sapply(c("mal", "sa", "uga", "zim"), function(x) {
    dt[country == x,
       rbinom(n = nrow(dt[country == x]), size = 1, prob = male_prev[m_age, x])]
  }))
  
  return(hiv)
}