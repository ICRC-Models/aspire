#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   30 March 2018
# assign_m_hiv_rr: Function to assign a reduced risk of contacting an HIV-positive male partner given female age. Function is similar to an inverse logistic function, with c representing 1 - "carrying capacity" (i.e., 1 minus the minimum risk) and s representing a scale parameter (i.e., steepness of slope in declining risk by age).
#
# input:  Female age.
# output: Vector of relative risk values that reduce a woman's probability of having an HIV-positive male partner relative to observed prevalence distribution in the general population.
# 
#######################################################################################

assign_m_hiv_rr <- function(age, c, s) {
  if(age <= 26) {
    rr <- 0 # c/(1 + exp(-2 * (age - 27)))
  } else {
    rr <- c/(1 + exp(-s * (age - 27)))
  }
  return(1 - rr)
}
