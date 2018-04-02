#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   30 March 2018
# assign_m_hiv_rr: Function to assign a reduced risk of contacting an HIV-positive male partner given female age. Function is similar to an inverse logistic function, but scale parameter differs by age less than 28 and age greater than 28.
#
# input:  Female age.
# output: Vector of relative risk values that reduce a woman's probability of having an HIV-positive male partner relative to observed prevalence distribution in the general population.
# 
#######################################################################################

assign_m_hiv_rr <- function(age, l) {
  if(age <= 26) {
    rr <- 0 # l/(1 + exp(-2 * (age - 27)))
  } else {
    rr <- l/(1 + exp(-0.45 * (age - 27)))
  }
  return(1 - rr)
}
