#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   22 October 2019
# male_hiv_incidence: Function to assign new HIV infections
#
# input:  m_dt, baseline HIV incidence rate, relative risk of HIV positivity given female report of condom use in prior week, and current timestep
# output: Modified m_dt with updated infection status and time of infection.
# 
#######################################################################################

male_hiv_incidence <- function(m_dt, base_male_hiv_incidence, cond_rr, t_curr) {
  
  # Convert annual incidence rate to monthly probability
  monthly_inc_prob <- 1 - exp(-base_male_hiv_incidence * 1/12)
  
  # Create variable for probability of incident infection. Probabilities are country-specific and modified by age-dependent relative risk of having an HIV-positive partner and relative risk of having an HIV-positive partner given reported condom use in the prior week.
  m_dt[country == "mal", prob_inf := (monthly_inc_prob * params$mal_inc_ratio) * m_hiv_rr * ifelse(b_condom_lweek, cond_rr, 1)]
  m_dt[country == "sa",  prob_inf :=  monthly_inc_prob * m_hiv_rr * ifelse(b_condom_lweek, cond_rr, 1)]
  m_dt[country == "uga", prob_inf := (monthly_inc_prob * params$uga_inc_ratio) * m_hiv_rr * ifelse(b_condom_lweek, cond_rr, 1)]
  m_dt[country == "zim", prob_inf := (monthly_inc_prob * params$zim_inc_ratio) * m_hiv_rr * ifelse(b_condom_lweek, cond_rr, 1)]
  
  # Identify incident infections among men proportional to prob_inf, irrespective of active status nor current HIV status.
  m_dt[, inc_inf := rbinom(n    = nrow(m_dt),
                           size = 1,
                           prob = prob_inf)]
  
  # Among those currently uninfected, in an active relationship, and with a new infection, change HIV status to positive and assign time of infection as current time step.
  m_dt[active == 1 & hiv == 0 & inc_inf == 1, `:=`(hiv = as.integer(1), t_inf = t_curr)]
  
  # Assign viral load
  m_dt[t_inf == t_curr, vl := assign_male_vl_trt_unknown(dt = m_dt[t_inf == t_curr, .(m_age, hiv)])]
  
  # Remove unneeded columns
  m_dt[, `:=`(prob_inf = NULL, inc_inf = NULL)]
  
  return(m_dt)
}