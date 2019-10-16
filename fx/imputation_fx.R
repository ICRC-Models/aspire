#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 November 2018
# impute_unplw: Function to impute a value for unprotected sex in the last week at baseline using predictive model of baseline ASPIRE data.
#
# input:  Subset of f_dt.
# output: Vector of 0/1 values indicating if a woman used a condom for all sex acts in the previous week.
# 
#######################################################################################

impute_unplw <- function(dt) {
  
  dt[, log_odds_unplw := params$pm_unpsex_intercept +
                         params$pm_unpsex_basesppositive * (m_hiv == "positive") +
                         params$pm_unpsex_basespnegative * (m_hiv == "negative") +
                         params$pm_unpsex_basespdk * (m_hiv == "unknown") +
                         params$pm_unpsex_basestd * b_sti +
                         params$pm_unpsex_noalc * b_noalc +
                         params$pm_unpsex_age * f_age + 
                         params$pm_unpsex_married * b_married]
  
  dt[, prob_unplw := exp(log_odds_unplw)/(1 + exp(log_odds_unplw))]
  
  dt[, unplw := rbinom(n = nrow(dt), size = 1, prob = prob_unplw)]
  
  dt[, b_condom_lweek := as.numeric(!unplw)]
  
  return(dt[, b_condom_lweek])
}
