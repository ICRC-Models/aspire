#######################################################################################
# 
# Kathryn Peebles
# 2017-10-18
# assign_adh_t0: function to assign binary baseline adherence value according to parameters obtained in Hidden Markov model of adherence among ASPIRE trial participants
#
# input:  subset of data in f_dt, including variables id, bl_age, arm, blantyre, lilongwe, umkomaas, isipingo, tongaat, emavundleni
# output: vector of 0/1 values indicating adherence at baseline.
# 
# assign_adh_fup: function to assign binary follow-up adherence value according to parameters obtained in Hidden Markov model of adherence among ASPIRE trial participants
#
# input:  
#######################################################################################

## Ask EB: I assign adherence at time 0. Adherence, however, is used in the model for HIV transmission to indicate adherence at the previous 30 days. As such, adherence at t0 is not used except as it sets a trajectory for adherence at subsequent time steps. Is adherence at "baseline" actually adherence at the first follow-up visit (i.e., time == 30)?

assign_adh_t0 <- function(dt) {
  # Model uses logistic regression. At baseline: log(prob_adh/(1 - prob_adh)) = -0.567 + 0.0835*age - 0.919*blantyre + 0.978*lilongwe - 1.27*umkomaas - 3.22*isipingo - 2.82*tongaat - 1.22*emavundleni
  dt[arm == 1, log_odds_adh := -0.567 + 0.0835 * age - 0.919 * blantyre + 0.978 * lilongwe - 1.27 * umkomaas - 3.22 * isipingo - 2.82 * tongaat - 1.22 * emavundleni]
  dt[arm == 1, prob_adh := exp(log_odds_adh)/(1 + exp(log_odds_adh))]
  dt[arm == 1, adh := rbinom(n = nrow(dt[arm == 1]), size = 1, prob = prob_adh)]

  return(dt$adh)
}


assign_adh_fup <- function(dt, t) {
  ## TO DO: Make adherence model coefficients parameters to be passed to function.

  # Assign unprotected sex in the previous week as random binomial draw with probability = (number of monthly unprotected acts)/(monthly total acts). Because sex acts are in monthly increments, this assumes that the probability of unprotected sex is equal across the whole month (i.e., that protection at one act is independent of protection at subsequent acts).
  # If n_acts at t is 0, then probability of unprotected sex at t is 0
  prob_unprotect_sex <- dt[time == t, ifelse(test = round(n_acts) == 0,
                                             yes  = 0,
                                             no   = sum(n_acts_vi_no_condom, n_acts_ai_no_condom)/round(n_acts)), by = id]
  prob_unprotect_sex[V1 > 1, V1 := 1] # Rarely, rounding of acts results in prob_unprotect_sex > 1
  dt[time == t, unprotect_sex := rbinom(n = nrow(prob_unprotect_sex), size = 1, prob = prob_unprotect_sex$V1)]

  # Model for non-adherent to adherent: log_odds_adherent = b0 + b1 * (unprotected sex in previous week), where b0 = -2.45, [-2.95, -1.94] and b1 = -0.252, [-0.485, -0.022]. Per Jingyang's email dated 2017-10-02, no other predictors in model were significant.
  dt[which(time == t - 30 & adh == 0) + 1, log_odds_adh := -2.45 - 0.252 * unprotect_sex]
  
  # Model for remaining at adherent: log_odds_adherent = b0, where b0 = 3.44, [3.02, 3.91]. Per Jingyang's email dated 2017-10-02, no other predictors in model were significant.
  dt[which(time == t - 30 & adh == 1) + 1, log_odds_adh := 3.44]
  
  # Expit to convert log-odds of adherence to probability of adherence
  dt[time == t, prob_adh := exp(log_odds_adh)/(1 + exp(log_odds_adh))]
  
  # Adherence at time t determined as random binomial draw with probability equal to prob_adh
  dt[time == t & !is.na(prob_adh), adh := rbinom(n = nrow(dt[time == t & !is.na(prob_adh)]), size = 1, prob = prob_adh)]

  # Among pregnant women, adherence is 0 for months 2-9 and at first month following pregnancy. At second month following pregnancy, adherence is assigned as the same value as at month 1 of pregnancy (i.e., at month of pregnancy detection). This assumes that pregnancy does not alter a woman's adherence pattern.
  dt[which(time == t - 30 & preg == 1 & arm == 1) + 1, adh := 0]
  if(t > 300) {
    ids_reinit <- intersect(dt[time == t - 30 & preg == 0 & arm == 1, id], dt[time == t - 60 & preg == 1 & arm == 1, id])
    dt[time == t & id %in% ids_reinit, adh := dt[time == t - 300 & id %in% ids_reinit, adh]]
  }

  return(dt[time == t, adh])
}
