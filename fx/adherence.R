#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   18 October 2017
# assign_adh_t0: Function to assign binary baseline adherence value according to parameters obtained in Hidden Markov model of adherence among ASPIRE trial participants
#
# input:  Subset of data in f_dt, including variables id, bl_age, arm, blantyre, lilongwe, umkomaas, isipingo, tongaat, emavundleni
# output: Vector of 0/1 values indicating adherence at baseline.
# 
# assign_adh_fup: Function to assign binary follow-up adherence value according to parameters obtained in Hidden Markov model of adherence among ASPIRE trial participants
#
# input:  Subset of data in study_dt and f_dt
# output: Vector 0/1 values indicating adherence at time == t
#######################################################################################

## Ask EB: I assign adherence at time 0. Adherence, however, is used in the model for HIV transmission to indicate adherence at the previous 30 days. As such, adherence at t0 is not used except as it sets a trajectory for adherence at subsequent time steps. Is adherence at "baseline" actually adherence at the first follow-up visit (i.e., time == 30)?

assign_adh_t0 <- function(dt) {
  # Model uses logistic regression. At baseline: log(prob_adh/(1 - prob_adh)) = -0.567 + 0.0835*age - 0.919*blantyre + 0.978*lilongwe - 1.27*umkomaas - 3.22*isipingo - 2.82*tongaat - 1.22*emavundleni
  dt[arm == 1, log_odds_adh := -0.567 + 0.0835 * f_age - 0.919 * as.numeric(site == "Blantyre") + 0.978 * as.numeric(site == "Lilongwe") - 1.27 * as.numeric(site == "Umkomaas") - 3.22 * as.numeric(site == "Isipingo") - 2.82 * as.numeric(site == "Tongaat") - 1.22 * as.numeric(site == "Emavundleni Cent")]
  dt[arm == 1, prob_adh := exp(log_odds_adh)/(1 + exp(log_odds_adh))]
  dt[arm == 1, adh := rbinom(n = nrow(dt[arm == 1]), size = 1, prob = prob_adh)]

  return(dt$adh)
}


assign_adh_fup <- function(s_dt, f_dt, t) {
  ## TO DO: Make adherence model coefficients parameters to be passed to function.
  ## TO DO: Ask EB: Is this the best way to model adherence? Should we use observed adherence values from follow-up for participants?
  
  ## Below method of using coefficients from Jingyang's Markov model did not work properly. Rather, proportion adherent for all age groups converged to 0.7. Try modifying the function to use baseline values above as intercepts, and fixed probability of remaining adherent for all participants.

  # Assign unprotected sex in the previous week as random binomial draw with probability = condom_prob. This is disconnected from protection of acts calculated in hiv_transmission.
  dt <- merge(x = s_dt, y = f_dt, by = "id", all.x = T)
  # dt[time == t, unprotect_sex := rbinom(n = nrow(dt[time == t]), size = 1, prob = prob_condom)]
  # 
  # # Model for non-adherent to adherent: log_odds_adherent = b0 + b1 * (unprotected sex in previous week), where b0 = -2.45, [-2.95, -1.94] and b1 = -0.252, [-0.485, -0.022]. Per Jingyang's email dated 2017-10-02, no other predictors in model were significant.
  # dt[which(time == t - 30 & adh == 0) + 1, log_odds_adh := -2.45 - 0.252 * unprotect_sex]
  # 
  # # Model for remaining at adherent: log_odds_adherent = b0, where b0 = 3.44, [3.02, 3.91]. Per Jingyang's email dated 2017-10-02, no other predictors in model were significant.
  # dt[which(time == t - 30 & adh == 1) + 1, log_odds_adh := 3.44]
  # 
  # # Expit to convert log-odds of adherence to probability of adherence
  # dt[time == t, prob_adh := exp(log_odds_adh)/(1 + exp(log_odds_adh))]
  
  # Adherence at time t determined as random binomial draw with probability equal to prob_adh
  dt[which(time == t - 30 & adh == 0) + 1, prob_adh := 0]
  dt[which(time == t - 30 & adh == 1) + 1, prob_adh := 0.995]
  dt[time == t & arm == 1, adh := rbinom(n = nrow(dt[time == t & arm == 1]), size = 1, prob = prob_adh)]

  return(dt[time == t, adh])
}
