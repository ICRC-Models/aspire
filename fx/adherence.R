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
  dt[arm == 1 & days_pre_adh_int != 0 & f_age_cat != "27-45", prob_adh := prob_adh * params$pre_adh_int_rr_bl]
  dt[arm == 1, adh := rbinom(n = nrow(dt[arm == 1]), size = 1, prob = prob_adh)]

  return(dt$adh)
}


assign_adh_fup <- function(s_dt, f_dt, t) {
  ## TO DO: Make adherence model coefficients parameters to be passed to function.
  ## TO DO: EB is estimating time-varying coefficients for probability of adherence. For now, use baselines estimates from Jingyang's Markov model, with age-group specific probability of remaining adherent based on approximate interpretation of Figure 3 in Baeten, et al., 2016
  dt <- merge(x = s_dt, y = f_dt, by = "id", all.x = T)
  
  # In Figure 3, adherence remains relatively steady among women ages 22-26 and 27-45, but declines by about 10% over the duration of the trial among women ages 18-21. However, data in this figure are from after adherence interventions were implemented in August, 2013. Therefore, assume that women have a lower probability of remaining adherent prior to initiation of adherence interventions.
  # Assume that women who are initially non-adherent remain non-adherent
  dt[which(time == t - 30 & adh == 0) + 1, prob_adh := 0]

  dt[which(time == t - 30 & adh == 1 & (f_age_cat == "22-26" | f_age_cat == "27-45")) + 1, prob_adh := 1]
  dt[which(time == t - 30 & adh == 1 & f_age_cat == "18-21") + 1, prob_adh := 0.9945]
  
  # If time prior to adherence interventions, probability of adherence is reduced by params$pre_adh_int_rr_fu
  dt[time == t & pre_adh_int == 1 & f_age_cat != "27-45", prob_adh := prob_adh * params$pre_adh_int_rr_fu]
  
  dt[time == t & arm == 1, adh := rbinom(n = nrow(dt[time == t & arm == 1]), size = 1, prob = prob_adh)]

  return(dt[time == t, adh])
}
