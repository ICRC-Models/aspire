#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 March 2019
# assign_adh_q1: Function to assign binary baseline adherence value according to parameters obtained in predictive model of adherence among ASPIRE trial participants (see email from EB dated 7 September 2018)
#
# input:  Subset of data in f_dt
# output: Vector of 0/1 values indicating adherence at baseline.
# 
#######################################################################################

assign_adh_q1 <- function(dt) {
  # Assign adherence value at quarter 1 based on predictive model for plasma dapivirine levels collected at Q1.
  dt[arm == 1 & visit == 3, log_odds_adh := params$pm_adhb_intercept + 
                               params$pm_adhb_age * f_age +
                               params$pm_adhb_married * b_married +
                               params$pm_adhb_edu * b_edu + 
                               params$pm_adhb_base_gon * b_gon +
                               params$pm_adhb_base_ct * b_ct +
                               params$pm_adhb_base_tr * b_trr +
                               params$pm_adhb_pknow * b_pknow +
                               params$pm_adhb_circ * b_circ +
                               params$pm_adhb_noalc * b_noalc +
                               params$pm_adhb_npart_eq_two * b_npart_eq_two +
                               params$pm_adhb_npart_gt_two * b_npart_gt_two +
                               params$pm_adhb_enrolldt * enrolldt_after_Apr12013 +
                               params$pm_adhb_Lilongwe * (site == "Lilongwe") +
                               params$pm_adhb_BothaHill * (site == "Botha's Hill") +
                               params$pm_adhb_EmavundleniCent * (site == "Emavundleni Cent") +
                               params$pm_adhb_eThekwini * (site == "eThekwini") +
                               params$pm_adhb_Isipingo * (site == "Isipingo") +
                               params$pm_adhb_RKKhan * (site == "RK Khan") +
                               params$pm_adhb_Tongaat * (site == "Tongaat") +
                               params$pm_adhb_Umkomaas * (site == "Umkomaas") +
                               params$pm_adhb_Verulam * (site == "Verulam") +
                               params$pm_adhb_WHRI * (site == "WHRI") +
                               params$pm_adhb_MUJHU * (site == "MU-JHU") +
                               params$pm_adhb_SekeSouth * (site == "Seke South") +
                               params$pm_adhb_Spilhaus * (site == "Spilhaus") +
                               params$pm_adhb_Zengeza * (site == "Zengeza")]
  
  dt[arm == 1 & visit == 3, prob_adh := exp(log_odds_adh)/(1 + exp(log_odds_adh))]
  
  dt[arm == 1 & visit == 3, adh := rbinom(n = nrow(dt[arm == 1 & visit == 3]), size = 1, prob = prob_adh)]
  
  dt[arm == 1, adh := na.locf(adh, fromLast = T), by = id]

  return(dt[, adh])
}

#######################################################################################
#
# assign_adh_fup: Function to assign binary follow-up adherence value according to parameters obtained in predictive model of adherence among ASPIRE trial participants (see email from EB dated 7 September 2018)
#
# input:  Subset of data in f_dt
# output: Vector 0/1 values indicating adherence at time == t
#
#######################################################################################

assign_adh_fup <- function(dt, t) {
  # Assigned adherence value at visits 2+ is based on predictive model for plasma dapivirine levels collected after Q1.
  dt[, lag_adh := shift(adh, n = 3), by = id]
  
  dt[arm == 1 & visit == t, log_odds_adh := params$pm_adhfu_intercept + 
                               params$pm_adhfu_lagadh1 * lag_adh +
                               params$pm_adhfu_visitnum * visit +
                               params$pm_adhfu_lagvisit_date_ind * prev_visit_date_after_Aug12013 +
                               params$pm_adhfu_base_age * f_age +
                               params$pm_adhfu_base_married * b_married +
                               params$pm_adhfu_base_edu * b_edu + 
                               params$pm_adhfu_base_gon * b_gon +
                               params$pm_adhfu_base_ct * b_ct +
                               params$pm_adhfu_base_tr * b_trr +
                               params$pm_adhfu_base_pknow * b_pknow +
                               params$pm_adhfu_base_circ * b_circ +
                               params$pm_adhfu_base_noalc * b_noalc +
                               params$pm_adhfu_base_npart_eq_two * b_npart_eq_two +
                               params$pm_adhfu_base_npart_gt_two * b_npart_gt_two +
                               params$pm_adhfu_enrolldt * enrolldt_after_Apr12013 +
                               params$pm_adhfu_Lilongwe * (site == "Lilongwe") +
                               params$pm_adhfu_BothaHill * (site == "Botha's Hill") +
                               params$pm_adhfu_EmavundleniCent * (site == "Emavundleni Cent") +
                               params$pm_adhfu_eThekwini * (site == "eThekwini") +
                               params$pm_adhfu_Isipingo * (site == "Isipingo") +
                               params$pm_adhfu_RKKhan * (site == "RK Khan") +
                               params$pm_adhfu_Tongaat * (site == "Tongaat") +
                               params$pm_adhfu_Umkomaas * (site == "Umkomaas") +
                               params$pm_adhfu_Verulam * (site == "Verulam") +
                               params$pm_adhfu_WHRI * (site == "WHRI") +
                               params$pm_adhfu_MUJHU * (site == "MU-JHU") +
                               params$pm_adhfu_SekeSouth * (site == "Seke South") +
                               params$pm_adhfu_Spilhaus * (site == "Spilhaus") +
                               params$pm_adhfu_Zengeza * (site == "Zengeza")]
  
  dt[arm == 1 & visit == t, prob_adh := exp(log_odds_adh)/(1 + exp(log_odds_adh))]
  
  dt[arm == 1 & visit == t, adh := rbinom(n = nrow(dt[arm == 1 & visit == t]), size = 1, prob = prob_adh)]
  
  dt[arm == 1, adh := na.locf(adh, fromLast = T), by = id]
  
  return(dt[visit == t | visit == t - 1 | visit == t - 2, adh])
}

#######################################################################################
#
# assign_adh_secondary: Function to assign categorical adherence value at baseline and follow-up according to systematically varied proportion of women with consistent, inconsistent, and no adherence.
#
# input:  Subset of data in f_dt
# output: Vector 1/2/3 values indicating adherence at time == t
#
#######################################################################################

assign_adh_secondary <- function(dt, prop_full_adh, prop_partial_adh, prop_non_adh) {
  
  # 1 = consistent_adh, 2 = inconsistent_adh, 3 = no_adh
  dt[, adh := sample(x = 1:3, size = 1, replace = T, prob = c(prop_full_adh, prop_partial_adh, prop_non_adh)), by = id]
  
  return(dt[, adh])
}
