#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 November 2018
# impute_unplw: Function to impute a value for unprotected sex in the last week at baseline using predictive model of baseline ASPIRE data.
#
# input:  Subset of f_dt.
# output: Vector of 0/1 values indicating if a woman used a condom for all sex acts in the previous week.
# 
# impute_anal_n_acts: Function to impute the number of anal sex acts. First, impute any anal sex at months baseline and month 3, given baseline characteristics. Second, among those with imputed value of any anal sex, assign number of acts as a random draw from the distribution of acts reported by all other participants who engaged in AI.
#
# input:  Subset of f_dt
# output: Vector of values indicating the number of anal sex acts a woman would have reported at baseline and month 3.
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


impute_anal_n_acts <- function(dt) {
  # 18 observations where transactional sex is missing. Rather than creating another predictive model for transactional sex, assign value with marginal probability of transactional sex (6.2%, per primary Baeten, et al., 2016 paper).
  
  dt[is.na(b_trans_sex), b_trans_sex := as.double(rbinom(n = nrow(dt[is.na(b_trans_sex)]), size = 1, prob = 0.062))]

  dt[, log_odds_any_ai := params$pm_ai_intercept +
                          params$pm_ai_Lilongwe * (site == "Lilongwe") +
                          params$pm_ai_BothaHill * (site == "Botha's Hill") +
                          params$pm_ai_EmavundleniCent * (site == "Emavundleni Cent") +
                          params$pm_ai_eThekwini * (site == "eThekwini") +
                          params$pm_ai_Isipingo * (site == "Isipingo") +
                          params$pm_ai_RKKhan * (site == "RK Khan") +
                          params$pm_ai_Tongaat * (site == "Tongaat") +
                          params$pm_ai_Umkomaas * (site == "Umkomaas") +
                          params$pm_ai_Verulam * (site == "Verulam") +
                          params$pm_ai_WHRI * (site == "WHRI") +
                          params$pm_ai_MUJHU * (site == "MU-JHU") +
                          params$pm_ai_SekeSouth * (site == "Seke South") +
                          params$pm_ai_Spilhaus * (site == "Spilhaus") +
                          params$pm_ai_Zengeza * (site == "Zengeza") +
                          params$pm_ai_edu * b_edu +
                          params$pm_ai_bdepo * b_depo +
                          params$pm_ai_bneten * b_neten +
                          params$pm_ai_trans_sex * b_trans_sex]
  
  dt[, prob_any_ai := exp(log_odds_any_ai)/(1 + exp(log_odds_any_ai))]
  
  dt[, any_ai := rbinom(n = nrow(dt), size = 1, prob = prob_any_ai)]
  
  # Assign value of number of acts. For those with no AI, 0 AI acts. For those with any AI, randomly select a value from the distribution of number of AI acts reported by other women engaging in AI at the same time point.
  dt[any_ai == 0, anal_n_acts := 0]
  
  dt[any_ai == 1 & visit == 0, anal_n_acts := sample(x = f_dt[anal_n_acts > 0 & visit == 0, anal_n_acts], size = nrow(dt[any_ai == 1 & visit == 0]), replace = T)]
  
  dt[any_ai == 1 & visit == 3, anal_n_acts := sample(x = f_dt[anal_n_acts > 0 & visit == 3, anal_n_acts], size = nrow(dt[any_ai == 1 & visit == 3]), replace = T)]
  
  return(dt[, anal_n_acts])
}
