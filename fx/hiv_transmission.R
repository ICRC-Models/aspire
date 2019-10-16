#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   15 March 2018
# hiv_transmission: Function to calculate HIV transmission risk across all acts with all partners in monthly time intervals. HIV transmission is assigned according to binomial draw with probability equal to cumulative risk.
#
# input:  Subsets of m_dt, f_dt
# output: Vector of 0/1 values indicating HIV status
# 
#######################################################################################

hiv_transmission <- function(f_dt, m_dt, t, l, rr_ai, rr_ring, ri_dt = NULL, b_dt = NULL) {
  
  dt <- merge(x = m_dt, y = f_dt, by = "id", all.y = F)
  
  # Expand table to include one row per act
  dt_vi <- dt[rep(1:.N, pp_vi_acts)][, ai := 0]
  dt_ai <- dt[rep(1:.N, pp_ai_acts)][, ai := 1]
  dt <- rbindlist(l = list(dt_vi, dt_ai))
  setkey(dt, id)
  dt[, act_id := 1:.N, by = id]
  
  # Assign condom use (condom/no condom) to each act.
  dt[, condom := rbinom(n = 1:nrow(dt), size = 1, prob = prob_condom)]
  
  # Estimate per-act risk of acquiring HIV for all co-factors
  dt[, risk_inf := 1 - ((1 - l) ^ exp(log(params$rr_vl) * (vl - 4.0) +
                                      log(params$rr_age) * (f_age - 30) +
                                      log(params$rr_condom) * condom +
                                      log(rr_ai) * ai +
                                      log(params$rr_sti) * b_n_sti +
                                      log(params$rr_bv) * b_bv))]
  
  # Add column that copies risk_inf to indicate risk of infection per act prior to application of ring
  dt[, risk_inf_pre_ring := risk_inf]
  
  b_dt[time == t & arm == 0, c("median_age", "median_vl", "median_prob_condom", "mean_sti", "prop_bv", "n_partners_hiv") := dt[arm == 0, .(median(f_age), median(vl), mean(condom), mean(b_n_sti), mean(b_bv), as.double(length(unique(id))))]]
  b_dt[time == t & arm == 0, n_acts := dt[arm == 0, as.double(sum(max(act_id))), by = id][, sum(V1)]]
  
  b_dt[time == t & arm == 1, c("median_age", "median_vl", "median_prob_condom", "mean_sti", "prop_bv", "n_partners_hiv") := dt[arm == 1, .(median(f_age), median(vl), mean(condom), mean(b_n_sti), mean(b_bv), as.double(length(unique(id))))]]
  b_dt[time == t & arm == 1, n_acts := dt[arm == 1, as.double(sum(max(act_id))), by = id][, sum(V1)]]
  
  # For each vaginal act among women in the active arm with high adherence, multiply risk of infection by ring RR
  dt[ai == 0 & adh == 1 & arm == 1, risk_inf := risk_inf * rr_ring]
  
  # Determine HIV infection at each act as a random binomial draw with probability equal to the risk of infection in that act. Also calculate the number of HIV infections in the absence of the ring, with risk specified by risk_inf_pre_ring. 
  dt[, hiv_act     := rbinom(n = nrow(dt), size = 1, prob = risk_inf)]
  dt[, hiv_no_ring := rbinom(n = nrow(dt), size = 1, prob = risk_inf_pre_ring)]
  
  # Determine if a woman acquired HIV: yes if she acquired HIV at any act, no if she did not acquire HIV at any act
  dt[, hiv_t := as.integer(any(hiv_act == 1)), by = id]
  
  # By arm, calculate the total number of HIV infections and total number of acts. Also calculate the number of HIV infections prior to ring application to evaluate the imbalance in arms.
  ri_dt[time == t & arm == 0, c("n_inf", "n_inf_pre_ring", "n_acts") := .(dt[arm == 0, sum(hiv_act)], dt[arm == 0, sum(hiv_no_ring)], sum(dt[arm == 0, .(n_acts = max(act_id)), by = id][, n_acts]))]
  ri_dt[time == t & arm == 1, c("n_inf", "n_inf_pre_ring", "n_acts") := .(dt[arm == 1, sum(hiv_act)], dt[arm == 1, sum(hiv_no_ring)], sum(dt[arm == 1, .(n_acts = max(act_id)), by = id][, n_acts]))]
  
  # By arm, calculate the total number of HIV infections and total number of acts among women who engage exclusively in VI. Calculate among all women and among only fully adherent women.
  ri_dt[time == t & arm == 0, c("n_inf_vi_excl", "n_acts_vi_excl", "n_inf_vi_excl_adh", "n_acts_vi_excl_adh") := .(dt[arm == 0 & !(any_ai), sum(hiv_act)], sum(dt[arm == 0 & !(any_ai), .(n_acts = max(act_id)), by = id][, n_acts]), dt[arm == 0 & !(any_ai) & adh == 1, sum(hiv_act)], sum(dt[arm == 0 & !(any_ai) & adh == 1, .(n_acts = max(act_id)), by = id][, n_acts]))]
  ri_dt[time == t & arm == 1, c("n_inf_vi_excl", "n_acts_vi_excl", "n_inf_vi_excl_adh", "n_acts_vi_excl_adh") := .(dt[arm == 1 & !(any_ai), sum(hiv_act)], sum(dt[arm == 1 & !(any_ai), .(n_acts = max(act_id)), by = id][, n_acts]), dt[arm == 1 & !(any_ai) & adh == 1, sum(hiv_act)], sum(dt[arm == 1 & !(any_ai) & adh == 1, .(n_acts = max(act_id)), by = id][, n_acts]))]
  
  # By arm, calculate the total number of HIV infections and total number of acts among women who engage in any AI. If no women engage in AI, assign values of 0 to number of HIV infections and acts.
  if(dt[arm == 0, any(any_ai)]) {
    ri_dt[time == t & arm == 0, c("n_inf_any_ai", "n_acts_any_ai") := .(dt[arm == 0 & any_ai, sum(hiv_act)], sum(dt[arm == 0 & any_ai, .(n_acts = max(act_id)), by = id][, n_acts]))]
  } else {
    ri_dt[time == t & arm == 0, c("n_inf_any_ai", "n_acts_any_ai") := .(0, 0)]
  }
  if(dt[arm == 1, any(any_ai)]) {
    ri_dt[time == t & arm == 1, c("n_inf_any_ai", "n_acts_any_ai") := .(dt[arm == 1 & any_ai, sum(hiv_act)], sum(dt[arm == 1 & any_ai, .(n_acts = max(act_id)), by = id][, n_acts]))]
  } else {
    ri_dt[time == t & arm == 1, c("n_inf_any_ai", "n_acts_any_ai") := .(0, 0)]
  }
  
  # Calculate among all women and among only fully adherent women.
  if(dt[arm == 0 & adh == 1, any(any_ai)]) {
    ri_dt[time == t & arm == 0, c("n_inf_any_ai_adh", "n_acts_any_ai_adh") := .(dt[arm == 0 & any_ai & adh == 1, sum(hiv_act)], sum(dt[arm == 0 & any_ai & adh == 1, .(n_acts = max(act_id)), by = id][, n_acts]))]
  } else {
    ri_dt[time == t & arm == 0, c("n_inf_any_ai_adh", "n_acts_any_ai_adh") := .(0, 0)]
  }
  if(dt[arm == 1 & adh == 1, any(any_ai)]) {
    ri_dt[time == t & arm == 1, c("n_inf_any_ai_adh", "n_acts_any_ai_adh") := .(dt[arm == 1 & any_ai & adh == 1, sum(hiv_act)], sum(dt[arm == 1 & any_ai & adh == 1, .(n_acts = max(act_id)), by = id][, n_acts]))]
  } else {
    ri_dt[time == t & arm == 1, c("n_inf_any_ai_adh", "n_acts_any_ai_adh") := .(0, 0)]
  }
  
  # Create a data table censored at the first HIV infection per woman to calculate the number of acts with an HIV-positive partner and cumulative risk per woman.
  censor_dt <- dt[, .SD[cumsum(cumsum(hiv_act)) <= 1], by = id]
  exposure_dt <- censor_dt[, .(n_acts_hiv_pos = max(act_id), cum_risk_pre_ring = 1 - Reduce(f = `*`, x = (1 - risk_inf_pre_ring))), by = id]
  
  # Create table to track act type at which each infection occurred (arm and adh are already tracked in f_dt)
  act_dt <- censor_dt[hiv_act == 1, .(id, inf_act_type = ifelse(ai == 1, "ai", "vi"))]
  
  # Merge act type with f_dt
  f_dt <- merge(x = f_dt, y = act_dt, by = "id", all.x = T)
  
  # Merge exposure information with f_dt (number of acts and cumulative risk)
  f_dt <- merge(x = f_dt, y = exposure_dt, by = "id", all.x = T)
  f_dt[is.na(n_acts_hiv_pos), `:=`(n_acts_hiv_pos = as.integer(0), cum_risk_pre_ring = as.double(0))] # Women without an HIV-positive partner have 0 acts with an HIV-positive partner and 0 exposure.
  
  # Merge HIV infection results with f_dt
  f_dt <- merge(x = f_dt, y = unique(dt[, .(id, hiv_t)]), by = "id", all.x = T)
  f_dt[n_acts_hiv_pos == 0, hiv_t := as.integer(0)] # Women without an HIV-positive partner have 0 risk for the month, and correspondingly do not acquire HIV.
  
  if(nrow(f_dt) < 2614) { stop("Number of participants fewer than 2,614.") }
  
  return(list(f_dt$hiv_t, f_dt$n_acts_hiv_pos, f_dt$cum_risk_pre_ring, f_dt$inf_act_type))
}
