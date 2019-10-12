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

hiv_transmission <- function(f_dt, m_dt, t, l, rr_ai, rr_ring_full_adh, rr_ring_partial_adh, ri_dt = NULL, b_dt = NULL) {
  
  dt <- merge(x = m_dt, y = f_dt, by = "id", all.y = F)
  
  # Expand table to include one row per act
  dt_vi <- dt[rep(1:.N, pp_vi_acts)][, ai := 0]
  dt_ai <- dt[rep(1:.N, pp_ai_acts)][, ai := 1]
  dt <- rbindlist(l = list(dt_vi, dt_ai))
  setkey(dt, id)
  dt[, act_id := 1:.N, by = id]
  
  # Assign condom use (condom/no condom) to each act. For estimation of risk, replace NA values of adh (currently assigned to those in placebo arm) with value of 3 (no adherence)
  dt[, `:=`(condom = rbinom(n = 1:nrow(dt), size = 1, prob = prob_condom),
            adh    = ifelse(test = is.na(adh), yes = 3, no = adh))]
  
  # Estimate per-act risk of acquiring HIV for all co-factors
  dt[, risk_inf := 1 - ((1 - l) ^ exp(log(params$rr_vl) * (vl - 4.0) +
                                        log(params$rr_age) * (f_age - 30) +
                                        log(params$rr_condom) * condom +
                                        log(rr_ai) * ai +
                                        log(params$rr_sti) * b_n_sti +
                                        log(params$rr_bv) * b_bv))]
  
  # Add column that copies risk_inf to indicate risk of infection per-act prior to application of ring
  dt[, risk_inf_pre_ring := risk_inf]
  
  b_dt[time == t & arm == 0, c("median_age", "median_vl", "median_prob_condom", "mean_sti", "prop_bv", "n_partners_hiv") := dt[arm == 0, .(median(f_age), median(vl), mean(condom), mean(b_n_sti), mean(b_bv), as.double(length(unique(id))))]]
  b_dt[time == t & arm == 0, n_acts := dt[arm == 0, as.double(sum(max(act_id))), by = id][, sum(V1)]]
  
  b_dt[time == t & arm == 1, c("median_age", "median_vl", "median_prob_condom", "mean_sti", "prop_bv", "n_partners_hiv") := dt[arm == 1, .(median(f_age), median(vl), mean(condom), mean(b_n_sti), mean(b_bv), as.double(length(unique(id))))]]
  b_dt[time == t & arm == 1, n_acts := dt[arm == 1, as.double(sum(max(act_id))), by = id][, sum(V1)]]
  
  # Multiply risk of infection by corresponding ring RR (not adherent, partially adherent, or fully adherent; AI or VI act)
  dt[ai == 0 & adh == 1 & arm == 1, risk_inf := risk_inf * rr_ring_full_adh]
  dt[ai == 0 & adh == 2 & arm == 1, risk_inf := risk_inf * rr_ring_partial_adh]
  
  # Determine HIV infection at each act as a random binomial draw with probability equal to the risk of infection in that act. Also calculate the number of HIV infections in the absence of the ring, with risk specified by risk_inf_pre_ring. 
  dt[, hiv_act     := rbinom(n = nrow(dt), size = 1, prob = risk_inf)]
  dt[, hiv_no_ring := rbinom(n = nrow(dt), size = 1, prob = risk_inf_pre_ring)]
  
  # Determine if a woman acquired HIV: yes if she acquired HIV at any act, no if she did not acquire HIV at any act
  dt[, hiv_t := as.integer(any(hiv_act == 1)), by = id]
  
  # By arm, calculate the total number of HIV infections and total number of acts. Also calculate the number of HIV infections prior to ring application to evaluate the imbalance in arms.
  ri_dt[time == t & arm == 0, c("n_inf_total", "n_inf_pre_ring", "n_acts_total") := .(dt[arm == 0, sum(hiv_act)], dt[arm == 0, sum(hiv_no_ring)], sum(dt[arm == 0, .(n_acts = max(act_id)), by = id][, n_acts]))]
  ri_dt[time == t & arm == 1, c("n_inf_total", "n_inf_pre_ring", "n_acts_total") := .(dt[arm == 1, sum(hiv_act)], dt[arm == 1, sum(hiv_no_ring)], sum(dt[arm == 1, .(n_acts = max(act_id)), by = id][, n_acts]))]
  
  # By arm, calculate the number of HIV infections and number of acts when women are censored at their first HIV infection
  ri_dt[time == t & arm == 0, c("n_inf_one_inf_per_woman", "n_acts_one_inf_per_woman") := .(dt[arm == 0 & hiv_t == 1, length(unique(id))], sum(dt[arm == 0 & hiv_act == 1, .(n_acts = min(act_id)), by = id][, n_acts]) + sum(dt[arm == 0 & hiv_t == 0, .(n_acts = max(act_id)), by = id][, n_acts]))]
  # If there are no infections in the active arm (happens occasionally), counting acts from an empty data table will return a warning. In those instances, count acts only from those with no 
  if(dt[arm == 1 & hiv_act == 1, length(unique(id)) == 0]) {
    ri_dt[time == t & arm == 1, c("n_inf_one_inf_per_woman", "n_acts_one_inf_per_woman") := .(dt[arm == 1 & hiv_t == 1, length(unique(id))], sum(dt[arm == 1 & hiv_t == 0, .(n_acts = max(act_id)), by = id][, n_acts]))]
  } else {
    ri_dt[time == t & arm == 1, c("n_inf_one_inf_per_woman", "n_acts_one_inf_per_woman") := .(dt[arm == 1 & hiv_t == 1, length(unique(id))], sum(dt[arm == 1 & hiv_act == 1, .(n_acts = min(act_id)), by = id][, n_acts]) + sum(dt[arm == 1 & hiv_t == 0, .(n_acts = max(act_id)), by = id][, n_acts]))]
  }
  
  # Create a data table censored at the first HIV infection per woman to calculate the number of acts with an HIV-positive partner and cumulative risk per woman.
  censor_dt <- dt[, .SD[cumsum(cumsum(hiv_act)) <= 1], by = id]
  exposure_dt <- censor_dt[, .(n_acts_hiv_pos = max(act_id), cum_risk_pre_ring = 1 - Reduce(f = `*`, x = (1 - risk_inf_pre_ring))), by = id]
  
  # To calculate the proportion of risk attributable to AI and non-adherence, first calculate the per-act risk of *not* being infected
  dt[, risk_no_inf := 1 - risk_inf]
  
  # Assign cumulative risk of infection for vi acts and ai acts by id. Assign to appropriate time step in study_dt.
  cum_risk_inf_vi <- dt[arm == 1 & ai == 0, .(cum_risk_inf_vi = 1 - Reduce(f = `*`, x = risk_no_inf)), by = id]
  cum_risk_inf_ai <- dt[arm == 1 & ai == 1, .(cum_risk_inf_ai = 1 - Reduce(f = `*`, x = risk_no_inf)), by = id]
  cum_risk_inf_dt <- dt[, .(risk = 1 - Reduce(f = `*`, x = risk_no_inf)), by = id]
  
  # Estimate proportion of HIV risk attributable to VI and AI
  cum_risk_inf_dt <- merge(x = cum_risk_inf_dt, y = cum_risk_inf_vi, by = "id", all = T)
  cum_risk_inf_dt <- merge(x = cum_risk_inf_dt, y = cum_risk_inf_ai, by = "id", all = T)
  
  # Merge f_dt and cum_risk_inf_dt
  f_dt <- merge(x = f_dt, y = cum_risk_inf_dt[, .(id, risk, cum_risk_inf_ai, cum_risk_inf_vi)], by = "id", all.x = T)
  f_dt[is.na(cum_risk_inf_ai), cum_risk_inf_ai := 0]
  f_dt[is.na(cum_risk_inf_vi), cum_risk_inf_vi := 0]
  f_dt[, `:=`(prop_risk_vi = cum_risk_inf_vi/(cum_risk_inf_vi + cum_risk_inf_ai),
              prop_risk_ai = cum_risk_inf_ai/(cum_risk_inf_vi + cum_risk_inf_ai))]
  f_dt[is.na(prop_risk_vi), prop_risk_vi := 0]
  f_dt[is.na(prop_risk_ai), prop_risk_ai := 0]
  
  # Merge exposure information with f_dt (number of acts and cumulative risk)
  f_dt <- merge(x = f_dt, y = exposure_dt, by = "id", all.x = T)
  f_dt[is.na(n_acts_hiv_pos), `:=`(n_acts_hiv_pos = as.integer(0), cum_risk_pre_ring = as.double(0))] # Women without an HIV-positive partner have 0 acts with an HIV-positive partner and 0 exposure.
  
  # Merge HIV infection results with f_dt
  f_dt <- merge(x = f_dt, y = unique(dt[, .(id, hiv_t)]), by = "id", all.x = T)
  f_dt[is.na(risk), hiv_t := as.integer(0)] # Women without an HIV-positive partner have 0 risk for the month, and correspondingly do not acquire HIV.
  
  if(nrow(f_dt) < 2614) { stop("Number of participants fewer than 2,614.") }
  
  return(list(f_dt$hiv_t, f_dt$n_acts_hiv_pos, f_dt$cum_risk_pre_ring, f_dt$prop_risk_vi, f_dt$prop_risk_ai))
}