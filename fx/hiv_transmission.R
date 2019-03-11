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

hiv_transmission <- function(f_dt, m_dt, t, l, rr_ai, rr_ring_full_adh, rr_ring_partial_adh) {
  
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

  # Estimate per-act risk of not acquiring HIV
  dt[, risk_no_inf := (1 - l) ^ exp(log(params$rr_vl) * (vl - 4.0) +
                                    log(params$rr_age) * (f_age - 30) +
                                    log(rr_ring_full_adh) * (adh == 1) * arm * (1 - ai) +
                                    log(rr_ring_partial_adh) * (adh == 2) * arm * (1 - ai) +
                                    log(params$rr_condom) * condom +
                                    log(rr_ai) * ai +
                                    log(params$rr_sti) * b_n_sti +
                                    log(params$rr_bv) * b_bv)]
  
  # Assign cumulative risk of infection for vi acts and ai acts by id. Assign to appropriate time step in study_dt.
  cum_risk_inf_vi <- dt[arm == 1 & ai == 0, 1 - Reduce(f = `*`, x = risk_no_inf), by = id]
  cum_risk_inf_ai <- dt[arm == 1 & ai == 1, 1 - Reduce(f = `*`, x = risk_no_inf), by = id]
  
  # Estimate proportion of HIV risk attributable to VI and AI
  cum_risk_inf_dt <- merge(x = cum_risk_inf_vi, y = cum_risk_inf_ai, by = "id", all = T)
  setnames(cum_risk_inf_dt, old = c("V1.x", "V1.y"), new = c("cum_risk_inf_vi", "cum_risk_inf_ai"))
  cum_risk_inf_dt[is.na(cum_risk_inf_ai), cum_risk_inf_ai := 0]
  cum_risk_inf_dt[is.na(cum_risk_inf_vi), cum_risk_inf_vi := 0]
  cum_risk_inf_dt[, `:=`(prop_risk_vi = cum_risk_inf_vi/(cum_risk_inf_vi + cum_risk_inf_ai),
                         prop_risk_ai = cum_risk_inf_ai/(cum_risk_inf_vi + cum_risk_inf_ai))]
  
  # Estimate cumulative risk of acquiring HIV over all acts
  risk <- dt[, 1 - Reduce(f = `*`, x = risk_no_inf), by = id]
  
  # Merge f_dt and risk data.table
  f_dt <- merge(x = f_dt, y = risk, by = "id", all.x = T)
  setnames(x = f_dt, old = "V1", new = "risk")
  
  # Merge f_dt and cum_risk_inf_dt (for those in placebo arm, there will be a value for risk, but not for prop_risk_ai)
  f_dt <- merge(x = f_dt, y = cum_risk_inf_dt[, .(id, prop_risk_ai, prop_risk_vi)], by = "id", all.x = T)
  f_dt[is.na(prop_risk_vi), prop_risk_vi := 0]
  f_dt[is.na(prop_risk_ai), prop_risk_ai := 0]

  # Assign HIV transmission as a binomial draw with probability equal to total monthly risk
  f_dt[is.na(risk), hiv_t := as.integer(0)] # Women without an HIV-positive partner have 0 risk for the month, and correspondingly do not acquire HIV.
  f_dt[is.na(hiv_t), hiv_t := rbinom(n = nrow(f_dt[is.na(hiv_t)]), size = 1, prob = risk)]
  
  return(list(f_dt$hiv_t, f_dt$prop_risk_vi, f_dt$prop_risk_ai))
}