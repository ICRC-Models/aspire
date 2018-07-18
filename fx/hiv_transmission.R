#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   15 March 2018
# hiv_transmission: Function to calculate HIV transmission risk across all acts with all partners in monthly time intervals. HIV transmission is assigned according to binomial draw with probability equal to cumulative risk.
#
# input:  Subsets of study_dt, m_dt, f_dt
# output: Vector of 0/1 values indicating HIV status
# 
#######################################################################################

hiv_transmission <- function(s_dt, m_dt, f_dt, t, l) {
  ## TO DO: Make all of the relative risks parameters that are loaded from elsewhere.
  dt <- merge(x = m_dt, y = f_dt, by = "id", all.y = F)
  
  dt <- merge(x = dt, y = s_dt, by = "id", all.y = F)
  
  # Assign number of acts by randomly taking floor or ceiling of average monthly per-partner acts (pp_acts)
  dt[, acts := ifelse(test = rbinom(n = nrow(dt), size = 1, prob = 0.5) == 1,
                      yes  = ceiling(pp_acts),
                      no   = floor(pp_acts))]
  
  # Expand table to include one row per act
  dt <- dt[rep(1:.N, acts)][, act_id := 1:.N, by = id]
  
  # Assign to each act: type of act (ai/vi) and condom use (condom/no condom). For estimation of risk, replace NA values of adh (currently assigned to those in placebo arm) with 0
  dt[, `:=`(ai     = rbinom(n = 1:nrow(dt), size = 1, prob = prob_ai),
            condom = rbinom(n = 1:nrow(dt), size = 1, prob = prob_condom),
            adh    = ifelse(test = is.na(adh), yes = 0, no = adh))]

  # Estimate per-act risk of not acquiring HIV
  dt[, risk_no_inf := (1 - l) ^ exp(log(2.89)  * (vl - 4.0) +
                                    log(0.96)  * (f_age - 30) +
                                    log(0.25)  * adh * arm * (1 - ai) +
                                    log(0.22)  * condom +
                                    log(17.25) * ai +
                                    log(2.50)  * n_sti +
                                    log(3.63)  * bv)]
  
  # Assign cumulative risk of infection for vi acts and ai acts by id. Assign to appropriate time step in study_dt.
  cum_risk_inf_vi <- dt[arm == 1 & ai == 0, 1 - Reduce(f = `*`, x = risk_no_inf), by = id]
  cum_risk_inf_ai <- dt[arm == 1 & ai == 1, 1 - Reduce(f = `*`, x = risk_no_inf), by = id]
  
  # Estimate proportion of HIV risk attributable to vi and ai
  cum_risk_inf_dt <- merge(x = cum_risk_inf_vi, y = cum_risk_inf_ai, by = "id", all.x = T, all.y = T)
  cum_risk_inf_dt[is.na(V1.x), V1.x := 0]
  cum_risk_inf_dt[is.na(V1.y), V1.y := 0]
  cum_risk_inf_dt[, `:=`(prop_risk_vi = V1.x/(V1.x + V1.y),
                         prop_risk_ai = V1.y/(V1.x + V1.y))]
  
  # Estimate cumulative risk of acquiring HIV over all acts
  risk <- dt[, 1 - Reduce(f = `*`, x = risk_no_inf), by = id]
  
  # Merge s_dt and risk data.table
  s_dt <- merge(x = s_dt, y = risk, by = "id", all.x = T)
  setnames(x = s_dt, old = "V1", new = "risk")
  
  # Merge s_dt and cum_risk_inf_dt
  s_dt <- merge(x = s_dt, y = cum_risk_inf_dt[, .(id, prop_risk_ai)], by = "id", all.x = T)
  s_dt[is.na(prop_risk_ai), prop_risk_ai := 0]

  # If HIV-positive at previous timestep, assign HIV status of 1 at all subsequent time steps
  s_dt[hiv == 1, hiv_t := as.integer(1)]
  
  # If no HIV-positive partners at current timestep, assign HIV status from previous time step
  s_dt[is.na(risk), hiv_t := hiv]

  # Assign HIV transmission as a binomial draw with probability equal to total monthly risk
  s_dt[is.na(hiv_t), hiv_t := rbinom(n = nrow(s_dt[is.na(hiv_t)]), size = 1, prob = risk)]
  
  return(list(s_dt$hiv_t, s_dt$prop_risk_ai))
}