#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 March 2018
# run_sim: function to run a specified number of simulations of the ASPIRE trial
#
# input:  None
# output: A longitudinal dataset
#
#######################################################################################

run_sim <- function() {

  ## Create vector of participant ids repeated once for each partner
  id <- f_dt[, rep(id, n_part)]
  
  ## Create table with one row for each male partner
  m_dt <- data.table(id      = id,
                     country = f_dt$country[id],
                     f_age   = f_dt$f_age[id],
                     active  = as.integer(1))
  
  ## Assign male partner age, HIV status, and viral load values.
  m_dt[, m_age := assign_male_age(m_ids = id, f_age = f_age)]
  m_dt[, hiv   := assign_male_hiv_status(dt = m_dt[, .(country, m_age)])]
  m_dt[, vl    := assign_male_vl(hiv = m_dt$hiv, age = m_dt$m_age)]
  
  ## Append fixed characteristics to f_dt that will vary in simulations
  f_dt[, `:=`(prob_ai     = assign_prob_ai(n = nrow(f_dt), prop_ai = params$prop_ai, alpha = params$alpha, beta = params$beta),
              prob_condom = assign_prob_condom(countries = f_dt[, country]))]
  
  setkey(f_dt, id)
  
  ## Create study data table
  time_steps <- seq(0, 30 * 12 * 3, 30)
  study_dt <- setorder(as.data.table(expand.grid(f_dt$id, time_steps)))
  setnames(study_dt, old = c("Var1", "Var2"), new = c("id", "time"))
  
  ## Add baseline variable values
  study_dt[time == 0, `:=`(f_age = f_dt$f_age,
                           adh   = as.integer(assign_adh_t0(dt = f_dt[, .(id, f_age, arm, site)])),
                           hiv   = as.integer(0))]
  
  setkey(study_dt, id, time)

  study_dt <- merge(x = study_dt, y = f_dt[, .SD, .SDcols = names(f_dt)[which(names(f_dt) != "f_age")]])
  
  ## Assign study retention status according to observed number of days of follow-up
  ## TO DO: Ask EB: values of fu_days do not align with study visit schedule. How was the exact number of days of follow-up determined when censoring occurred between study visits?
  study_dt[, on_study := ifelse(time <= fu_days, 1, 0)]
  
  ## Run simulation
  for(t in time_steps[2:37]) {
    
    study_dt[time == t, `:=`(f_age = aging(prev_ages = study_dt[time == t - 30, f_age]),
                             adh   = assign_adh_fup(s_dt = study_dt[time %in% c(t, t - 30), .(id, time, arm, adh)], f_dt = f_dt[, .(id, prob_condom)], t = t),
                             hiv   = hiv_transmission(s_dt = study_dt[time == t - 30, .(id, adh, hiv)],
                                                      m_dt = m_dt[active == 1 & hiv == 1, .(id, vl)],
                                                      f_dt = f_dt[, .(id, f_age, pp_acts, sti, arm, prob_ai, prob_condom)],
                                                      t   = t,
                                                      l   = params$lambda))]
    
    m_dt <- partner_change(m_dt = m_dt, f_dt = f_dt[, .(id, country, f_age, max_part)])
    
    # make run_sim a function that saves long dataset to be used in Cox proportional hazards models. Then can call run_sim within other functions that will save other output (e.g., for ABC, sum infections by arm by age)
  }
  
  ## Clean up: Remove all rows where on_study == 0
  study_dt <- study_dt[on_study == 1]
  
  ## Clean up: Among seroconverters, remove rows subsequent to visit at which HIV is first detected
  study_dt <- study_dt[, .SD[cumsum(hiv) <= 1], by = id]
  
  return(study_dt)
  
  ## TO DO: Check that male partner HIV prevalence distribution assigned in model matches observed age- and country-specific prevalence distribution
}
