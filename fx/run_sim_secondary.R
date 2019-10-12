#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 March 2019
# run_sim_secondary: function to run a simulation of the ASPIRE trial
#
# input:  Specified parameters
# output: Likelihood, proportion of acts protected by adherence among women in ring arm, proportion of acts that are AI (prevalence * per-person proportion), proportion efficacy dilution attributable to non-adherence and AI
#
#######################################################################################

run_sim_secondary <- function(lambda = params$lambda, cond_rr = params$cond_rr, c = params$c, s = params$s, rr_ai, prev_ai, prop_ai, prop_full_adh, prop_partial_adh, prop_non_adh, rr_ring_full_adh, rr_ring_partial_adh, i, output_dir, debug_sim_hiv = F, debug_sim_hiv_acts = F, debug_sim_acts = F, debug_sim_hiv_risk = F, debug_rerandomize = F) {
  
  load(file = paste0(wd, "/data-private/f_dt.RData"))
  
  if(debug_rerandomize) {
    # Create new variable to indicate if an individual reports a number of average monthly acts that is greater than the median reported number of monthly acts.
    f_dt[, gt_median_acts := as.numeric(pp_vi_acts_avg * b_n_part > median(pp_vi_acts_avg * b_n_part))]
    
    dpv_ids <- re_randomize_arm(dt = f_dt[visit == 0, .(id, site, b_bv, gt_median_acts)])
    
    f_dt[id %in% dpv_ids,    arm := 1]
    f_dt[!(id %in% dpv_ids), arm := 0]
  }
  
  # Impute missing values for condom use in the last week and number of anal sex acts
  f_dt[is.na(b_condom_lweek) & visit == 0, b_condom_lweek := impute_unplw(dt = f_dt[is.na(b_condom_lweek) & visit == 0, .(id, m_hiv, b_sti, f_age, b_noalc, b_married)])]
  f_dt[, b_condom_lweek := na.locf(b_condom_lweek), by = id]
  
  f_dt[is.na(anal_n_acts) & (visit == 0 | visit == 3), anal_n_acts := impute_anal_n_acts(dt = f_dt[is.na(anal_n_acts) & (visit == 0 | visit == 3), .(id, visit, site, b_edu, b_depo, b_neten, b_trans_sex)])]
  
  # Estimate per-partner average number of monthly anal acts (total average number of monthly AI acts divided by number of partners).
  f_dt[, n_ai_monthly_avg := round(mean(anal_n_acts, na.rm = T)/3, 2), by = id]
  f_dt[, pp_ai_acts_avg := round(n_ai_monthly_avg/b_n_part, 2)]
  
  # Estimate total per-partner average number of acts (VI + AI)
  f_dt[, pp_acts_avg := pp_vi_acts_avg + pp_ai_acts_avg]
  
  # Assign proportion of AI acts to each woman
  f_dt <- merge(x = f_dt, y = assign_prob_ai(ids = f_dt[, unique(id)], prev_ai = prev_ai, prop_ai = prop_ai), by = "id", all = T)
  
  # Redistribute VI and AI acts so that the total number of acts remains approximately constant across simulations with varying prevalence and frequency of AI
  f_dt[, pp_vi_acts_avg := pp_acts_avg * (1 - prob_ai)]
  f_dt[, pp_ai_acts_avg := pp_acts_avg * prob_ai]
  
  # Assign relative probability of having an HIV-positive male partner
  f_dt[, m_hiv_rr := sapply(f_age, function(x) { assign_m_hiv_rr(age = x, c = c, s = s) })]
  
  setkey(f_dt, id, visit)
  
  # Create vector of participant ids repeated once for each partner
  id <- f_dt[visit == 0, rep(id, b_n_part)]
  
  # Create table with one row for each male partner
  m_dt <- data.table(id             = id,
                     arm            = f_dt[visit == 0]$arm[id],
                     country        = f_dt[visit == 0]$country[id],
                     f_age          = f_dt[visit == 0]$f_age[id],
                     active         = as.integer(1),
                     max_part       = f_dt[visit == 0]$max_part[id],
                     b_condom_lweek = f_dt[visit == 0]$b_condom_lweek[id],
                     m_hiv_rr       = f_dt[visit == 0]$m_hiv_rr[id])
  
  if(nrow(m_dt[max_part == 1]) != nrow(f_dt[visit == 0 & b_married == 1])) { browser() }
  
  # Assign male partner age, HIV status, and viral load values.
  m_dt[, m_age := assign_male_age(m_ids = id, f_age = f_dt[visit == 0]$f_age)]
  
  if(debug_sim_hiv) {
    # If debugging, assign equal HIV risk to all women by assigning all partners to be HIV-positive and with the same viral load. 
    m_dt[, `:=`(hiv = as.integer(1), vl = 3.0)]
    
    f_dt[,  `:=`(prob_condom = assign_prob_condom(countries = f_dt[, country]),
                 pp_vi_acts  = pp_vi_acts_avg,
                 pp_ai_acts  = pp_ai_acts_avg)]
  } else if(debug_sim_hiv_acts) { 
    # If debugging, assign equal HIV risk to all women by assigning all partners to be HIV-positive and with the same viral load. 
    m_dt[, `:=`(hiv = as.integer(1), vl = 3.0)]
    
    # The number of sexual acts and condom use should also be equal.
    f_dt[,  `:=`(prob_condom = 0,
                 pp_vi_acts  = pp_vi_acts_avg,
                 pp_ai_acts  = pp_ai_acts_avg)]
  } else if(debug_sim_acts) {
    # If debugging, assign equal HIV risk to all women by assigning all partners to be HIV-positive and with the same viral load. 
    m_dt[, `:=`(hiv = as.integer(1), vl = 3.0)]
    
    # Condom use should also be equal, but number of sex acts can vary per timestep.
    f_dt[,  `:=`(prob_condom = 0,
                 pp_vi_acts  = assign_mvi_acts(dt = f_dt[, .(id, pp_vi_acts_avg)]),
                 pp_ai_acts  = assign_mai_acts(dt = f_dt[, .(id, pp_ai_acts_avg)]))]
  } else {
    m_dt[, hiv   := assign_male_hiv_status(dt = m_dt[, .(country, m_age, b_condom_lweek, m_hiv_rr)], cond_rr = cond_rr)]
    m_dt[, vl    := assign_male_vl_trt_unknown(dt = m_dt[, .(m_age, hiv)])]
    
    # Replace randomly assigned baseline hiv and viral load values for women who report an HIV-positive partner
    m_dt[id %in% f_dt[m_hiv == "positive", id], hiv := 1]
    m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == T, id], vl := log10(50)]
    m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == F, id], 
         vl := assign_male_vl_trt_known(dt = m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == F, id], .(m_age, hiv)])]
    
    # Among married women, if reported partner HIV status is negative, set partner status to negative
    m_dt[id %in% f_dt[b_married == 1 & m_hiv == "negative" & f_age_cat == "27-45", id], `:=`(hiv = as.integer(0), vl = as.integer(0))]
    
    # Append fixed characteristics to f_dt that will vary in simulations
    # f_dt[,  `:=`(prob_condom = assign_prob_condom(countries = f_dt[, country]),
    #              pp_vi_acts  = assign_mvi_acts(dt = f_dt[, .(id, pp_vi_acts_avg)]),
    #              pp_ai_acts  = assign_mai_acts(dt = f_dt[, .(id, pp_ai_acts_avg)]))]
    
    f_dt[,  `:=`(prob_condom = prop_condom,
                 pp_vi_acts  = assign_mvi_acts(dt = f_dt[, .(id, pp_vi_acts_avg)]),
                 pp_ai_acts  = assign_mai_acts(dt = f_dt[, .(id, pp_ai_acts_avg)]))]
  }
  
  # Add variables to track through simulations: hiv and proportion of risk attributable to VI/AI
  f_dt[visit == 0, `:=`(hiv          = as.integer(0),
                        prop_risk_vi = 0,
                        prop_risk_ai = 0)]
  
  # Assign adherence
  f_dt[arm == 1, adh := assign_adh_secondary(dt = f_dt[arm == 1, .(id, visit)], prop_full_adh = prop_full_adh, prop_partial_adh = prop_partial_adh, prop_non_adh = prop_non_adh)]
  
  # Create table to track cumulative risk in the placebo and active arms. The number of time steps will be variable, so make the table very large, and then reduce size below. 
  risk_inf_dt <- as.data.table(expand.grid(time = 1:150, arm = as.double(0:1), n_inf_total = NA_real_, n_inf_pre_ring = NA_real_, n_acts_total = NA_real_, n_inf_one_inf_per_woman = NA_real_, n_acts_one_inf_per_woman = NA_real_))
  setkey(risk_inf_dt, time)
  
  # Create data table to track co-factors for HIV infection by arm and risk, prior to application of RR of ring.
  balanced_dt <- as.data.table(expand.grid(time = 1:150, arm = as.double(0:1), median_age = NA_real_, median_vl = NA_real_, median_prob_condom = NA_real_, mean_sti = NA_real_, prop_bv = NA_real_, n_partners_hiv = NA_real_, n_acts = NA_real_))
  setkey(balanced_dt, time, arm)
  
  # Initialize number of HIV infections to 0 and timestep to 1.
  n_hiv_inf <- 0
  t         <- 1
  
  # Create vector to track the number of active partnerships at each timestep
  n_active_partnerships     <- numeric(150)
  n_active_partnerships_dpv <- numeric(150)
  n_active_partnerships_pbo <- numeric(150)
  
  # Run simulation until 168 infections have occurred
  while(n_hiv_inf < 168) {
  
    # Calculate the number of HIV-positive partners each woman has
    setkey(m_dt, id) # Remove this when code to track the number of HIV-positive partners each woman has is no longer in use
    f_dt[visit == t - 1, n_part_hiv := m_dt[, sum(hiv == 1 & active == 1), by = id]$V1]
    
    # HIV transmission
    f_dt[visit == t, c("hiv", "n_acts_hiv_pos", "cum_risk_pre_ring", "prop_risk_vi", "prop_risk_ai") := hiv_transmission(
      f_dt  = f_dt[visit == t, .(id, adh, f_age, pp_vi_acts, pp_ai_acts, b_n_sti, b_bv, arm, prob_condom, prop_condom)],
      m_dt  = m_dt[active == 1 & hiv == 1, .(id, vl)],
      t     = t,
      l     = lambda,
      rr_ai = rr_ai,
      rr_ring_full_adh    = rr_ring_full_adh,
      rr_ring_partial_adh = rr_ring_partial_adh,
      ri_dt     = risk_inf_dt,
      b_dt      = balanced_dt)]
    
    # Calculate number of acts of each type (VI and AI) and protection (adherent, partially adherent, not adherent)
    f_dt[visit == t, c("n_acts_vi", "n_acts_ai", "n_acts_full_adh", "n_acts_partial_adh", "n_acts_non_adh") := calculate_acts(f_dt = f_dt[visit == t, .(id, adh, pp_vi_acts, pp_ai_acts)], m_dt = m_dt[active == 1, id])]
    
    # HIV transmission is calculated at each time step independently of HIV status at previous time steps. Track if HIV transmission has occurred at any time step in order to complete incidence and survival analyses correctly.
    f_dt[visit == t, any_hiv := as.numeric(f_dt[visit %in% 1:t, any(hiv == 1), by = id]$V1)]
    
    # Track number of active partnerships at each timestep
    n_active_partnerships[t]     <- m_dt[, sum(active)]
    n_active_partnerships_dpv[t] <- m_dt[arm == 1, sum(active)]
    n_active_partnerships_pbo[t] <- m_dt[arm == 0, sum(active)]
    
    # Partner change
    m_dt <- partner_change(m_dt = m_dt, f_dt = f_dt[visit == t, .(id, arm, country, f_age, max_part, b_condom_lweek, m_hiv_rr)], cond_rr = cond_rr, debug_sim_hiv_acts = debug_sim_hiv_acts)
    
    # Count and update number of HIV infections
    n_hiv_inf <- sum(f_dt[, any(any_hiv == 1, na.rm = T), by = id]$V1)
    
    # Increment timestep
    t <- t + 1
    
    if(t > 150) { browser() }
  }
  
  # Clean up: Remove all rows where on_study == 0
  # f_dt <- f_dt[on_study == 1] # While debugging, don't censor women at observed censor time. Keep all observation time until specified number of HIV infections occurs.
  
  # Calculate the number of HIV-positive partners by arm at each timestep, prior to censoring
  n_part_hiv_dt <- f_dt[!is.na(n_part_hiv), .(n_part_hiv = sum(n_part_hiv), pre_censor = 1), by = c("visit", "arm")]
  
  # Create table to track number of HIV-positive partners over time by id
  check_hiv_assignment_dt <- f_dt[!is.na(n_part_hiv), .(id, visit, n_part_hiv, f_age, b_condom_lweek, m_hiv_rr)]
  
  # Clean up: Among seroconverters, remove rows subsequent to visit at which HIV is first detected (hiv_transmission module takes only the values for the current time step, and so does not take account of previous HIV status values).
  f_dt[visit == 0, any_hiv := 0]
  ids_hiv <- f_dt[any_hiv == 1, unique(id)]
  f_dt <- f_dt[, .SD[cumsum(any_hiv) <= 1], by = id]
  
  if(any(f_dt[, cumsum(hiv), by = id]$V1 > 1)) { stop("HIV censoring not done correctly.") }
  if(!all(f_dt[hiv == 1, unique(id)] == ids_hiv)) { stop("HIV censoring not done correctly.") }
  
  # Simulations should run until 168 infections have occurred. Identify visit with number of cumulative infections that is closest to 168, and remove visits after that visit.
  inf_by_visit <- f_dt[, .(n_inf = sum(hiv)), by = visit]
  inf_by_visit[, cum_inf := cumsum(n_inf)]
  censor_visit <- inf_by_visit[which.min(abs(cum_inf - 168)), visit]
  f_dt <- f_dt[visit <= censor_visit]
  
  # Also remove rows from risk_inf_dt subsequent to the censor_visit and columns from n_active_partnerships
  risk_inf_dt     <- risk_inf_dt[time <= censor_visit]
  balanced_dt     <- balanced_dt[time <= censor_visit]
  n_active_partnerships     <- n_active_partnerships[1:censor_visit]
  n_active_partnerships_dpv <- n_active_partnerships_dpv[1:censor_visit]
  n_active_partnerships_pbo <- n_active_partnerships_pbo[1:censor_visit]
  
  exp_dt <- f_dt[n_acts_hiv_pos > 0, .(exposures = sum(n_acts_hiv_pos), pre_ring_risk = mean(cum_risk_pre_ring)), by = c("visit", "arm")]
  
  # Calculate the number of HIV-positive partners by arm at each timestep, after censoring
  n_part_hiv_dt <- rbindlist(l = list(n_part_hiv_dt, f_dt[!is.na(n_part_hiv), .(n_part_hiv = sum(n_part_hiv), pre_censor = 0), by = c("visit", "arm")]))
  
  # Calculate proportion of acts of each type
  prop_total_acts_ai <- f_dt[, sum(n_acts_ai, na.rm = T)/sum(n_acts_ai, n_acts_vi, na.rm = T)]
  prev_ai_sim <- sum(f_dt[, any(n_acts_ai > 0, na.rm = T), by = id]$V1)/f_dt[, length(unique(id))]
  prop_ai_sim <- f_dt[id %in% f_dt[, any(n_acts_ai > 0, na.rm = T), by = id][V1 == T, id], sum(n_acts_ai, na.rm = T)/(sum(n_acts_ai, n_acts_vi, na.rm = T))]
  
  # Calculate proportion of ring-protected acts among women in dapivirine arm as the proportion acts in which women are fully adherent + the proportion of acts in which women are partially adherent * the ratio of protection from partial adherence to protection from full adherence.
  prop_acts_full_adh    <- f_dt[arm == 1, sum(n_acts_full_adh, na.rm = T)/(sum(n_acts_full_adh, n_acts_partial_adh, n_acts_non_adh, na.rm = T))]
  prop_acts_partial_adh <- f_dt[arm == 1, sum(n_acts_partial_adh, na.rm = T)/(sum(n_acts_full_adh, n_acts_partial_adh, n_acts_non_adh, na.rm = T))]
  prop_acts_ring_protected <- prop_acts_full_adh + prop_acts_partial_adh * (rr_ring_partial_adh/rr_ring_full_adh)
  
  # Calculate proportion of HIV infections in the dapivirine arm that are attributable to each source of lack of protection (ring failure to provide protection, non-adherence, and AI)
  n_inf_ring_failure <- f_dt[arm == 1 & hiv == 1 & adh == 1, sum(prop_risk_vi)] +
                        f_dt[arm == 1 & hiv == 1 & adh == 2, sum(prop_risk_vi) * rr_ring_partial_adh/rr_ring_full_adh]
  n_inf_non_adh      <- f_dt[arm == 1 & hiv == 1 & adh == 2, sum(prop_risk_vi) * (1 - rr_ring_partial_adh/rr_ring_full_adh)] +
                        f_dt[arm == 1 & hiv == 1 & adh == 3, sum(prop_risk_vi)]
  n_inf_ai           <- f_dt[arm == 1 & hiv == 1, sum(prop_risk_ai)]
  
  prop_eff_dilution_non_adh <- n_inf_non_adh/sum(n_inf_ring_failure, n_inf_non_adh, n_inf_ai)
  prop_eff_dilution_ai      <- n_inf_ai/sum(n_inf_ring_failure, n_inf_non_adh, n_inf_ai)
  
  # Clean up: Keep only last observation for each participant
  f_dt <- f_dt[, .SD[.N], by = id]
  
  # Calculate the hazard ratio in simulated data
  hr <- unname(summary(f_dt[, coxph(formula = Surv(time = as.numeric(visit_date - enrolldt)/365, event = hiv) ~ arm + strata(site))])$conf.int[1])
  hr_lb <- unname(summary(f_dt[, coxph(formula = Surv(time = as.numeric(visit_date - enrolldt)/365, event = hiv) ~ arm + strata(site))])$conf.int[3])
  hr_ub <- unname(summary(f_dt[, coxph(formula = Surv(time = as.numeric(visit_date - enrolldt)/365, event = hiv) ~ arm + strata(site))])$conf.int[4])
  
  if(debug_rerandomize) {
    # If re-randomization occured, calculate the hazard ratio stratified by site, BV status, and dichotomous average monthly number of sex acts
    hr_strat <- unname(summary(f_dt[, coxph(formula = Surv(time = as.numeric(visit_date - enrolldt)/365, event = hiv) ~ arm + strata(site) + strata(b_bv) + strata(gt_median_acts))])$conf.int[1])
    hr_strat_lb <- unname(summary(f_dt[, coxph(formula = Surv(time = as.numeric(visit_date - enrolldt)/365, event = hiv) ~ arm + strata(site) + strata(b_bv) + strata(gt_median_acts))])$conf.int[3])
    hr_strat_ub <- unname(summary(f_dt[, coxph(formula = Surv(time = as.numeric(visit_date - enrolldt)/365, event = hiv) ~ arm + strata(site) + strata(b_bv) + strata(gt_median_acts))])$conf.int[4])
    
    output <- list(hr = hr, hr_lb = hr_lb, hr_ub = hr_ub, hr_strat = hr_strat, hr_strat_lb = hr_strat_lb, hr_strat_ub = hr_strat_ub, prop_total_acts_ai = prop_total_acts_ai, prev_ai_sim = prev_ai_sim, prop_ai_sim = prop_ai_sim, prop_acts_ring_protected = prop_acts_ring_protected, prop_eff_dilution_non_adh = prop_eff_dilution_non_adh, prop_eff_dilution_ai = prop_eff_dilution_ai, i = i, ri_dt = risk_inf_dt, study_dt = f_dt, n_active_part = n_active_partnerships, n_active_part_dpv = n_active_partnerships_dpv, n_active_part_pbo = n_active_partnerships_pbo, male_dt = m_dt, b_dt = balanced_dt, n_part_hiv_pos_dt = n_part_hiv_dt, hiv_assign_dt = check_hiv_assignment_dt, exposure_dt = exp_dt)
    
  } else {
    
    output <- list(hr = hr, hr_lb = hr_lb, hr_ub = hr_ub, prop_total_acts_ai = prop_total_acts_ai, prev_ai_sim = prev_ai_sim, prop_ai_sim = prop_ai_sim, prop_acts_ring_protected = prop_acts_ring_protected, prop_eff_dilution_non_adh = prop_eff_dilution_non_adh, prop_eff_dilution_ai = prop_eff_dilution_ai, i = i, ri_dt = risk_inf_dt, study_dt = f_dt, n_active_part = n_active_partnerships, n_active_part_dpv = n_active_partnerships_dpv, n_active_part_pbo = n_active_partnerships_pbo, male_dt = m_dt, b_dt = balanced_dt, n_part_hiv_pos_dt = n_part_hiv_dt, hiv_assign_dt = check_hiv_assignment_dt, exposure_dt = exp_dt)
    
  }
  
  rm(list = c("f_dt", "m_dt"))
  
  return(output)
  # save(output, file = paste0(output_dir, "/output_", i, ".RDATA"))
}
