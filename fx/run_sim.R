#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 March 2019
# run_sim: function to run a simulation of the ASPIRE trial
#
# input:  Specified parameters
# output: If calibrating, particle of parameters and likelihood and incidence by age in the placebo arm. Otherwise, quantities of interest in analysis (hazard ratio, proportion of acts protected by adherence among women in active arm, proportion of acts that are AI, proportion efficacy dilution attributable to non-adherence and AI, and risk tracking over simulations.
#
#######################################################################################

run_sim <- function(lambda, cond_rr, c, s, rr_ai, p_rate_rr, base_male_hiv_incidence, prev_ai, prop_ai = 0.078, prop_adh, rr_ring, re_randomize = F, i = NULL, calibrate = F, reduce_output = F) {
  
  load(file = paste0(wd, "/data-private/f_dt.RData"))
  
  if(re_randomize) {
    # Create new variable to indicate if an individual reports a number of average monthly acts that is greater than the median reported number of monthly acts.
    f_dt[, gt_median_acts := as.numeric(pp_vi_acts_avg * b_n_part > median(pp_vi_acts_avg * b_n_part))]
    
    dpv_ids <- re_randomize_arm(dt = f_dt[visit == 0, .(id, site, b_bv, gt_median_acts)])
    
    f_dt[id %in% dpv_ids,    arm := 1]
    f_dt[!(id %in% dpv_ids), arm := 0]
  }
  
  # Impute missing values for condom use in the last week
  f_dt[is.na(b_condom_lweek) & visit == 0, b_condom_lweek := impute_unplw(dt = f_dt[is.na(b_condom_lweek) & visit == 0, .(id, m_hiv, b_sti, f_age, b_noalc, b_married)])]
  f_dt[, b_condom_lweek := na.locf(b_condom_lweek), by = id]
  
  # Assign proportion of AI acts to each woman given study site
  f_dt <- merge(x = f_dt, y = assign_prob_ai(dt = unique(f_dt[, .(site, id)]), prev_ai = prev_ai, prop_ai = prop_ai), by = "id", all = T)
  
  # Redistribute VI acts so that the total number of acts remains approximately constant across simulations with varying prevalence and frequency of AI
  f_dt[, pp_vi_acts_avg := pp_vi_acts_avg * (1 - prob_ai)]
  f_dt[, pp_ai_acts_avg := pp_vi_acts_avg * prob_ai]
  
  # Assign relative probability of having an HIV-positive male partner
  f_dt[, m_hiv_rr := sapply(f_age, function(x) { assign_m_hiv_rr(age = x, c = c, s = s) })]
  
  # Among women who became pregnant, there were no seroconversions. Set m_hiv_rr to 0 for these women to preclude having an HIV-positive partner and prevent incident infections in their partners.
  f_dt[any_preg == 1, m_hiv_rr := 0]
  
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
                     m_hiv_rr       = f_dt[visit == 0]$m_hiv_rr[id],
                     t_inf          = -3) # Assign time of infection as -3 to all initial partners (assumes no partners at baseline have acute infection)
  
  if(nrow(m_dt[max_part == 1]) != nrow(f_dt[visit == 0 & b_married == 1])) { browser() }
  
  # Assign male partner age, HIV status, and viral load values.
  m_dt[, m_age := assign_male_age(m_ids = id, f_age = f_dt[visit == 0]$f_age)]
  m_dt[, hiv   := assign_male_hiv_status(dt = m_dt[, .(country, m_age, b_condom_lweek, m_hiv_rr)], cond_rr = cond_rr)]
  m_dt[, vl    := assign_male_vl_trt_unknown(dt = m_dt[, .(m_age, hiv)])]
  
  # Replace randomly assigned baseline hiv and viral load values for women who report an HIV-positive partner
  m_dt[id %in% f_dt[m_hiv == "positive", id], hiv := 1]
  m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == T, id], vl := log10(50)]
  m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == F, id], 
       vl := assign_male_vl_trt_known(dt = m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == F, id], .(m_age, hiv)])]
  
  # One woman who became pregnant reported an HIV-positive partner and reported that her partner was not on ARVs. As this participant did not seroconvert, assume minimum viral load in her partner.
  m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == F & any_preg == 1, id], vl := log10(50)]
  
  # Among married women, if reported partner HIV status is negative, set partner status to negative
  m_dt[id %in% f_dt[b_married == 1 & m_hiv == "negative" & f_age_cat == "27-45", id], `:=`(hiv = as.integer(0), vl = as.integer(0))]
  
  f_dt[,  `:=`(prob_condom = prop_condom,
               pp_vi_acts  = assign_mvi_acts(dt = f_dt[, .(id, pp_vi_acts_avg)]),
               pp_ai_acts  = assign_mai_acts(dt = f_dt[, .(id, pp_ai_acts_avg)]))]
  
  # Add variable to track HIV infection through simulations
  f_dt[visit == 0, hiv := as.integer(0)]
  
  # Assign adherence to both placebo and active arms. Assign to placebo in order to create a comparable group for analyses comparing women across arms by adherence status. f_dt has key of id and visit before and after assign_adh function, so below code will merge correctly by id and visit.
  f_dt[, adh := assign_adh(dt = f_dt[, .(id, visit, obs_on_study, site, f_age, prob_adh, incr_prob_adh, post_adh_intervention, time_study_gt_three_mos, adh_period)], target_adh = prop_adh)]
  
  # Create table to track cumulative risk in the placebo and active arms. The number of time steps will be variable, so make the table very large, and then reduce size below. 
  risk_inf_dt <- as.data.table(expand.grid(time = 1:60, arm = as.double(0:1), n_inf = NA_real_, n_inf_pre_ring = NA_real_, n_acts = NA_real_, n_inf_vi_excl = NA_real_, n_acts_vi_excl = NA_real_, n_inf_vi_excl_adh = NA_real_, n_acts_vi_excl_adh = NA_real_, n_inf_any_ai = NA_real_, n_acts_any_ai = NA_real_, n_inf_any_ai_adh = NA_real_, n_acts_any_ai_adh = NA_real_))
  setkey(risk_inf_dt, time)
  
  # Create empty data table to track risk of infection at all HIV-exposed sex acts. Push to global environment to be available in hiv_transmission function.
  acts_dt <<- data.table(id = NA_real_, time = NA_real_, arm = NA_real_, risk_inf = NA_real_, hiv_act = NA_real_, act_id = NA_real_, post_enrollment = NA_real_, early_term_censor = NA_real_)
  
  # Create data table to track co-factors for HIV infection by arm and risk, prior to application of RR of ring.
  balanced_dt <- as.data.table(expand.grid(time = 1:60, arm = as.double(0:1), median_age = NA_real_, median_vl = NA_real_, median_prob_condom = NA_real_, mean_sti = NA_real_, prop_bv = NA_real_, n_partners_hiv = NA_real_, n_acts = NA_real_))
  setkey(balanced_dt, time, arm)
  
  # Initialize number of HIV infections to 0 and timestep to 1.
  n_hiv_inf <- 0
  t         <- 1
  
  # Create vector to track the number of active partnerships at each timestep
  n_active_partnerships <- numeric(60)
  
  # Run simulation until 168 infections have occurred. Simulations should run for a minimum of 31 timesteps.
  while(n_hiv_inf < 168 | t <= 31) {
  
    # Calculate the number of HIV-positive partners each woman has
    setkey(m_dt, id) # Remove this when code to track the number of HIV-positive partners each woman has is no longer in use
    f_dt[visit == t - 1, n_part_hiv := m_dt[, sum(hiv == 1 & active == 1), by = id]$V1]
    
    # HIV transmission
    f_dt[visit == t, c("hiv", "n_acts_hiv_pos", "cum_risk_pre_ring", "inf_act_type") := hiv_transmission(
      f_dt    = f_dt[visit == t, .(id, adh, f_age, pp_vi_acts, pp_ai_acts, b_n_sti, b_bv, arm, prob_condom, any_ai = pp_ai_acts_avg > 0, post_enrollment, early_term_censor)],
      m_dt    = m_dt[active == 1 & hiv == 1, .(id, vl, t_inf)],
      t       = t,
      l       = lambda,
      rr_ai   = rr_ai,
      rr_ring = rr_ring,
      ri_dt   = risk_inf_dt,
      b_dt    = balanced_dt)]
    
    # Append acts from timestep t to acts_dt
    acts_dt <<- rbindlist(l = list(acts_dt, per_act_risk_dt))
    
    # Calculate number of acts of each type (VI and AI) and protection (adherent or not adherent)
    f_dt[visit == t, c("n_acts_vi", "n_acts_ai", "n_acts_adh", "n_acts_non_adh") := calculate_acts(f_dt = f_dt[visit == t, .(id, adh, pp_vi_acts, pp_ai_acts)], m_dt = m_dt[active == 1, id])]
    
    # HIV transmission is calculated at each time step independently of HIV status at previous time steps. Track if HIV transmission has occurred at any time step in order to complete incidence and survival analyses correctly. An HIV infection can only be observed while a woman is in the study (i.e., post-enrollment and prior to early termination).
    f_dt[visit == t, any_hiv := as.numeric(f_dt[visit %in% 1:t, any(hiv == 1 & post_enrollment == 1 & early_term_censor == 0), by = id]$V1)]
    
    # Track number of active partnerships at each timestep
    n_active_partnerships[t] <- m_dt[, sum(active)]
    
    # Partner change
    m_dt <- partner_change(m_dt = m_dt, f_dt = f_dt[visit == t, .(id, arm, country, f_age, max_part, b_condom_lweek, m_hiv_rr)], cond_rr = cond_rr, p_rate_rr = p_rate_rr)
    
    # Incident infections among men
    m_dt <- male_hiv_incidence(m_dt = m_dt, base_male_hiv_incidence = base_male_hiv_incidence, cond_rr = cond_rr, t_curr = t)
    
    # Count and update number of HIV infections
    n_hiv_inf <- sum(f_dt[, any(any_hiv == 1, na.rm = T), by = id]$V1)
    
    # Increment timestep
    t <- t + 1
    
    if(t > 60) { break } # If time to achieve 168 infections is greater than 60 months, end simulation.
  }
  
  # Calculate the number of HIV-positive partners by arm at each timestep, prior to censoring
  n_part_hiv_dt <- f_dt[!is.na(n_part_hiv), .(n_part_hiv = sum(n_part_hiv), pre_censor = 1), by = c("visit", "arm")]
  
  # Create adherence data table to estimate the proportion of women adherent before and after adherence interventions
  #adh_dt <- unique(f_dt[, .(id, post_adh_intervention, time_study_gt_three_mos, site, f_age_cat, adh)])
  adh_dt <- unique(f_dt[, .(id, arm, visit, post_adh_intervention, time_study_gt_three_mos, site, f_age_cat, adh)])
  
  # Clean up: Remove rows prior to a participant enrolling and after early termination.
  f_dt <- f_dt[post_enrollment == 1 & early_term_censor == 0]
  
  # Clean up: Among seroconverters, remove rows subsequent to visit at which HIV is first detected (hiv_transmission module takes only the values for the current time step, and so does not take account of previous HIV status values).
  f_dt[visit == 0, any_hiv := 0]
  ids_hiv <- f_dt[any_hiv == 1, unique(id)]
  f_dt <- f_dt[, .SD[cumsum(any_hiv) <= 1], by = id]
  
  if(any(f_dt[, cumsum(hiv), by = id]$V1 > 1)) { stop("HIV censoring not done correctly.") }
  if(!all(f_dt[hiv == 1, unique(id)] == ids_hiv)) { stop("HIV censoring not done correctly.") }
  
  # Simulations will run until 168 infections have occurred, and censoring code above will remove all visits subsequent to the last-occurring HIV infection. Identify the max visit number and remove rows/indices subsequent to max visit number from risk_inf_dt, balanced_dt, and n_active_partnerships.
  censor_visit <- f_dt[, max(visit)]
  
  risk_inf_dt           <- risk_inf_dt[time <= censor_visit]
  balanced_dt           <- balanced_dt[time <= censor_visit]
  n_active_partnerships <- n_active_partnerships[1:censor_visit]
  
  # Reset visit number to begin from 1 for each participant.
  f_dt[, visit := 1:.N, by = id]
  
  # Keep longitudinal dataset to estimate HR over time and to evaluate adherence
  long_dt <- f_dt[, .(id, visit, hiv, site, f_age_cat, arm, adh, prob_adh, post_adh_intervention, time_study_gt_three_mos, ring_worries, pknow_study, same_psp, bl_ocp)]
  
  exp_dt <- f_dt[, .(exposures = sum(n_acts_hiv_pos, na.rm = T)), by = c("visit", "arm")]
  
  # From acts_dt, calculate a woman's cumulative risk at each time step, up until her first HIV infection
  acts_dt <- acts_dt[post_enrollment == 1 & early_term_censor == 0]
  acts_dt[, min_time := min(time), by = id]
  acts_dt[, visit := time - min_time]
  acts_dt <- acts_dt[, .SD[cumsum(cumsum(hiv_act)) <= 1], by = id]
  acts_dt[, max_visit := max(visit), by = id]
  cum_risk_dt <- rbindlist(l = lapply(0:(acts_dt[, max(visit)]), function(x) acts_dt[visit %in% 0:x, .(visit = x, arm = unique(arm), max_visit = unique(max_visit), cum_risk = 1 - Reduce(f = `*`, x = (1 - risk_inf))), by = id]))
  setkey(cum_risk_dt, id, visit)
  cum_risk_dt <- cum_risk_dt[visit <= max_visit]
  cum_risk_dt <- as.data.table(cum_risk_dt %>% group_by(arm, visit) %>% summarise(mean_cum_risk = mean(cum_risk)))
  
  # Calculate the number of HIV-positive partners by arm at each timestep, after censoring
  n_part_hiv_dt <- rbindlist(l = list(n_part_hiv_dt, f_dt[!is.na(n_part_hiv), .(n_part_hiv = sum(n_part_hiv), pre_censor = 0), by = c("visit", "arm")]))
  
  # Calculate proportion of acts of each type
  prop_total_acts_ai <- f_dt[, sum(n_acts_ai, na.rm = T)/sum(n_acts_ai, n_acts_vi, na.rm = T)]
  prev_ai_sim <- sum(f_dt[, any(n_acts_ai > 0, na.rm = T), by = id]$V1)/f_dt[, length(unique(id))]
  prop_ai_sim <- f_dt[id %in% f_dt[, any(n_acts_ai > 0, na.rm = T), by = id][V1 == T, id], sum(n_acts_ai, na.rm = T)/(sum(n_acts_ai, n_acts_vi, na.rm = T))]
  
  # Calculate proportion of ring-protected acts among women in dapivirine arm
  prop_acts_adh <- f_dt[arm == 1, sum(n_acts_adh, na.rm = T)/(sum(n_acts_adh, n_acts_non_adh, na.rm = T))]
  
  # Calculate proportion of HIV infections in the dapivirine arm that are attributable to non-adherence and AI
  n_inf_non_adh <- f_dt[arm == 1 & adh == 0 & inf_act_type == "vi", sum(hiv)]
  n_inf_ai      <- f_dt[arm == 1 & inf_act_type == "ai", sum(hiv)]
  
  prop_eff_dilution_non_adh <- n_inf_non_adh/(n_inf_non_adh + n_inf_ai)
  prop_eff_dilution_ai      <- n_inf_ai/(n_inf_non_adh + n_inf_ai)
  
  # Clean up: Keep only last observation for each participant
  f_dt <- f_dt[, .SD[.N], by = id]
  
  if(calibrate) {
    # Sum infections, sample size, and person-time by age group for the placebo arm (ABC will fit parameters only to placebo arm).
    inf_sim_age <- as.data.table(f_dt[arm == 0] %>% group_by(f_age_cat) %>% summarise(hiv_inf = sum(hiv), n = length(unique(id)), py = sum((visit * 30)/365)))
    
    # Calculate incidence in placebo arm for months 15-23 and 24-31
    inf_sim_time <- data.table(f_age_cat = c("18-45", "18-45"),
                               hiv_inf   = c(long_dt[arm == 0 & visit %in% 15:23, sum(hiv)],
                                             long_dt[arm == 0 & visit %in% 24:31, sum(hiv)]),
                               n         = c(long_dt[arm == 0 & visit %in% 15:23, length(unique(id))],
                                             long_dt[arm == 0 & visit %in% 24:31, length(unique(id))]),
                               py        = c(sum(((long_dt[arm == 0 & visit %in% 15:23, max(visit), by = id]$V1 - 14) * 30)/365),
                                             sum(((long_dt[arm == 0 & visit %in% 24:31, max(visit), by = id]$V1 - 23) * 30)/365)))
    
    inf_sim <- rbindlist(l = list(inf_sim_age, inf_sim_time))
    
    # Calculate the distance criterion as the probability of the observed data under the poisson distribution paramaterized with lambda from simulated data (i.e., L(parameters | data)). Multiplying likelihood of all age categories assumes independent poisson distributions.
    rho <- prod(sapply(1:nrow(inf_sim), function(x) { dpois(x = inf_obs[x, hiv_inf], lambda = inf_sim[x, hiv_inf/py] * inf_obs[x, py]) }))
    
    particle <- list(lambda = lambda, cond_rr = cond_rr, c = c, s = s, rr_ai = rr_ai, p_rate_rr = p_rate_rr, base_male_hiv_incidence = base_male_hiv_incidence, rho = rho, i = i)
    
    save(particle, file = paste0(getwd(), "/particle_", particle[["i"]], ".RData"))
    save(inf_sim,  file = paste0(getwd(), "/inf_sim_", i, ".RData"))
    
    rm(list = c("f_dt", "m_dt", "acts_dt", "balanced_dt", "cum_risk_dt", "exp_dt", "n_part_hiv_dt", "risk_inf_dt"))
    rm(acts_dt, envir = globalenv())
  } else {
    # Sum infections, sample size, and person-time by age group for the placebo arm (ABC will fit parameters only to placebo arm).
    inf_sim_age <- as.data.table(f_dt[arm == 0] %>% group_by(f_age_cat) %>% summarise(hiv_inf = sum(hiv), n = length(unique(id)), py = sum((visit * 30)/365)))
    
    # Calculate incidence in placebo arm for months 15-23 and 24-31
    inf_sim_time <- data.table(f_age_cat = c("Months 15-23", "Months 24-31"),
                               hiv_inf   = c(long_dt[arm == 0 & visit %in% 15:23, sum(hiv)],
                                             long_dt[arm == 0 & visit %in% 24:31, sum(hiv)]),
                               n         = c(long_dt[arm == 0 & visit %in% 15:23, length(unique(id))],
                                             long_dt[arm == 0 & visit %in% 24:31, length(unique(id))]),
                               py        = c(sum(((long_dt[arm == 0 & visit %in% 15:23, max(visit), by = id]$V1 - 14) * 30)/365),
                                             sum(((long_dt[arm == 0 & visit %in% 24:31, max(visit), by = id]$V1 - 23) * 30)/365)))
    
    inf_sim <- rbindlist(l = list(inf_sim_age, inf_sim_time))
    
    # Calculate the hazard ratio in simulated data
    hr <- unname(summary(f_dt[, coxph(formula = Surv(time = (visit * 30)/365, event = hiv) ~ arm + strata(site))])$conf.int[1])
    hr_lb <- unname(summary(f_dt[, coxph(formula = Surv(time = (visit * 30)/365, event = hiv) ~ arm + strata(site))])$conf.int[3])
    hr_ub <- unname(summary(f_dt[, coxph(formula = Surv(time = (visit * 30)/365, event = hiv) ~ arm + strata(site))])$conf.int[4])
    
    if(reduce_output) {
      output <- list(hr = hr, prop_eff_dilution_non_adh = prop_eff_dilution_non_adh, prop_eff_dilution_ai = prop_eff_dilution_ai, ri_dt = risk_inf_dt, study_dt = f_dt, exposure_dt = exp_dt, inc_sim = inf_sim)
    } else {
      output <- list(hr = hr, hr_lb = hr_lb, hr_ub = hr_ub, prop_total_acts_ai = prop_total_acts_ai, prev_ai_sim = prev_ai_sim, prop_ai_sim = prop_ai_sim, prop_eff_dilution_non_adh = prop_eff_dilution_non_adh, prop_eff_dilution_ai = prop_eff_dilution_ai, ri_dt = risk_inf_dt, study_dt = f_dt, n_active_part = n_active_partnerships, male_dt = m_dt, b_dt = balanced_dt, n_part_hiv_pos_dt = n_part_hiv_dt, exposure_dt = exp_dt, cumulative_risk_dt = cum_risk_dt, risk_per_act_dt = acts_dt, adh_dt = adh_dt, long_dt = long_dt, inc_sim = inf_sim)
    }
    
    rm(list = c("f_dt", "m_dt", "acts_dt", "balanced_dt", "cum_risk_dt", "exp_dt", "n_part_hiv_dt", "risk_inf_dt", "adh_dt", "long_dt"))
    rm(acts_dt, envir = globalenv())
    
    return(output)
  }
}
