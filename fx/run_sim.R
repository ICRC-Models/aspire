#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 March 2018
# run_sim: function to run a simulation of the ASPIRE trial
#
# input:  None
# output: A longitudinal dataset
#
#######################################################################################

run_sim <- function(lambda = params$lambda, cond_rr = params$cond_rr, c = params$c, s = params$s, rr_ai = params$rr_ai, rr_ring_full_adh = params$rr_ring_full_adh, rr_ring_partial_adh = params$rr_ring_partial_adh, debug_sim_hiv_acts = F) {
  
  load(file = paste0(wd, "/data-private/f_dt.RData"))
  
  # Impute missing values for condom use in the last week and number of anal sex acts
  f_dt[is.na(b_condom_lweek) & visit == 0, b_condom_lweek := impute_unplw(dt = f_dt[is.na(b_condom_lweek) & visit == 0, .(id, m_hiv, b_sti, f_age, b_noalc, b_married)])]
  f_dt[, b_condom_lweek := na.locf(b_condom_lweek), by = id]
  
  f_dt[is.na(anal_n_acts) & (visit == 0 | visit == 3), anal_n_acts := impute_anal_n_acts(dt = f_dt[is.na(anal_n_acts) & (visit == 0 | visit == 3), .(id, visit, site, b_edu, b_depo, b_neten, b_trans_sex)])]
  
  # Estimate per-partner average number of monthly anal acts (total average number of monthly AI acts divided by number of partners).
  f_dt[, n_ai_monthly_avg := round(mean(anal_n_acts, na.rm = T)/3, 2), by = id]
  f_dt[, pp_ai_acts_avg := round(n_ai_monthly_avg/b_n_part, 2)]
  
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
  m_dt[, m_age := assign_male_age(m_ids = id, f_age = f_age)]
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
  f_dt[,  `:=`(prob_condom = prop_condom,
               pp_vi_acts  = pp_vi_acts_avg,
               pp_ai_acts  = pp_ai_acts_avg)]
  
  # f_dt[,  `:=`(prob_condom = assign_prob_condom(countries = f_dt[, country]),
  #              pp_vi_acts  = assign_mvi_acts(dt = f_dt[, .(id, pp_vi_acts_avg)]),
  #              pp_ai_acts  = assign_mai_acts(dt = f_dt[, .(id, pp_ai_acts_avg)]))]
  
  # Add variables to track through simulations - hiv and proportion of risk attributable to AI
  f_dt[visit == 0, `:=`(hiv          = as.integer(0),
                        prop_risk_vi = 0,
                        prop_risk_ai = 0)]
  
  # Assign adherence in quarterly intervals
  for(t in seq(3, 150, 3)) {
    if(t == 3) {
      f_dt[visit %in% c(t, t - 1, t - 2), adh := assign_adh_q1(dt = f_dt[visit %in% c(t, t - 1, t - 2), .(id, visit, f_age, b_married, b_edu, b_gon, b_ct, b_trr, b_pknow, b_circ, b_noalc, b_npart_eq_two, b_npart_gt_two, enrolldt_after_Apr12013, arm, site)])]
    } else {
      f_dt[visit %in% c(t, t - 1, t - 2), adh := assign_adh_fup(dt = f_dt[visit %in% c(t, t - 1, t - 2, t - 3), .(id, visit, adh, f_age, b_married, b_edu, b_gon, b_ct, b_trr, b_pknow, b_circ, b_noalc, b_npart_eq_two, b_npart_gt_two, enrolldt_after_Apr12013, prev_visit_date_after_Aug12013, arm, site)], t = t)]
    }
  }
  
  # Initialize number of HIV infections to 0 and timestep to 1.
  n_hiv_inf <- 0
  t         <- 1
  
  # Run simulation until 168 infections have occurred
  while(n_hiv_inf < 168) {
    
    # HIV transmission
    f_dt[visit == t, c("hiv", "prop_risk_vi", "prop_risk_ai") := hiv_transmission(
      f_dt  = f_dt[visit == t, .(id, adh, f_age, pp_vi_acts, pp_ai_acts, b_n_sti, b_bv, arm, prob_condom, prop_condom)],
      m_dt  = m_dt[active == 1 & hiv == 1, .(id, vl)],
      t     = t,
      l     = lambda,
      rr_ai = rr_ai,
      rr_ring_full_adh   = rr_ring_full_adh,
      rr_ring_partial_adh = rr_ring_partial_adh)]
    
    # Calculate number of acts of each type (VI and AI) and protection (adherent, partially adherent, not adherent)
    f_dt[visit == t, c("n_acts_vi", "n_acts_ai", "n_acts_full_adh", "n_acts_partial_adh", "n_acts_non_adh") := calculate_acts(f_dt = f_dt[visit == t, .(id, adh, pp_vi_acts, pp_ai_acts)], m_dt = m_dt[active == 1, id])]
    
    # HIV transmission is calculated at each time step independently of HIV status at previous time steps. Track if HIV transmission has occurred at any time step in order to complete incidence and survival analyses correctly.
    f_dt[visit == t, any_hiv := as.numeric(f_dt[visit %in% 1:t, any(hiv == 1), by = id]$V1)]
    if(!(all(f_dt[visit %in% t, .(id, visit, hiv)]$hiv[f_dt[visit %in% t, any(hiv == 1), by = id]$V1]) == 1)) { browser() }
    
    # Partner change
    m_dt <- partner_change(m_dt = m_dt, f_dt = f_dt[visit == t, .(id, arm, country, f_age, max_part, b_condom_lweek, m_hiv_rr)], cond_rr = cond_rr, debug_sim_hiv_acts = debug_sim_hiv_acts)
    
    # Count and update number of HIV infections
    n_hiv_inf <- sum(f_dt[, any(any_hiv == 1, na.rm = T), by = id]$V1)
    
    # Increment timestep
    t <- t + 1
    
    if(t > 150) { browser() }
  }
  
  # Clean up: Remove all rows where on_study == 0
  # f_dt <- f_dt[on_study == 1]
  
  # Clean up: Among seroconverters, remove rows subsequent to visit at which HIV is first detected (hiv_transmission module takes only the values for the current time step, and so does not take account of previous HIV status values).
  f_dt[visit == 0, any_hiv := 0]
  ids_hiv <- f_dt[any_hiv == 1, unique(id)]
  f_dt <- f_dt[, .SD[cumsum(any_hiv) <= 1], by = id]
  
  if(any(f_dt[, cumsum(hiv), by = id]$V1 > 1)) { stop("HIV censoring not done correctly.") }
 # if(!all(f_dt[hiv == 1, unique(id)] == ids_hiv)) { stop("HIV censoring not done correctly.") }
  if(!all(f_dt[hiv == 1, unique(id)] == ids_hiv)) { browser() }
  
  # Simulations should run until 168 infections have occurred. Identify visit with number of cumulative infections that is closest to 168, and remove visits after that visit.
  inf_by_visit <- f_dt[, .(n_inf = sum(hiv)), by = visit]
  inf_by_visit[, cum_inf := cumsum(n_inf)]
  censor_visit <- inf_by_visit[which.min(abs(cum_inf - 168)), visit]
  f_dt <- f_dt[visit <= censor_visit]
  
  output <- list(study_dt = f_dt)
  
  rm(list = c("f_dt", "m_dt"))
  
  return(output)
}
