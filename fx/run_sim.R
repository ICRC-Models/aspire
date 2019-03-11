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

run_sim <- function(lambda = params$lambda, cond_rr = params$cond_rr, c = params$c, s = params$s, rr_ai = params$rr_ai) {
  
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
  m_dt[id %in% f_dt[b_married == 1 & m_hiv == "negative" & b_f_age_cat == "27-45", id], `:=`(hiv = as.integer(0), vl = as.integer(0))]
  
  # Append fixed characteristics to f_dt that will vary in simulations
  f_dt[,  `:=`(prob_condom = assign_prob_condom(countries = f_dt[, country]),
               pp_vi_acts  = assign_mvi_acts(dt = f_dt[, .(id, pp_vi_acts_avg)]),
               pp_ai_acts  = assign_mai_acts(dt = f_dt[, .(id, pp_ai_acts_avg)]))]
  
  # Add variables to track through simulations - hiv and proportion of risk attributable to AI
  f_dt[visit == 0, `:=`(hiv          = as.integer(0),
                        prop_risk_ai = 0)]
  
  # Assign adherence in quarterly intervals
  for(t in seq(3, 33, 3)) {
    if(t == 3) {
      f_dt[visit %in% c(t, t - 1, t - 2), adh := assign_adh_q1(dt = f_dt[visit %in% c(t, t - 1, t - 2), .(id, visit, b_f_age, b_married, b_edu, b_gon, b_ct, b_trr, b_pknow, b_circ, b_noalc, b_npart_eq_two, b_npart_gt_two, enrolldt_after_Apr12013, arm, site)])]
    } else {
      f_dt[visit %in% c(t, t - 1, t - 2), adh := assign_adh_fup(dt = f_dt[visit %in% c(t, t - 1, t - 2, t - 3), .(id, visit, adh, b_f_age, b_married, b_edu, b_gon, b_ct, b_trr, b_pknow, b_circ, b_noalc, b_npart_eq_two, b_npart_gt_two, enrolldt_after_Apr12013, prev_visit_date_after_Aug12013, arm, site)], t = t)]
    }
  }
  
  # Run simulation
  for(t in 1:33) { # No participants are on study beyond visit 33
    
    f_dt[visit == t, c("hiv", "prop_risk_ai") := hiv_transmission(f_dt  = f_dt[visit == t, .(id, adh, f_age, pp_vi_acts, pp_ai_acts, b_n_sti, b_bv, arm, prob_condom)],
                                                                  m_dt  = m_dt[active == 1 & hiv == 1, .(id, vl)],
                                                                  t     = t,
                                                                  l     = lambda,
                                                                  rr_ai = rr_ai)]
    
    f_dt[visit == t, any_hiv := as.numeric(f_dt[visit %in% 1:t, any(hiv == 1), by = id]$V1)]
    
    m_dt <- partner_change(m_dt = m_dt, f_dt = f_dt[, .(id, country, f_age, max_part, b_condom_lweek, m_hiv_rr)], cond_rr = cond_rr)
  }
  
  # Clean up: Remove all rows where on_study == 0
  f_dt <- f_dt[on_study == 1]
  
  # Clean up: Among seroconverters, remove rows subsequent to visit at which HIV is first detected (hiv_transmission module takes only the values for the current time step, and so does not take account of previous HIV status values).
  f_dt[visit == 0, any_hiv := 0]
  f_dt <- f_dt[, .SD[cumsum(any_hiv) <= 1], by = id]
  
  return(f_dt)
}
