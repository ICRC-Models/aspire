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

run_sim_secondary <- function(lambda = params$lambda, cond_rr = params$cond_rr, c = params$c, s = params$s, rr_ai, prev_ai, prop_ai, prop_full_adh, prop_partial_adh, prop_non_adh, rr_ring_full_adh, rr_ring_partial_adh, i) {
  
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
  
  # Add variables to track through simulations - hiv and proportion of risk attributable to VI/AI
  f_dt[visit == 0, `:=`(hiv          = as.integer(0),
                        prop_risk_vi = 0,
                        prop_risk_ai = 0)]
  
  # Assign adherence
  f_dt[arm == 1, adh := assign_adh_secondary(dt = f_dt[arm == 1, .(id, visit)], prop_full_adh = prop_full_adh, prop_partial_adh = prop_partial_adh, prop_non_adh = prop_non_adh)]
  
  # Run simulation
  for(t in 1:33) { # No participants are on study beyond visit 33
    
    # HIV transmission
    f_dt[visit == t, c("hiv", "prop_risk_vi", "prop_risk_ai") := hiv_transmission(
      f_dt  = f_dt[visit == t, .(id, adh, f_age, pp_vi_acts, pp_ai_acts, b_n_sti, b_bv, arm, prob_condom)],
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
    
    # Partner change
    m_dt <- partner_change(m_dt = m_dt, f_dt = f_dt[, .(id, country, f_age, max_part, b_condom_lweek, m_hiv_rr)], cond_rr = cond_rr)
  }
  
  # Clean up: Remove all rows where on_study == 0
  f_dt <- f_dt[on_study == 1]
  
  # Clean up: Among seroconverters, remove rows subsequent to visit at which HIV is first detected (hiv_transmission module takes only the values for the current time step, and so does not take account of previous HIV status values).
  f_dt[visit == 0, any_hiv := 0]
  f_dt <- f_dt[, .SD[cumsum(any_hiv) <= 1], by = id]
  
  if(any(f_dt[, cumsum(hiv), by = id]$V1 > 1)) { stop("HIV censoring not done correctly.") }
  
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
  
  # Sum infections, sample size, and person-time by arm and age group
  inf_sim <- as.data.table(f_dt %>% group_by(arm, b_f_age_cat) %>% summarise(hiv_inf = sum(hiv), n = length(unique(id)), py = sum(as.numeric(visit_date - enrolldt)/365)))
  
  # Calculate the likelihood as the probability of the observed data under the poisson distribution paramaterized with lambda from simulated data (i.e., L(parameters | data)). Multiplying likelihood of all arm and age categories assumes independent poisson distributions. TO DO: Review this with EB.
  likelihood <- prod(sapply(1:nrow(inf_sim), function(x) { dpois(x = inf_obs_arm[x, hiv_inf], lambda = inf_sim[x, hiv_inf/py] * inf_obs_arm[x, py])}))
  
  output <- list(likelihood = likelihood, prop_total_acts_ai = prop_total_acts_ai, prev_ai_sim = prev_ai_sim, prop_ai_sim = prop_ai_sim, prop_acts_ring_protected = prop_acts_ring_protected, prop_eff_dilution_non_adh = prop_eff_dilution_non_adh, prop_eff_dilution_ai = prop_eff_dilution_ai, i = i)
  
  save(output, file = paste0(getwd(), "/output_", i, ".RDATA"))
}
