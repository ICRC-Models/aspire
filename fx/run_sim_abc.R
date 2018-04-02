#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   20 March 2018
# run_sim_abc: function to run a single simulation with AI parameters drawn in iterative Approximate Bayesian Computation model-fitting procedure.
#
# input:  Parameters for proportion of women engaging in AI and alpha and beta parameters determining frequency of AI among those who engage in AI.
# output: AI parameters used in simulation and distance criterion rho
#
#######################################################################################

run_sim_abc <- function(prop_ai, a, b, i) {
  
  setkey(f_dt, id)
  
  ## Create vector of participant ids repeated once for each partner
  id <- f_dt[, rep(id, n_part)]
  
  ## Create table with one row for each male partner
  m_dt <- data.table(id           = id,
                     country      = f_dt$country[id],
                     f_age        = f_dt$f_age[id],
                     active       = as.integer(1),
                     max_part     = f_dt$max_part[id],
                     condom_lweek = f_dt$condom_lweek[id],
                     m_hiv_rr     = f_dt$m_hiv_rr[id])
  
  if(nrow(m_dt[max_part == 1]) != nrow(f_dt[married == T])) { browser() }
  
  ## Assign male partner age, HIV status, and viral load values.
  m_dt[, m_age := assign_male_age(m_ids = id, f_age = f_age)]
  m_dt[, hiv   := assign_male_hiv_status(dt = m_dt[, .(country, m_age, f_age, condom_lweek, m_hiv_rr)], cond_rr = params$cond_rr)]
  m_dt[, vl    := assign_male_vl(hiv = m_dt$hiv, age = m_dt$m_age)]
  
  ## Replace randomly assigned baseline hiv and viral load values for women who report an HIV-positive partner
  m_dt[id %in% f_dt[m_hiv == "positive", id], hiv := 1]
  m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == T, id], vl := log10(50)]
  m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == F, id], 
       vl := assign_male_vl(hiv = m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == F, id]]$hiv,
                            age = m_dt[id %in% f_dt[m_hiv == "positive" & m_arv == F, id]]$m_age)] ## TO DO: this distribution includes treated males, but these males are not on treatment. Should sample only from distribution that excludes viral suppresion.
  
  ## Among married women, if reported partner HIV status as negative, set partner status to negative
  m_dt[id %in% f_dt[married == T & m_hiv == "negative", id], `:=`(hiv = as.integer(0), vl = as.integer(0))]
  
  ## Append fixed characteristics to f_dt that will vary in simulations
  f_dt[, `:=`(prob_ai     = assign_prob_ai(n = nrow(f_dt), prop_ai = prop_ai, alpha = a, beta = b),
              prob_condom = assign_prob_condom(countries = f_dt[, country]))]
  
  ## Create study data table
  time_steps <- seq(0, 30 * 12 * 3, 30)
  study_dt <- setorder(as.data.table(expand.grid(f_dt$id, time_steps)))
  setnames(study_dt, old = c("Var1", "Var2"), new = c("id", "time"))
  
  ## Add baseline variable values
  study_dt[time == 0, `:=`(f_age = f_dt$f_age,
                           adh   = as.integer(assign_adh_t0(dt = f_dt[, .(id, f_age, f_age_cat, arm, site, days_pre_adh_int)])),
                           hiv   = as.integer(0))]
  
  setkey(study_dt, id, time)
  
  study_dt <- merge(x = study_dt, y = f_dt[, .SD, .SDcols = names(f_dt)[which(names(f_dt) != "f_age")]])
  
  ## Assign study retention status according to observed number of days of follow-up
  ## TO DO: Ask EB: values of fu_days do not align with study visit schedule. How was the exact number of days of follow-up determined when censoring occurred between study visits?
  study_dt[, on_study := ifelse(time <= fu_days, 1, 0)]
  
  ## Specify days that are pre-adherence intervention
  study_dt[, pre_adh_int := ifelse(time < days_pre_adh_int, 1, 0)]
  
  ## Run simulation
  for(t in time_steps[2:37]) {
    
    study_dt[time == t, `:=`(f_age = aging(prev_ages = study_dt[time == t - 30, f_age]),
                             adh   = assign_adh_fup(s_dt = study_dt[time %in% c(t, t - 30), .(id, time, arm, adh, pre_adh_int)], f_dt = f_dt[, .(id, f_age_cat, days_pre_adh_int)], t = t),
                             hiv   = hiv_transmission(s_dt = study_dt[time == t - 30, .(id, adh, hiv)],
                                                      m_dt = m_dt[active == 1 & hiv == 1, .(id, vl)],
                                                      f_dt = f_dt[, .(id, f_age, pp_acts, n_sti, bv, arm, prob_ai, prob_condom)],
                                                      t   = t,
                                                      l   = params$lambda))]
    
    m_dt <- partner_change(m_dt = m_dt, f_dt = f_dt[, .(id, country, f_age, max_part, condom_lweek, m_hiv_rr)])
  }
  
  ## Clean up: Remove all rows where on_study == 0
  study_dt <- study_dt[on_study == 1]
  
  ## Clean up: Among seroconverters, remove rows subsequent to visit at which HIV is first detected
  study_dt <- study_dt[, .SD[cumsum(hiv) <= 1], by = id]  
  
  ## Clean up: Keep only last observation for each participant
  study_dt <- study_dt[, .SD[.N], by = id]
  
  ## TO DO: Ask EB: Sum of simulated infections by age group uses age group at baseline, as does summary of observed infections. Is this correct?
  ## Sum infections, sample size, and person-time by arm and age group
  inf_sim <- as.data.table(study_dt %>% group_by(arm, f_age_cat) %>% summarise(hiv_inf = sum(hiv), n = length(unique(id)), py = sum(time/365)))
   
  if(!all(inf_sim[, n] == inf_obs[, n])) { stop("Number of participants in simulation not equal to number of participants in trial.") }
  
  ## Calculate the distance criterion as the probability of the observed data under the poisson distribution paramaterized with lambda from simulated data (i.e., L(parameters | data)). Multiplying likelihood of all arm and age categories assumes independent poisson distributions. TO DO: Review this with EB.
  rho <- prod(sapply(1:nrow(inf_sim), function(x) { dpois(x = inf_obs[x, hiv_inf], lambda = inf_sim[x, hiv_inf/py] * inf_obs[x, py])}))
  
  particle <- list(prop_ai = prop_ai, a = a, b = b, rho = rho, i = i)
  
  save(particle, file = paste0(getwd(), "/particle_", particle[["i"]], ".RDATA"))
  
  ## TO DO: Include code that tracks # ai acts, # vi acts at each time step. Plot as distribution with expected proportions for code testing.
}
