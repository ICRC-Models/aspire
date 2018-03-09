#######################################################################################
# 
# Kathryn Peebles
# 2017-12-13
# run_sim: function to run a specified number of simulations of the ASPIRE trial
#
# input:  None
# output: A longitudinal dataset
#
#######################################################################################

run_sim <- function() {

  ## Identify females with multiple partners
  ids_mp <- unname(unlist(sapply(c("mal", "sa", "uga", "zim"), function(x) {
    f_dt[country == x, sample(x = id, size = part_dt[country == x, two_p], replace = F)]
  })))
  ids_p  <- sort(c(f_dt[, id], ids_mp))
  
  ## Assign male partner age. Assumes male-female partner age distribution is the same by country.
  f_ages_cat <- colnames(age_mix_mat_cond)[findInterval(x = f_dt$age[ids_p], vec = c(15, 20, 25, 30, 35, 40, 45, 50))] # vec creates bins [15, 20), [20, 25), [25, 30), [30, 35), [35, 40), [40, 45), [45, 50)
  m_partner_age <- unname(sapply(f_ages_cat, function(age_cat) {
    sample(x = rownames(age_mix_mat_cond), size = 1, replace = F, prob = age_mix_mat_cond[, age_cat])
  }))
  
  # ## Check distribution of female ages by male partner age category
  # boxplot(formula = age_country_dt$age[ids_p] ~ m_partner_age, ylab = "female age", xlab = "male partner age category")
  
  ## Given male partner age, assign male partner HIV status as a binomial draw with probability equal to age- and country-specific prevalence
  m_dt <- data.table(ids_p = ids_p, country = f_dt$country[ids_p], age = m_partner_age, hiv = as.integer(0), vl = 0)
  sapply(c("mal", "sa", "uga", "zim"), function(x) {
    m_dt[country == x,
         hiv := rbinom(n = nrow(m_dt[country == x]), size = 1, prob = male_prev[age, x])]
  })
  
  ## Assign viral load values to HIV+ males
  m_dt[, vl := ifelse(hiv == 0, 0,
                      log10(unlist(sapply(sample(1:5,
                                                 size = nrow(m_dt[hiv == 1]),
                                                 replace = T,
                                                 prob = vl_dist[, "proportion"]),
                                          function(x) { runif(n = 1,
                                                              min = vl_dist[x, "min_vl"],
                                                              max = vl_dist[x, "max_vl"]) }))))]
  
  # ## Check distribution of viral load values
  # ggplot(data = m_dt[hiv == 1], aes(x = vl, y = 1:nrow(m_dt[hiv == 1]))) +
  #   geom_point() + labs(x = "Male viral load (log10)", y = "")
  # 
  # ## Check that mean VL is approximately equal across countries
  # m_dt %>% filter(hiv == 1) %>% group_by(country) %>% summarise(mean = mean(vl))
  
  vl_dt <- as.data.table(m_dt %>% select(ids_p, vl) %>% group_by(ids_p) %>% mutate(i = 1:n()))
  vls_p1 <- vl_dt[i == 1, vl]
  vls_p2 <- sapply(unique(vl_dt$id), function(x) { ifelse(x %in% ids_mp, vl_dt[ids_p == x & i == 2]$vl, NA) })
  
  ## Append additional fixed characteristics to f_dt
  f_dt[, `:=`(prob_ai     = assign_prob_ai(n = nrow(f_dt), prop_ai = params$prop_ai, alpha = params$alpha, beta = params$beta),
              prob_condom = assign_prob_condom(countries = f_dt[, country]),
              coital_freq = assign_avg_acts(countries = f_dt[, country]),
              censor_time = runif(n = nrow(f_dt), min = 30 * 12, max = 30 * 12 * 3))]
  
  ## TO DO: Include STIs. Talk to EB about which ones to include.
  
  ## Randomly select women from active arm of Malawi and South Africa to be members of sites with significantly different adherence, per email dated 2017-10-17 from Jingyang.
  blantyre    <- sample(f_dt[country == "mal" & arm == 1]$id, size = 64)
  lilongwe    <- f_dt[country == "mal" & arm == 1 & !(id %in% blantyre)]$id
  umkomaas    <- sample(f_dt[country == "sa" & arm == 1]$id, size = 51)
  isipingo    <- sample(f_dt[country == "sa" & arm == 1 & !(id %in% umkomaas)]$id, size = 58)
  tongaat     <- sample(f_dt[country == "sa" & arm == 1 & !(id %in% c(umkomaas, isipingo))]$id, size = 52)
  emavundleni <- sample(f_dt[country == "sa" & arm == 1 & !(id %in% c(umkomaas, isipingo, tongaat))]$id, size = 80)
  
  f_dt[, `:=`(blantyre    = ifelse(id %in% blantyre, 1, 0),
              lilongwe    = ifelse(id %in% lilongwe, 1, 0),
              umkomaas    = ifelse(id %in% umkomaas, 1, 0),
              isipingo    = ifelse(id %in% isipingo, 1, 0),
              tongaat     = ifelse(id %in% tongaat, 1, 0),
              emavundleni = ifelse(id %in% emavundleni, 1, 0))]
  
  setkey(f_dt, id)
  
  ## Create study data table
  time_steps <- seq(0, 30 * 12 * 3, 30)
  study_dt <- setorder(as.data.table(expand.grid(f_dt$id, time_steps)))
  setnames(study_dt, old = c("Var1", "Var2"), new = c("id", "time"))
  
  ## Add baseline variable values
  study_dt[time == 0, `:=`(age      = f_dt$age,
                           preg     = as.integer(0),
                           vl_p1    = vls_p1,
                           vl_p2    = vls_p2,
                           adh      = as.integer(assign_adh_t0(dt = f_dt[, .(id, age, arm, blantyre, lilongwe, umkomaas, isipingo, tongaat, emavundleni)])),
                           hiv      = as.integer(0),
                           on_study = as.integer(1))]
  
  setkey(study_dt, id, time)
  
  study_dt <- merge(x = study_dt, y = f_dt[, .SD, .SDcols = names(f_dt)[which(names(f_dt) != "age")]])
  
  ## Assign number of monthly acts of each type (ai and vi) and protection status (condom and no condom)
  study_dt[, n_acts := coital_freq + rnorm(n = nrow(study_dt), mean = 0, sd = 0.5) * 2]
  study_dt[n_acts < 0, n_acts := 0] # rnorm above sometimes results in negative acts. Reassign to 0.
  study_dt[, `:=`(n_acts_vi_condom    = round(n_acts * (1 - prob_ai) * prob_condom),
                  n_acts_vi_no_condom = round(n_acts * (1 - prob_ai) * (1 - prob_condom)),
                  n_acts_ai_condom    = round(n_acts * prob_ai * prob_condom),
                  n_acts_ai_no_condom = round(n_acts * prob_ai * (1 - prob_condom)))]
  
  ## Run simulation
  for(t in time_steps[2:37]) {
    
    study_dt[time == t, `:=`(age      = aging(prev_ages = study_dt[time == t - 30, age]),
                             preg     = pregnancy(dt = study_dt[time %in% (t - 270):t, .(id, time, preg)], t = t),
                             adh      = assign_adh_fup(dt = study_dt[time %in% c(t, t - 30, t - 60, t - 300), .(id, time, arm, n_acts, n_acts_vi_condom, n_acts_ai_condom, n_acts_vi_no_condom, n_acts_ai_no_condom, adh, preg)], t = t),
                             on_study = study_status(dt = study_dt[time == t - 30, .(id, censor_time, time, on_study)], t = t, parms = params),
                             vl_p1    = viral_load(vl = study_dt[time == t - 30, vl_p1]),
                             vl_p2    = viral_load(vl = study_dt[time == t - 30, vl_p2]),
                             hiv      = hiv_transmission(dt = study_dt[time == t - 30, .(id, time, age, vl_p1, vl_p2, adh, hiv, arm, n_acts_vi_condom, n_acts_vi_no_condom, n_acts_ai_condom, n_acts_ai_no_condom)], t = t, l = params$lambda))]
    
    # make run_sim a function that saves long dataset to be used in Cox proportional hazards models. Then can call run_sim within other functions that will save other output (e.g., for ABC, sum infections by arm by age)
  }
  
  ## Clean up: Remove all rows where on_study == 0
  study_dt <- study_dt[on_study == 1]
  
  ## Clean up: Among seroconverters, remove rows subsequent to visit at which HIV is first detected
  study_dt <- study_dt[, .SD[cumsum(hiv) <= 1], by = id]
  
  return(study_dt)
  
  ## TO DO: Check that male partner HIV prevalence distribution assigned in model matches observed age- and country-specific prevalence distribution
}
