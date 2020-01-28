#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 March 2019
#
# assign_adh: Function to assign categorical adherence value at baseline and follow-up according to systematically varied proportion of women with consistent and no adherence.
#
# input:  Subset of data in f_dt
# output: Vector 0/1 values indicating adherence at time == t
#
#######################################################################################

assign_adh <- function(dt, target_adh) {
  
  # If target adherence coverage is 100%, assign all visits as adherent. Otherwise, use predictive model to assign adherence.
  if(target_adh == 1) {
    dt[, adh := 1]
  } else {
    
    # Calculate adherence for visit 0 for all women
    dt[visit == 0, adh := rbinom(n = nrow(dt[visit == 0]), size = 1, prob = prob_adh)]
    
    # For women with assigned adherence of 1 at visit 0, carry forward value of 1 for remaining visits.
    ids_adh <- dt[adh == 1, unique(id)]
    dt[id %in% ids_adh, adh := na.locf(adh, na.rm = F), by = id]
    
    # For women with assigned adherence of 0 at visit 0, there are two opportunities to become adherent: following implementation of adherence interventions and after 3 months of study participation. Carry forward value of 0 for duration of period 1, then conduct new random binomial draw for period 2.
    ids_non_adh <- dt[adh == 0, unique(id)]
    dt[id %in% ids_non_adh & adh_period == 1, adh := na.locf(adh, na.rm = F), by = id]
    
    dt[id %in% ids_non_adh & adh_period == 2, temp_n := 1:.N, by = id]
    dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1, adh := rbinom(n = nrow(dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1]), size = 1, prob = incr_prob_adh)]
    
    # Carry forward adherence values of 1 and 0 newly assigned to period 2
    ids_adh <- dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1 & adh == 1, unique(id)]
    dt[id %in% ids_adh, adh := na.locf(adh, na.rm = F), by = id]
    
    ids_non_adh <- dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1 & adh == 0, unique(id)]
    dt[id %in% ids_non_adh & adh_period == 2, adh := na.locf(adh, na.rm = F), by = id]
    
    # Assign adherence values for period 3 among those who were still non-adherent in period 2
    dt[, temp_n := NULL]
    dt[id %in% ids_non_adh & adh_period == 3, temp_n := 1:.N, by = id]
    dt[id %in% ids_non_adh & adh_period == 3 & temp_n == 1, adh := rbinom(n = nrow(dt[id %in% ids_non_adh & adh_period == 3 & temp_n == 1]), size = 1, prob = incr_prob_adh)]
    
    # Carry forward adherence values of 1 and 0 newly assigned to period 3
    dt[, adh := na.locf(adh, na.rm = F), by = id]
    
    # Calculate proportion of visits at which women are adherent (approximate, since observed time on study per obs_on_study will vary from simulation to simulation)
    prop_adh_current <- dt[obs_on_study == 1, mean(adh)]
    
    if(target_adh > prop_adh_current) {
      
      # Assign the new probability of adherence as the difference between the current proportion adherent and the increase in probablity to meet target adherence coverage. 
      dt[, new_prob_adh := (prob_adh * (target_adh/prop_adh_current)) - prob_adh]
      
      while(prop_adh_current < (target_adh - 0.015)) {
        # Reassign adherence among those initially assigned non-adherent
        dt[adh == 0, adh := NA]
        
        # Assign adherence at visit 0
        dt[visit == 0 & is.na(adh), adh := rbinom(n = nrow(dt[visit == 0 & is.na(adh)]), size = 1, prob = new_prob_adh)]
        
        # For women with assigned adherence of 1 at visit 0, carry forward value of 1 for remaining visits.
        ids_adh <- dt[adh == 1, unique(id)]
        dt[id %in% ids_adh, adh := na.locf(adh, na.rm = F), by = id]
        
        # For women with assigned adherence of 0 at visit 0, there are two opportunities to become adherent: following implementation of adherence interventions and after 3 months of study participation. Carry forward value of 0 for duration of period 1, then conduct new random binomial draw for period 2.
        ids_non_adh <- dt[adh == 0, unique(id)]
        dt[id %in% ids_non_adh & adh_period == 1, adh := na.locf(adh, na.rm = F), by = id]
        
        dt[, temp_n := NULL]
        dt[id %in% ids_non_adh & adh_period == 2, temp_n := 1:.N, by = id]
        dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1, adh := rbinom(n = nrow(dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1]), size = 1, prob = new_prob_adh)]
        
        # Carry forward adherence values of 1 and 0 newly assigned to period 2
        ids_adh <- dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1 & adh == 1, unique(id)]
        dt[id %in% ids_adh, adh := na.locf(adh, na.rm = F), by = id]
        
        ids_non_adh <- dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1 & adh == 0, unique(id)]
        dt[id %in% ids_non_adh & adh_period == 2, adh := na.locf(adh, na.rm = F), by = id]
        
        # Assign adherence values for period 3 among those who were still non-adherent in period 2
        dt[, temp_n := NULL]
        dt[id %in% ids_non_adh & adh_period == 3, temp_n := 1:.N, by = id]
        dt[id %in% ids_non_adh & adh_period == 3 & temp_n == 1, adh := rbinom(n = nrow(dt[id %in% ids_non_adh & adh_period == 3 & temp_n == 1]), size = 1, prob = new_prob_adh)]
        
        # Carry forward adherence values of 1 and 0 newly assigned to period 3
        dt[, adh := na.locf(adh, na.rm = F), by = id]
        
        prop_adh_current <- dt[obs_on_study == 1, mean(adh)]
      }
    } else {
      
      # Assign the new probability of adherence as the difference between the current proportion adherent and the increase in probablity to meet target adherence coverage. 
      dt[, new_prob_adh := prob_adh * (target_adh/prop_adh_current)]
      
      # Re-assign value of increased probability of adherence at periods 2 and 3
      dt[, incr_prob_adh := c(0, diff(new_prob_adh)), by = id]
      
      # Reset prop_adh_current to 0  for the while loop below to function correctly
      prop_adh_current <- 0
      
      while(prop_adh_current < (target_adh - 0.015)) {
        # Reassign adherence for all participants
        dt[, adh := NULL]
        
        # Calculate adherence for visit 0 for all women
        dt[visit == 0, adh := rbinom(n = nrow(dt[visit == 0]), size = 1, prob = new_prob_adh)]
        
        # For women with assigned adherence of 1 at visit 0, carry forward value of 1 for remaining visits.
        ids_adh <- dt[adh == 1, unique(id)]
        dt[id %in% ids_adh, adh := na.locf(adh, na.rm = F), by = id]
        
        # For women with assigned adherence of 0 at visit 0, there are two opportunities to become adherent: following implementation of adherence interventions and after 3 months of study participation. Carry forward value of 0 for duration of period 1, then conduct new random binomial draw for period 2.
        ids_non_adh <- dt[adh == 0, unique(id)]
        dt[id %in% ids_non_adh & adh_period == 1, adh := na.locf(adh, na.rm = F), by = id]
        
        dt[id %in% ids_non_adh & adh_period == 2, temp_n := 1:.N, by = id]
        dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1, adh := rbinom(n = nrow(dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1]), size = 1, prob = incr_prob_adh)]
        
        # Carry forward adherence values of 1 and 0 newly assigned to period 2
        ids_adh <- dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1 & adh == 1, unique(id)]
        dt[id %in% ids_adh, adh := na.locf(adh, na.rm = F), by = id]
        
        ids_non_adh <- dt[id %in% ids_non_adh & adh_period == 2 & temp_n == 1 & adh == 0, unique(id)]
        dt[id %in% ids_non_adh & adh_period == 2, adh := na.locf(adh, na.rm = F), by = id]
        
        # Assign adherence values for period 3 among those who were still non-adherent in period 2
        dt[, temp_n := NULL]
        dt[id %in% ids_non_adh & adh_period == 3, temp_n := 1:.N, by = id]
        dt[id %in% ids_non_adh & adh_period == 3 & temp_n == 1, adh := rbinom(n = nrow(dt[id %in% ids_non_adh & adh_period == 3 & temp_n == 1]), size = 1, prob = incr_prob_adh)]
        
        # Carry forward adherence values of 1 and 0 newly assigned to period 3
        dt[, adh := na.locf(adh, na.rm = F), by = id]
        
        prop_adh_current <- dt[obs_on_study == 1, mean(adh)]
      }
    }
  }
  
  return(dt[, adh])
}
