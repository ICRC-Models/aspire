#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   15 March 2018
# partner_change: Function to form and terminate relationships. Formation probabilities obtained from McGrath, 2017 are applied as dissolution probabilities so that number of partners per participant is about the same throughout the simulation.
#
# input:  m_dt and subset of f_dt including id, country, age, and maximum number of partners.
# output: Modified m_dt with updated active status and appended rows for each newly formed relationship.
# 
#######################################################################################

partner_change <- function(m_dt, f_dt, cond_rr, debug_sim_hiv = F, debug_sim_hiv_acts = F) {
  # browser()
  n_part_dt <- as.data.table(m_dt[active == 1, table(id)])
  n_part_dt[, id := as.integer(id)]
  
  # If any study participants have no active partnerships, they are not included in n_part_dt above. Add id back to table so that they are included in partnership formation below.
  if(length(setdiff(f_dt$id, n_part_dt$id)) != 0) {
    n_part_dt <- rbind(n_part_dt, data.table(id = setdiff(f_dt$id, n_part_dt$id), N = as.integer(0)))
  }
  
  f_dt <- merge(x = n_part_dt, y = f_dt, by = "id")
  #print(paste0("Number unique ids: ", f_dt[, length(unique(id))]))
  # Partnership formation. If number of partners is less than maximum number of partners, form a new partnership with female age-specific probability
  f_dt <- merge(x = f_dt, y = p_rates_dt, by = "f_age")
  #print(paste0("Number possible rels to gain: ", nrow(f_dt[N < max_part])))
  f_dt[N < max_part, new_part := rbinom(n    = nrow(f_dt[N < max_part]),
                                        size = 1, 
                                        prob = prob)]
  #print(paste0("Number rels gained: ", f_dt[, sum(new_part, na.rm = T)]))
  
  # Partnership dissolution. For each woman in an active relationship, end a single relationship with female age-specific probability.
  m_dt[active == 1, rel_n := 1:.N, by = id] # Count all active relationships
  m_dt <- merge(x = m_dt, y = p_rates_dt, by = "f_age")
  #print(paste0("Number unique ids in m_dt: ", m_dt[, length(unique(id))]))
  #print(paste0("Number possible rels to lose: ", nrow(m_dt[active == 1 & max_part != 1 & rel_n == 1])))
  #n_active_rels <- m_dt[, sum(active)]
  m_dt[active == 1 & max_part != 1 & rel_n == 1, 
       active := rbinom(n    = nrow(m_dt[active == 1 & max_part != 1 & rel_n == 1]),
                        size = 1,
                        prob = 1 - prob)]
  #print(paste0("Number rels lost: ", m_dt[, sum(active)] - n_active_rels))
  
  # Add new partners to m_dt and assign country, age, HIV status, and viral load. New partners are added to m_dt after determination of relationship formation and dissolution so that newly-formed relationships are not dissolved in the same time step (i.e., without contributing any risk).
  m_dt <- rbind(m_dt[, .SD, .SDcols = names(m_dt)[!(names(m_dt) %in% c("prob", "rel_n"))]], f_dt[new_part == 1, .(id, arm, f_age, country, max_part, b_condom_lweek, m_hiv_rr)], fill = T)
  
  m_dt[is.na(active), active := as.integer(1)]
  #print(paste0("Number active rels at end of timestep: ", m_dt[, sum(active)]))
  
  m_dt[is.na(m_age), m_age := assign_male_age(m_ids = m_dt[is.na(m_age), id], f_age = f_dt$f_age)]
  
  if(debug_sim_hiv | debug_sim_hiv_acts) {
    m_dt[is.na(hiv), hiv := as.integer(1)]
    m_dt[is.na(vl),  vl  := 3.0]
  } else {
    m_dt[is.na(hiv),   hiv   := assign_male_hiv_status(dt = m_dt[is.na(hiv), .(country, m_age, f_age, b_condom_lweek, m_hiv_rr)], cond_rr = cond_rr)]
    m_dt[is.na(vl),    vl    := assign_male_vl_trt_unknown(dt = m_dt[is.na(vl), .(m_age, hiv)])]
  }
  
  return(m_dt)
}