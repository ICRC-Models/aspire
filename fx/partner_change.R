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

partner_change <- function(m_dt, f_dt, cond_rr) {
  n_part_dt <- as.data.table(m_dt[active == 1, table(id)])
  n_part_dt[, id := as.integer(id)]
  
  # If any study participants have no active partnerships, they are not included in n_part_dt above. Add id back to table so that they are included in partnership formation below.
  if(length(setdiff(f_dt$id, n_part_dt$id)) != 0) {
    n_part_dt <- rbind(n_part_dt, data.table(id = setdiff(f_dt$id, n_part_dt$id), N = as.integer(0)))
  }
  
  f_dt <- merge(x = n_part_dt, y = f_dt, by = "id")
  
  # Partnership formation. If number of partners is less than maximum number of partners, form a new partnership with female age-specific probability
  f_dt <- merge(x = f_dt, y = p_rates_dt, by = "f_age")
  f_dt[N < max_part, new_part := rbinom(n    = nrow(f_dt[N < max_part]),
                                        size = 1, 
                                        prob = prob)]
  
  # Partnership dissolution. For each active relationship, end relationship with female age-specific probability
  m_dt <- merge(x = m_dt, y = p_rates_dt, by = "f_age")
  m_dt[active == 1 & max_part != 1, 
       active := rbinom(n    = nrow(m_dt[active == 1 & max_part != 1]),
                        size = 1,
                        prob = 1 - prob)]
  
  # Add new partners to m_dt and assign country, age, HIV status, and viral load. New partners are added to m_dt after determination of relationship formation and dissolution so that newly-formed relationships are not dissolved in the same time step (i.e., without contributing any risk).
  m_dt <- rbind(m_dt[, .SD, .SDcols = names(m_dt)[which(names(m_dt) != "prob")]], f_dt[new_part == 1, .(id, f_age, country, max_part, b_condom_lweek, m_hiv_rr)], fill = T)
  
  m_dt[is.na(active), active := as.integer(1)]
  
  m_dt[is.na(m_age), m_age := assign_male_age(m_ids = m_dt[is.na(m_age), id], f_age = f_dt$f_age)]
  m_dt[is.na(hiv),   hiv   := assign_male_hiv_status(dt = m_dt[is.na(hiv), .(country, m_age, f_age, b_condom_lweek, m_hiv_rr)], cond_rr = cond_rr)]
  m_dt[is.na(vl),    vl    := assign_male_vl_trt_unknown(dt = m_dt[is.na(vl), .(m_age, hiv)])]
  
  return(m_dt)
}