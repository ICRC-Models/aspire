#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 March 2019
# calculate_acts: Function to calculate the number and type (VI, AI, adherent, partially adherent, non-adherent) of acts per woman per month.
#
# input:  Subsets of m_dt, f_dt
# output: Count of number of acts of each type per woman
# 
#######################################################################################

calculate_acts <- function(f_dt, m_dt) {
  
  n_rel_dt <- data.table(table(m_dt))
  setnames(n_rel_dt, old = c("m_dt", "N"), new = c("id", "n_part"))
  n_rel_dt[, id := as.numeric(id)]
  
  dt <- merge(x = n_rel_dt, y = data.table(id = 1:2614), by = "id", all = T)
  
  dt <- merge(x = n_rel_dt, y = f_dt, by = "id", all = T)
  dt[is.na(n_part), n_part := 0]
  
  # Women in the placebo arm have values of NA for adherence. Set these to 3 for calculation purposes.
  dt[is.na(adh), adh := 3]
  
  # Calculate quantities of interest
  dt[, `:=`(n_vi_acts = pp_vi_acts * n_part,
            n_ai_acts = pp_ai_acts * n_part,
            n_acts_full_adh    = ifelse(adh == 1, (pp_vi_acts + pp_ai_acts) * n_part, 0),
            n_acts_partial_adh = ifelse(adh == 2, (pp_vi_acts + pp_ai_acts) * n_part, 0),
            n_acts_non_adh     = ifelse(adh == 3, (pp_vi_acts + pp_ai_acts) * n_part, 0))]
  
  return(dt[, .(n_vi_acts, n_ai_acts, n_acts_full_adh, n_acts_partial_adh, n_acts_non_adh)])
}
