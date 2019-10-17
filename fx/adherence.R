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

assign_adh <- function(dt, prop_adh) {
  
  # Create age category to merge on
  dt[, age_cat := ifelse(f_age <= 26, "18-26", "27-45")]
  
  # Merge simulated data were proportion adherent by site and age group
  dt <- merge(x = dt, y = prop_adh_dt[, .SD, .SDcols = c("site", "age_cat", prop_adh)], by = c("site", "age_cat"))
  setnames(dt, old = prop_adh, new = "prop_adh")
  
  setkey(dt, id, visit)
  
  dt[, adh := rbinom(n = length(unique(id)), size = 1, prob = prop_adh), by = id]
  
  return(dt[, adh])
}
