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
  
  # Merge simulated data with proportion adherent by site and age group
  dt <- merge(x = dt, y = prop_adh_dt[, .SD, .SDcols = c("site", "age_cat", "post_adh_intervention", prop_adh)], by = c("site", "age_cat", "post_adh_intervention"))
  setnames(dt, old = prop_adh, new = "prop_adh")
  
  setkey(dt, id, visit)
  
  # Assign adherence in the pre-adherence intervention time period with a binomial draw
  dt[post_adh_intervention == 0, adh := rbinom(n = length(unique(id)), size = 1, prob = prop_adh), by = id]
  
  # By site and age group, calculate the absolute percentage difference in adherence from pre-adherence intervention time period to post-adherence intervention time period
  adh_diff_dt <- dcast(data = prop_adh_dt[, .SD, .SDcols = c("site", "age_cat", "post_adh_intervention", prop_adh)], formula = site + age_cat ~ post_adh_intervention, value.var = prop_adh)
  adh_diff_dt[, increase_adh := `1` - `0`]
  
  adh_diff_dt <- merge(x = adh_diff_dt, y = dt[, .(n_total = length(unique(id)), n_adh = sum(visit == 0 & adh == 1)), by = c("site", "age_cat")], by = c("site", "age_cat"))
  adh_diff_dt[, n_to_sample := ifelse(round(n_total * increase_adh) + n_adh > n_total, n_total - n_adh, round(n_total * increase_adh))]
  
  ids_newly_adh <- unname(unlist(sapply(dt[, unique(site)], function(x) {
    sapply(dt[, unique(age_cat)], function(y) {
      if(adh_diff_dt[site == x & age_cat == y, n_to_sample > 0]) { sample(x = dt[site == x & age_cat == y & adh == 0, unique(id)], size = adh_diff_dt[site == x & age_cat == y, n_to_sample]) }
    })
  })))
  
  # For women identified to be newly adherent, give their post-adherence intervention time period a value of 1 for adherent
  dt[id %in% ids_newly_adh & post_adh_intervention == 1, adh := 1]
  
  # For all other women, carry forward their value assigned for the pre-adherence intervention time period
  dt[, adh := na.locf(adh), by = id]
  
  return(dt[, adh])
}
