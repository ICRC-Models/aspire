## Purpose: Function to re-randomize arm stratified by behavioral and demographic risk characteristics
## Author:  Kathryn Peebles
## Date:    1 October 2019

re_randomize_arm <- function(dt) {
  # Initialize active arm IDs with an empty vector
  dpv_ids <- numeric(0)
  
  # For each combination of site, BV status, and median acts indicator, randomize half of women to receive the dapivirine ring.
  for(i in dt[, unique(site)]) {
    for(j in 0:1) {
      for(k in 0:1) {
        n_to_randomize <- nrow(dt[site == i & b_bv == j & gt_median_acts == k])
        
        # If there are an odd number of women to randomize, randomly assign whether rounding will use the floor or ceiling function.
        if(n_to_randomize %% 2 != 0) {
          rnd_fx <- rbinom(n = 1, size = 1, prob = 0.5)
          dpv_n <- ifelse(rnd_fx == 0, floor(n_to_randomize/2), ceiling(n_to_randomize/2))
        } else {
          dpv_n <- n_to_randomize/2
        }
        
        dpv_ids <- c(dpv_ids, dt[site == i & b_bv == j & gt_median_acts == k, sample(x = id, size = dpv_n)])
      }
    }
  }
  
  return(dpv_ids)
}
