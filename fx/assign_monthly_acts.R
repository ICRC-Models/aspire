#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   11 November 2018
#
# assign_mvi_acts: Function to assign a monthly number of vaginal sex acts as a Poisson random draw with mean equal to average monthly vaginal acts reported at baseline.
#
# input:  Subset of f_dt.
# output: Vector of number of vaginal sex acts per woman per partner per month to be used in simulation.
#
# assign_mai_acts: Function to assign a monthly number of anal sex acts as a Poisson random draw with mean equal to reported average monthly acts.
#
# input:  Subset of f_dt.
# output: Vector of number of anal sex acts per woman per partner per month to be used in simulation.
# 
#######################################################################################

assign_mvi_acts <- function(dt) {
  dt[, pp_mvi_acts := rpois(n = nrow(dt), lambda = pp_vi_acts_avg)]
  return(dt[, pp_mvi_acts])
}


assign_mai_acts <- function(dt) {
  dt[pp_ai_acts_avg == 0, pp_mai_acts := 0]
  dt[pp_ai_acts_avg != 0, pp_mai_acts := as.double(rpois(n = nrow(dt[pp_ai_acts_avg != 0]), lambda = pp_ai_acts_avg))]
  return(dt[, pp_mai_acts])
}