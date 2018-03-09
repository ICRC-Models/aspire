#######################################################################################
# 
# Kathryn Peebles
# 2017-11-01
# study_status: Assign study status, determined by study status at previous time point, pre-determined censoring time, and homogeneous risk of loss-to-follow-up.
#
# input:  subset of study_dt that includes time step t - 30 and variables id, censor_time, time, on_study.
# output: vector of 0/1 values indicating study status.
#
#######################################################################################

study_status <- function(dt, t, parms) {
  dt$on_study_t <- ifelse(dt$on_study == 0, 0, NA)
  dt$on_study_t <- ifelse(dt$censor_time < t, 0, dt$on_study_t)
  dt$on_study_t <- ifelse(is.na(dt$on_study_t), rbinom(n = sum(which(is.na(dt$on_study_t))), size = 1, prob = (1 - parms$loss_to_fup_prob)), dt$on_study_t)
  
  return(as.integer(dt$on_study_t))
}