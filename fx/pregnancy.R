## TO DO: We estimate HIV as per usual, but upon re-joining study, censor if HIV = 1. Talk to EB about this function- did HIV testing continue at monthly intervals for pregnant women?
## Were pregnancy rates different by age? If more common among younger women, could contribute to lower efficacy in that group.

#######################################################################################
# 
# Kathryn Peebles
# 2017-11-10
# pregnancy: Pregnancy occurred at rate of 3.9 per 100 person-years in the dapivirine group and 4.0 per 100 person-years in the placebo group. As this difference is statistically not significant, we model a pregnancy rate of 3.95 per 100 person-years. 30-day risk of pregnancy is 1 - exp(-1 * (3.95/100) * (30/365)).
#
# input:  Subset of study_dt that includes current time step and previous nine time steps.
# output: Vector of 0/1 values indicating pregnancy status.
#
#######################################################################################

pregnancy <- function(dt, t) {
  
  # If not pregnant at previous time step, assign pregnancy with probability equal to 1 - exp(-1 * (3.95/100) * (30/365))
  dt[time == t & id %in% dt[time == t - 30 & preg == 0, id], 
     preg := rbinom(n = length(id %in% dt[time == t - 30 & preg == 0, id]), size = 1, prob = (1 - exp(-1 * (3.95/100) * (30/365))))]
  
  # If pregnant at previous time step, assign pregnancy on basis of number of months pregnant. If less than 9, pregnancy status is 1. If equal to 9, pregnancy status is 0. This assumes that all women with incident pregnancies deliver term infants.
  dt[time == t & id %in% dt[time == t - 30 & preg == 1, id],
     preg := sapply(dt[time == t - 30 & preg == 1, id], function(x) {
       ifelse(test = dt[time %in% (t - 270):(t - 30) & id == x, sum(preg) == 9],
              yes  = as.integer(0),
              no   = as.integer(1))
     })]
  
  return(dt[time == t, preg])
}
