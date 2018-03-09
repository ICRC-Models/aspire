#######################################################################################
# 
# Kathryn Peebles
# 2017-11-19
# aging: function to increment age of study participants by 30 days at each time step.
#
# input:  Vector of female ages at time t - 1.
# output: Vector of female ages at time t.
#
#######################################################################################

aging <- function(prev_ages) {
  ages <- prev_ages + 30/365
  return(ages)
}