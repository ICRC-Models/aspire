#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   10 November 2017
# viral_load: function to assign viral load to HIV-positive male partners. This function should eventually account for viral load progression and acute HIV infections. However, as a first pass, viral load is assumed fixed, and the baseline viral load is carried forward to all subsequent time points.
#
# input:  viral load values at time t - 30
# output: viral load values at time t
#
#######################################################################################

## TO DO: Function currently carries forward viral load value from previous time step. Could incorporate viral load progression in future version of model.

viral_load <- function(vl) {
  return(vl)
}