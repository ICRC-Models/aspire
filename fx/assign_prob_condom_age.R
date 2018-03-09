#######################################################################################
# 
# Kathryn Peebles
# 2017-11-19
# assign_prob_condom_age: function to assign individual probability of condom use. Called in creation of data table to store individual participant characteristics.
#
# input:  Subset of f_dt containing only "age_cat" variable
# output: Vector of values bounded in [0, 1] indicating individual-specific per-act probability of condom use
#
#######################################################################################

assign_prob_condom_age <- function(age_cats,
                                   a_pc_18_21, b_pc_18_21,
                                   a_pc_22_26, b_pc_22_26,
                                   a_pc_27_45, b_pc_27_45) {
  
  prob_condom <- vector(mode = "numeric", length = length(age_cats))
  prob_condom[age_cats == "18-21"] <- rbeta(n      = sum(age_cats == "18-21"),
                                            shape1 = a_pc_18_21,
                                            shape2 = b_pc_18_21)
  prob_condom[age_cats == "22-26"] <- rbeta(n      = sum(age_cats == "22-26"),
                                            shape1 = a_pc_22_26,
                                            shape2 = b_pc_22_26)
  prob_condom[age_cats == "27-45"] <- rbeta(n      = sum(age_cats == "27-45"),
                                            shape1 = a_pc_27_45,
                                            shape2 = b_pc_27_45)
  
  return(prob_condom)
}