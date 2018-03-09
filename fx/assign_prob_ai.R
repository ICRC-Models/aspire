#######################################################################################
# 
# Kathryn Peebles
# 2017-08-15
# assign_prob_ai: function to assign individual probability of engaging in an act of AI at each sex act. Called in creation of data table to store individual participant characteristics.
#
#######################################################################################

assign_prob_ai <- function(n, prop_ai, alpha, beta) {
  ai <- rbinom(n = n, size = 1, prob = prop_ai)
  ai[ai == 1] <- rbeta(sum(ai), alpha, beta)
  return(ai)
}