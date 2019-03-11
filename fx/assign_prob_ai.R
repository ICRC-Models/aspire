#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 March 2019
# assign_prob_ai: function to assign individual probability of engaging in an act of AI at each sex act given a population prevalence of AI and mean proportion of acts that are AI among women engaging in AI.
#
#######################################################################################

assign_prob_ai <- function(ids, prev_ai, prop_ai) {
  ai <- rbinom(n = length(ids), size = 1, prob = prev_ai)
  
  ab_params <- est_beta_params(mu = prop_ai, var = 0.005)
  
  ai[ai == 1] <- rbeta(sum(ai), ab_params$alpha, ab_params$beta)
  
  ai_dt <- data.table(id = ids, prob_ai = ai)
  
  return(ai_dt)
}

#######################################################################################
#
# est_beta_params: function to estimate the parameters of a beta distribution given the mean and variance of the distribution. From: https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
#
#######################################################################################

est_beta_params <- function(mu, var) {
  alpha <- ((1 - mu)/var - 1/mu) * mu ^ 2
  beta  <- alpha * (1/mu - 1)
  
  return(params = list(alpha = alpha, beta = beta))
}