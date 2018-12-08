# Purpose: Estimate beta parameters to be used in ABC model-fitting procedure for per-act risk of vaginal intercourse, using point estimate and confidence interval reported in Boily, et. al.'s 2009 meta-analysis.
# Author:  Kathryn Peebles (adapted from script from Anna Bershteyn)
# Date:    3 December 2018

# setwd("~/Dropbox/KenyaPrevData")
# 
# library(data.table)

lambda_mean <- 0.0030
lambda_lb   <- 0.0014
lambda_ub   <- 0.0063

alpha_initial <- 1000

to_minimize <- function(alpha) {
  alpha <- abs(alpha)
  beta_computed <- alpha * ((1 - lambda_mean)/lambda_mean)
  penalty <- (lambda_lb - qbeta(0.025, alpha, beta_computed, ncp = 0, lower.tail = T, log.p = F)) ^ 2 + (lambda_ub - qbeta(0.975, alpha, beta_computed, ncp = 0, lower.tail = T, log.p = F)) ^ 2
  return(penalty)
}

to_minimize(alpha_initial)

alpha_optimized <- optim(alpha_initial, to_minimize, control = list(maxiter = 10000, abstol = 0.00000001))$par

beta_optimized <- alpha_optimized * ((1 - lambda_mean)/lambda_mean)

plot(dbeta(x = seq(0, 0.01, 0.0001), shape1 = alpha_optimized, shape2 = beta_optimized), type = "l", xaxt = "n", xlab = "Per-vaginal-act probability of HIV transmission", ylab = "Density")
axis(side = 1, at = seq(0, length(seq(0, 0.01, 0.0001)), length.out = length(seq(0, 0.01, 0.002))), labels = seq(0, 0.01, 0.002))
