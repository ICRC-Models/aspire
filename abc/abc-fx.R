#######################################################################################
# 
# Author:  Kathryn Peebles
# Date:    22 March 2018
# Purpose: Functions to be used in ABC model-fitting procedure.
#
#######################################################################################

#######################################################################################
# 
# select_particles: Function to select particles for iteration t. Particles are sampled from iteration t-1 with probability equal to their normalized weight.
#
# input:  data.table of N-alpha particles from iteration t-1.
# output: data.table of particles to be simulated in iteration t
#
#######################################################################################

select_particles <- function(dt, n_new_particles) {
  particles <- dt[, sample(x = .I, size = n_new_particles, replace = T, prob = w)]
  dt <- dt[particles, ]
  return(dt)
}

#######################################################################################
# 
# perturb_particles: Function to perturb each of the particles selected from iteration t-1 using a multivariate Gaussian perturbation particle with covariance equal to two times the weighted empirical covariance of included particles in t-1. Ensure that all perturbations are within the specified priors to ensure weights of subsequent particles > 0.
#
# input:  data.table of particles selected in call to select_particles and covariance matrix estimated from particles in iteration t-1.
# output: data.table of particles to be simulated in iteration t.
#
#######################################################################################

perturb_particles <- function(dt, cov_mat) {
  dt[, in_prior := F]
  
  while(dt[, any(!in_prior)]) {
    not_in_prior <- dt[, .I[!(in_prior)]]
    
    for(i in not_in_prior) {
      dt[i, c("lambda_new", "cond_rr_new", "c_new", "s_new", "rr_ai_new") := as.list(rmnorm(n = 1, mean = as.vector(c(lambda, cond_rr, c, s, rr_ai)), varcov = cov_mat))]
      dt[i, in_prior := is_in_prior(lambda = lambda_new, cond_rr = cond_rr_new, c = c_new, s = s_new, rr_ai = rr_ai_new)]
    }
  }
  
  return(dt)
}

#######################################################################################
# 
# is_in_prior: Function to test is perturbed particle is within the bounds of the prior distributions. 
#
# input:  Values of perturbed values of prop_ai, a, and b parameters
# output: Logical value (T/F) indicating whether particle is in bounds of priors
#
#######################################################################################

is_in_prior <- function(lambda, cond_rr, c, s, rr_ai) {
  d <- dbeta(x = lambda, shape1 = priors[param == "lambda", alpha], shape2 = priors[param == "lambda", beta]) *
       dunif(x = cond_rr, min = priors[param == "cond_rr", min], max = priors[param == "cond_rr", max]) * 
       dunif(x = c,       min = priors[param == "c",       min], max = priors[param == "c",       max]) *
       dunif(x = s,       min = priors[param == "s",       min], max = priors[param == "s",       max]) *
       dunif(x = rr_ai,   min = priors[param == "rr_ai",   min], max = priors[param == "rr_ai", max])
  
  return(as.logical(d))
}

#######################################################################################
# 
# calc_weights: Function to calculate the weights of all particles in iteration t
#
# input:  data.table of particles from iteration t-1, vector of weights for particles in iteration t-1, and data.table of particles from iteration t
# output: vector of weights of length equal to number of particles in iteration t
#
#######################################################################################

calc_weights <- function(prev_iter_particles, weights, curr_iter_particles) {
  numerator   <- calc_numerator_weights(curr_iter_particles = curr_iter_particles)
  denominator <- calc_denominator_weights(prev_iter_particles = as.matrix(prev_iter_particles), 
                                          weights = weights, 
                                          curr_iter_particles = as.matrix(curr_iter_particles[, .(lambda, cond_rr, c, s, rr_ai)]))
  return(numerator/denominator)
}


#######################################################################################
# 
# calc_numerator_weights: Function to calculate the denominator of weight calculation for all particles in iteration t
#
# input:  data.table of particles from iteration t
# output: vector of numerator values of length equal to number of particles in iteration t
#
#######################################################################################

calc_numerator_weights <- function(curr_iter_particles) {
  curr_iter_particles[, d_prior := dbeta(x = lambda, shape1 = priors[param == "lambda", alpha], shape2 = priors[param == "lambda", beta]) *
                        dunif(x = cond_rr, min = priors[param == "cond_rr", min], max = priors[param == "cond_rr", max]) * 
                        dunif(x = c, min = priors[param == "c", min], max = priors[param == "c", max]) *
                        dunif(x = s, min = priors[param == "s", min], max = priors[param == "s", max]) *
                        dunif(x = rr_ai, min = priors[param == "rr_ai", min], max = priors[param == "rr_ai", max])]
  return(curr_iter_particles[, d_prior])
}

#######################################################################################
# 
# calc_denominator_weights: Function to calculate the denominator of weight calculation for all particles in iteration t
#
# input:  data.table of particles from iteration t-1, vector of weights for particles in iteration t-1, and data.table of particles from iteration t
# output: vector of denominator values of length equal to number of particles in iteration t
#
#######################################################################################

calc_denominator_weights <- function(prev_iter_particles, weights, curr_iter_particles) {
  # Specify dimensions (number of parameters) of prev_iter_particles
  d <- ncol(prev_iter_particles)
  
  # Specify covariance matrix and inverse covariance matrix of particles in iteration t-1
  cov_mat <- as.matrix(2 * cov.wt(x = as.matrix(prev_iter_particles), wt = weights)$cov)
  inv_cov_mat <- solve(cov_mat)
  
  # Specify constant multiplier in multivariate normal distribution
  const <- 1/sqrt((2 * pi) ^ d * det(cov_mat))
  
  # Create empty vector in which to store denominator values for particles in iteration t
  denom <- numeric(nrow(curr_iter_particles))
  
  # For each particle simulated in iteration t, calculate the probability of moving to that particle from each accepted particle in iteration t-1
  for(i in 1:nrow(curr_iter_particles)) {
    for(j in 1:nrow(prev_iter_particles)) {
      denom[i] <- denom[i] + weights[j] * const * exp(-0.5 * (curr_iter_particles[i,] - prev_iter_particles[j,]) %*% inv_cov_mat %*% (curr_iter_particles[i,] - prev_iter_particles[j,]))
    }
  }
  
  return(denom)
}

