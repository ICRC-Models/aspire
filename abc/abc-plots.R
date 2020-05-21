## Purpose: Create plot that shows distribution of simulated particles
## Author:  Kathryn Peebles
## Date:    7 January 2019

# Attach packages
library(ggplot2)
library(ggtern)
library(data.table)
library(ggpubr)
library(scales)
library(viridis)

setwd("~/Documents/code/aspire/abc")

# Plot line showing density of particles for each iteration
load(paste0(getwd(), "/priors.RDATA"))

n_distributions <- 14 # number of distributions + 1 for iteration 0

for(i in 0:(n_distributions - 1)) { # n_distributions - 1 to start with 0
  load(paste0(getwd(), "/hyak/t", i, "/particles_acc_", i, ".RDATA"))
}

params_vec <- c("lambda", "cond_rr", "c", "s", "rr_ai", "p_rate_rr", "base_male_hiv_incidence")
x_lab_vec <- c("Baseline probability of HIV transmission", 
               "Relative risk of HIV-positive partner\namong women reporting condom use at last sex", 
               "Carrying capacity term of inverse logistic\ndistribution of age-dependent relative risk\nof HIV-positive partner", 
               "Scale term of inverse logistic distribution\nof age-dependent relative risk\nof HIV-positive partner", 
               "Relative risk of anal intercourse",
               "Relative increase in partner change rate",
               "Reference male HIV-1 incidence")

#### CALCULATE P_ACC FOR ALL SAMPLES ####
for(dist_n in 0:(n_distributions - 1)) {
  print(paste0("Distribution: ", dist_n))
  print(get(x = paste0("particles_acc_", dist_n))[, sum(t == dist_n)/.N])
}

#### DENSITY PLOTS OF ALL PARAMETERS ####
cols <- viridis(25)
cols <- rev(cols[6:20])

for(i in 1:length(params_vec)) {
  # Create data table with density of prior distribution
  if(params_vec[i] == "lambda") {
    prior_d_dt <- data.table(x = seq(particles_acc_0[, min(.SD), .SDcols = params_vec[i]], particles_acc_0[, max(.SD), .SDcols = params_vec[i]], length.out = 60000),
                             y = dbeta(x = seq(particles_acc_0[, min(.SD), .SDcols = params_vec[i]], particles_acc_0[, max(.SD), .SDcols = params_vec[i]], length.out = 60000), shape1 = priors[param == params_vec[i], alpha], shape2 = priors[param == params_vec[i], beta]),
                             dist = "Prior")
  } else {
    prior_d_dt <- data.table(x = seq(particles_acc_0[, min(.SD), .SDcols = params_vec[i]], particles_acc_0[, max(.SD), .SDcols = params_vec[i]], length.out = 60000),
                             y = dunif(x = seq(particles_acc_0[, min(.SD), .SDcols = params_vec[i]], particles_acc_0[, max(.SD), .SDcols = params_vec[i]], length.out = 60000), min = priors[param == params_vec[i], min], max = priors[param == params_vec[i], max]),
                             dist = "Prior")
  }
  
  # Create data table with weighted density of each intermediate distribution
  densities_dt_list <- lapply(0:(n_distributions - 1), function(iter) 
    data.table(x    = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w, adjust = 2), .SDcols = c(params_vec[i], "w")]$x,
               y    = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w, adjust = 2), .SDcols = c(params_vec[i], "w")]$y,
               dist = ifelse(iter == (n_distributions - 1), "Posterior distribution", paste0("Intermediate distribution ", iter))))
  
  # Combine above two data tables
  densities_dt <- rbindlist(l = densities_dt_list)
  densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))
  
  # Factor variable dist so that colors and legend values are plotted in desired order
  densities_dt[, dist := factor(x = dist, levels = c("Prior", sapply(0:(n_distributions - 2), function(x) paste0("Intermediate distribution ", x)), "Posterior distribution"))]
  
  densities_dt[dist == "Posterior distribution", prob_y := y/(sum(y))]
  
  cri <- round(quantile(x = sample(x = densities_dt[!is.na(prob_y), x], size = 1000, replace = T, prob = densities_dt[!is.na(prob_y), prob_y]), probs = c(0.05, 0.95)), 3)
  
  # Create ggplot
  density_plot <- ggplot() +
    geom_line(data = densities_dt, aes(x = x, y = y, group = dist, col = dist), size = 1) +
    geom_vline(xintercept = cri, linetype = "dashed") +
    scale_color_manual(values = cols) +
    labs(x = x_lab_vec[i], y = "Density", col = "") +
    theme(axis.text = element_text(size = 11, family = "Times"),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 11, family = "Times"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size = 11, family = "Times"))
  
  assign(x = paste0(params_vec[i], "_plot"), value = density_plot)
  
  rm(list = c("prior_d_dt", "densities_dt_list", "densities_dt", "density_plot"))
}

ggarrange(lambda_plot, rr_ai_plot, c_plot, s_plot, cond_rr_plot, p_rate_rr_plot, base_male_hiv_incidence_plot, ncol = 2, nrow = 4, common.legend = T, legend = "right")

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Manuscript supplement/Figures/ABC_Distributions_20200214.pdf", width = 8.5, height = 7)

#### IDENTIFY BEST-FIT PARAMETERS ####
# Load particles for iteration with p_acc = 0.100
load(paste0(getwd(), "/hyak/t13/particles_acc_13.RDATA"))

# Identify best-fit parameter from univariate weighted densities for each parameter
for(i in 1:length(params_vec)) {
  d <- particles_acc_13[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]
  
  val <- d$x[which.max(d$y)]
  
  print(paste0("Best-fitting value for ", x_lab_vec[i], " is: ", val))
}

#### CALCULATE 90% CREDIBLE INTERVAL ####
iter <- n_distributions - 1

for(i in 1:length(params_vec)) {
  densities_dt <- data.table(x = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]$x,
                             y = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]$y)
  
  densities_dt[, prob_y := y/(sum(y))]
  
  cri <- round(quantile(x = sample(x = densities_dt[, x], size = 1000, replace = T, prob = densities_dt[, prob_y]), probs = c(0.05, 0.95)), 3)
  
  print(paste0("90% CrI for ", x_lab_vec[i], ": ", cri[1], ", ", cri[2]))
}
