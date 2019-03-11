## Purpose: Create plot that shows distribution of simulated particles
## Author:  Kathryn Peebles
## Date:    7 January 2019

library(ggplot2)
library(ggtern)
library(data.table)
library(ggpubr)
library(scales)

setwd("~/Documents/code/aspire/abc")

# Plot line showing density of particles for each iteration
load(paste0(getwd(), "/priors.RDATA"))

n_distributions <- 24

for(i in 0:(n_distributions - 1)) { # n_distributions - 1 to start with 0
  load(paste0(getwd(), "/hyak/t", i, "/particles_acc_", i, ".RDATA"))
}

params_vec <- c("lambda", "cond_rr", "c", "s", "rr_ai")
x_lab_vec <- c("Baseline probability of HIV transmission", 
               "Relative risk of HIV-positive partner\namong women reporting condom use at last sex", 
               "Carrying capacity term of inverse logistic\ndistribution of age-dependent relative risk\nof HIV-positive partner", 
               "Scale term of inverse logistic distribution\nof age-dependent relative risk\nof HIV-positive partner", 
               "Relative risk of anal intercourse")

#### CALCULATE P_ACC FOR ALL SAMPLES ####
for(dist_n in 0:(n_distributions - 1)) {
  print(paste0("Distribution: ", dist_n))
  print(get(x = paste0("particles_acc_", dist_n))[, sum(t == dist_n)/.N])
}

#### P_ACC = 0.10 ####
# Fewer than 10% of particles were accepted at iteration 13, so so set n_iterations to 14 (13 iterations + iteration 0)
n_distributions <- 14

cols <- hue_pal()(n_distributions + 1) # Number of intermediate distributions + 1 (for prior distribution)

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
  densities_dt_list <- lapply(0:(n_distributions - 1), function(iter) data.table(x = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]$x,
                                                                                 y = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]$y,
                                                                                 dist = paste0("Intermediate distribution ", iter)))
  
  # Combine above two data tables
  densities_dt <- rbindlist(l = densities_dt_list)
  densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))
  
  # Factor variable dist so that colors and legend values are plotted in desired order
  densities_dt[, dist := factor(x = dist, levels = c("Prior", sapply(0:(n_distributions - 1), function(x) paste0("Intermediate distribution ", x))))]
  
  # Create ggplot
  density_plot <- ggplot() +
    geom_line(data = densities_dt, aes(x = x, y = y, group = dist, col = dist), size = 1) +
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

ggarrange(lambda_plot, rr_ai_plot, c_plot, s_plot, cond_rr_plot, ncol = 2, nrow = 3, common.legend = T, legend = "right")

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Manuscript/Supplement/Figures/ABC_Intermediate_Distributions_P_Acc_10.pdf", width = 8.5, height = 7)

#### P_ACC = 0.025 ####
n_distributions <- 24

cols <- hue_pal()(n_distributions + 1) # Number of intermediate distributions + 1 (for prior distribution)

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
  densities_dt_list <- lapply(0:(n_distributions - 1), function(iter) data.table(x = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]$x,
                                                                                 y = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]$y,
                                                                                 dist = paste0("Intermediate distribution ", iter)))
  
  # Combine above two data tables
  densities_dt <- rbindlist(l = densities_dt_list)
  densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))
  
  # Factor variable dist so that colors and legend values are plotted in desired order
  densities_dt[, dist := factor(x = dist, levels = c("Prior", sapply(0:(n_distributions - 1), function(x) paste0("Intermediate distribution ", x))))]
  
  # Create ggplot
  density_plot <- ggplot() +
    geom_line(data = densities_dt, aes(x = x, y = y, group = dist, col = dist), size = 1) +
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
          legend.text = element_text(size = 11, family = "Times")) +
    guides(col = guide_legend(ncol = 1))
  
  assign(x = paste0(params_vec[i], "_plot"), value = density_plot)
  
  rm(list = c("prior_d_dt", "densities_dt_list", "densities_dt", "density_plot"))
}

ggarrange(lambda_plot, rr_ai_plot, c_plot, s_plot, cond_rr_plot, ncol = 2, nrow = 3, common.legend = T, legend = "right")

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Manuscript/Supplement/Figures/ABC_Intermediate_Distributions_P_Acc_2.5.pdf", width = 8.5, height = 7)

#### INDIVIDUAL PARAMETERS ####
x_lab_vec <- c("Baseline probability of HIV transmission", 
               "Relative risk of HIV-positive partner among women reporting condom use at last sex", 
               "Carrying capacity term of inverse logistic distribution of age-dependent relative risk of HIV-positive partner", 
               "Scale term of inverse logistic distribution of age-dependent relative risk of HIV-positive partner", 
               "Relative risk of anal intercourse")

file_name_vec <- c("Lambda", "RR_Condom", "C", "S", "RR_AI")

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
  densities_dt_list <- lapply(0:(n_distributions - 1), function(iter) data.table(x = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]$x,
                                                                                 y = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]$y,
                                                                                 dist = paste0("Intermediate distribution ", iter)))
  
  # Combine above two data tables
  densities_dt <- rbindlist(l = densities_dt_list)
  densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))
  
  # Factor variable dist so that colors and legend values are plotted in desired order
  densities_dt[, dist := factor(x = dist, levels = c("Prior", sapply(0:(n_distributions - 1), function(x) paste0("Intermediate distribution ", x))))]
  
  # Create ggplot
  density_plot <- ggplot() +
    geom_line(data = densities_dt, aes(x = x, y = y, group = dist, col = dist), size = 1) +
    scale_color_manual(values = cols) +
    labs(x = x_lab_vec[i], y = "Density", col = "") +
    theme(axis.text = element_text(size = 11),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 11),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          text = element_text(family = "Times")) +
    facet_wrap(. ~ dist, ncol = 5)

  assign(x = paste0(params_vec[i], "_plot"), value = density_plot)
  
  ggsave(filename = paste0("/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Manuscript/Supplement/Figures/ABC_Intermediate_Distributions_", file_name_vec[i], ".pdf"), width = 8.5, height = 7)
  
  rm(list = c("prior_d_dt", "densities_dt_list", "densities_dt", "density_plot"))
}
