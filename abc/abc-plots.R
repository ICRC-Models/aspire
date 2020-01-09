## Purpose: Create plot that shows distribution of simulated particles
## Author:  Kathryn Peebles
## Date:    7 January 2019

library(ggplot2)
library(ggtern)
library(data.table)
library(ggpubr)
library(scales)
library(viridis)

setwd("~/Documents/code/aspire/abc")

# Plot line showing density of particles for each iteration
load(paste0(getwd(), "/priors.RDATA"))

n_distributions <- 13 # number of distributions + 1 for iteration 0

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
col_func <- colorRampPalette(colors = c("#D5FD83", "#5C6C39"))
# cols <- col_func(n_distributions + 1)
# cols <- hue_pal()(n_distributions + 1) # Number of intermediate distributions + 1 (for prior distribution)
cols <- rev(plasma(n_distributions + 1))

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
  densities_dt_list <- lapply(0:(n_distributions - 1), function(iter) data.table(x = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w, adjust = 2), .SDcols = c(params_vec[i], "w")]$x,
                                                                                 y = get(x = paste0("particles_acc_", iter))[, density(x = unname(unlist(.SD[, 1])), weights = w, adjust = 2), .SDcols = c(params_vec[i], "w")]$y,
                                                                                 dist = paste0("Intermediate distribution ", iter)))
  
  # Combine above two data tables
  densities_dt <- rbindlist(l = densities_dt_list)
  densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))
  
  # Factor variable dist so that colors and legend values are plotted in desired order
  densities_dt[, dist := factor(x = dist, levels = c("Prior", sapply(0:(n_distributions - 1), function(x) paste0("Intermediate distribution ", x))))]
  
  densities_dt[dist == paste0("Intermediate distribution ", (n_distributions - 1)), prob_y := y/(sum(y))]
  
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

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Manuscript supplement/Figures/ABC_Interim_Distributions_20191111.pdf", width = 8.5, height = 7)

#### BIVARIATE DENSITY PLOTS ####
# Plot relative increase in partner change rate and reference male HIV-1 incidence together
hm <- ggtern::kde2d.weighted(x = particles_acc_12[, p_rate_rr], y = particles_acc_12[, base_male_hiv_incidence], n = 100, w = particles_acc_12[, w])

hm_max <- hm$z == max(hm$z)
hm_row <- which.max(rowSums(hm_max))
hm_col <- which.max(colSums(hm_max))
bestfit <- c(hm$x[hm_row], hm$y[hm_col])
bestfit

pdf(file = "/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Manuscript supplement/Figures/ABC_Bivariate_Distribution_20191110.pdf", width = 8.5, height = 4)
image(x    = hm, 
      col  = plasma(n = 1000),
      xlab = "Relative increase in partner change rate", 
      ylab = "Reference male HIV-1 incidence",
      family = "serif")
points(x = bestfit[1], y = bestfit[2], col = "black", pch = 19)
dev.off()
  
# For dissertation defense
pdf(file = "/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/ABC_Bivariate_Distribution_20191110.pdf", width = 4, height = 4)
image(x    = hm, 
      col  = plasma(n = 1000),
      xlab = "Relative increase in partner change rate", 
      ylab = "Reference male HIV-1 incidence",
      family = "sans")
points(x = bestfit[1], y = bestfit[2], col = "black", pch = 19)
dev.off()

#### IDENTIFY BEST-FIT PARAMETERS ####
# Load particles for iteration with p_acc = 0.100
load(paste0(getwd(), "/hyak/t12/particles_acc_12.RDATA"))

# Identify best-fit parameter from univariate weighted densities for each parameter
for(i in 1:length(params_vec)) {
  d <- particles_acc_11[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]
  
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

ggplot(data = densities_dt) +
  geom_line(aes(x = x, y = y)) +
  geom_vline(xintercept = cri, linetype = "dashed")

#### INCIDENCE PLOTS FOR CALIBRATION SLIDES ####
n_iter <- 7 # 4 iterations plus prior and particles accepted from the prior

for(i in 1:n_iter) {
  if(i == 1 ) {
    i_to_sample <- sample(x = 1:100000, size = 100, replace = F)
    
    # Create empty data table in which to store incidence
    inc_dt <- data.table(f_age_cat = NA_real_, hiv_inf = NA_real_, n = NA_real_, py = NA_real_, i = NA_real_, panel = NA_real_)
    
    for(i in i_to_sample) {
      load(paste0("~/Documents/code/aspire/abc/hyak/t0/inf_sims/inf_sim_", i, ".RData"))
      inf_sim[, i := i]
      inf_sim[1:3, panel := "age"]
      inf_sim[4:5, panel := "time"]
      inf_sim[4, f_age_cat := "Months 15-23"]
      inf_sim[5, f_age_cat := "Months 24-31"]
      
      inc_dt <- rbindlist(l = list(inc_dt, inf_sim))
    }
    
    setnames(inc_dt, old = "f_age_cat", new = "age_cat")
    inc_dt <- inc_dt[!is.na(hiv_inf)]
    inc_dt[, inc_rate := (hiv_inf/py) * 100]
    
    assign(x = paste0("inc_prior"), value = inc_dt)
    rm(inc_dt)
  } else {
    rows_to_sample <- sample(x = 1:60000, size = 100, replace = F)
    temp_particles <- eval(parse(text = paste0("particles_acc_", i - 2)))
    particles_acc_sub <- temp_particles[rows_to_sample, .(i, t)]
    
    # Create empty data table in which to store incidence
    inc_dt <- data.table(f_age_cat = NA_real_, hiv_inf = NA_real_, n = NA_real_, py = NA_real_, i = NA_real_, panel = NA_real_)
    
    for(j in 1:nrow(particles_acc_sub)) {
      load(paste0("~/Documents/code/aspire/abc/hyak/t", particles_acc_sub[j, t], "/inf_sims/inf_sim_", particles_acc_sub[j, i], ".RData"))
      inf_sim[, i := particles_acc_sub[j, i]]
      inf_sim[1:3, panel := "age"]
      inf_sim[4:5, panel := "time"]
      inf_sim[4, f_age_cat := "Months 15-23"]
      inf_sim[5, f_age_cat := "Months 24-31"]
      
      inc_dt <- rbindlist(l = list(inc_dt, inf_sim))
    }
    
    setnames(inc_dt, old = "f_age_cat", new = "age_cat")
    inc_dt <- inc_dt[!is.na(hiv_inf)]
    inc_dt[, inc_rate := (hiv_inf/py) * 100]
    
    assign(x = paste0("inc_", i - 2), value = inc_dt)
    rm(inc_dt)
  }
}

load("~/Documents/code/aspire/data-public/inf_obs.RData")

inf_obs[1:3, panel := "age"]
inf_obs[4:5, panel := "time"]
inf_obs[4, age_cat := "Months 15-23"]
inf_obs[5, age_cat := "Months 24-31"]

# Colors
# Dark purple: "#5F2285"
# Medium purple: "#8C30C4"
# Light purple: "#AF3AF5"
# Dark green: "#5C6C39"
# Light green: "#D5FD83"
# Alt green: "#B2D25C"

max_y_val <- max(c(inc_0[, inc_rate], inc_1[, inc_rate], inc_2[, inc_rate], inc_3[, inc_rate], inc_4[, inc_rate], inc_prior[, inc_rate]))

for(i in 1:n_iter) {
  if(i == 1) {
    inc_dt <- eval(parse(text = "inc_prior"))
  } else {
    inc_dt <- eval(parse(text = paste0("inc_", i - 2)))
  }
  
  ggplot(data = inc_dt) +
    geom_jitter(data = inc_dt, aes(x = age_cat, y = inc_rate), col = "#AF3AF5", alpha = 0.5, width = 0.1) +
    geom_point(data = inf_obs, aes(x = age_cat, y = inc_rate), col = "black", size = 2.25) +
    geom_errorbar(data = inf_obs, aes(x = age_cat, ymin = inc_rate - 1.96 * se, ymax = inc_rate + 1.96 * se), width = 0.1, size = 0.75) +
    coord_cartesian(ylim = c(0, 20)) +
    labs(x = "", y = "Incidence rate (per 100 woman-years)", col = "") +
    theme(text = element_text(family = "Times", size = 17),
          panel.border = element_rect(color = "black", fill = NA),
          panel.background = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.05, colour = "black"),
          legend.key = element_blank(),
          strip.background = element_blank()) +
    facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = as_labeller(c("age" = "Age group", "time" = "Trial time period")))
  
  ggsave(filename = paste0("~/Documents/code/aspire/abc/figures/incidence_", i, ".pdf"), width = 10, height = 5)
  
}

#### INCIDENCE PLOTS FOR SUPPLEMENT ####
# Sample row numbers to subset from final set of accepted particles
rows_to_sample <- sample(x = 1:60000, size = 100, replace = F)

# Subset final set of accepted particles for plotting
particles_acc_sub <- particles_acc_4[rows_to_sample, .(i, t)]

# Create empty data table in which to store incidence
inc_dt <- data.table(f_age_cat = NA_real_, hiv_inf = NA_real_, n = NA_real_, py = NA_real_, i = NA_real_, t = NA_real_, panel = NA_real_)

for(j in 1:nrow(particles_acc_sub)) {
  load(paste0("~/Documents/code/aspire/abc/hyak/t", particles_acc_sub[j, t], "/inf_sims/inf_sim_", particles_acc_sub[j, i], ".RData"))
  inf_sim[, `:=`(i = particles_acc_sub[j, i],
                 t = particles_acc_sub[j, t])]
  inf_sim[1:3, panel := "age"]
  inf_sim[4:5, panel := "time"]
  inf_sim[4, f_age_cat := "Months 15-23"]
  inf_sim[5, f_age_cat := "Months 24-31"]
  
  inc_dt <- rbindlist(l = list(inc_dt, inf_sim))
}

setnames(inc_dt, old = "f_age_cat", new = "age_cat")
inc_dt <- inc_dt[!is.na(hiv_inf)]
inc_dt[, inc_rate := (hiv_inf/py) * 100]

load("~/Documents/code/aspire/data-public/inf_obs.RData")

inf_obs[1:3, panel := "age"]
inf_obs[4:5, panel := "time"]
inf_obs[4, age_cat := "Months 15-23"]
inf_obs[5, age_cat := "Months 24-31"]

# Colors
# Dark purple: "#5F2285"
# Medium purple: "#8C30C4"
# Light purple: "#AF3AF5"
# Dark green: "#5C6C39"
# Light green: "#D5FD83"
# Alt green: "#B2D25C"

ggplot(data = inc_dt) +
  geom_jitter(data = inc_dt, aes(x = age_cat, y = inc_rate), col = "#AF3AF5", alpha = 0.5, width = 0.1) +
  geom_point(data = inf_obs, aes(x = age_cat, y = inc_rate), col = "black", size = 2.25) +
  geom_errorbar(data = inf_obs, aes(x = age_cat, ymin = inc_rate - 1.96 * se, ymax = inc_rate + 1.96 * se), width = 0.1, size = 0.75) +
  #coord_cartesian(ylim = c(0, 9)) +
  labs(x = "", y = "Incidence rate (per 100 woman-years)", col = "") +
  theme(text = element_text(family = "Times", size = 15),
        panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.05, colour = "black"),
        legend.key = element_blank(),
        strip.background = element_blank()) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = as_labeller(c("age" = "Age group", "time" = "Trial time period")))

ggsave(filename = "~/Documents/code/aspire/abc/figures/incidence_supplement.pdf", width = 8, height = 4)


