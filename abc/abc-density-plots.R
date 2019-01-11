## Purpose: Create plot that shows distribution of simulated particles
## Author:  Kathryn Peebles
## Date:    7 January 2019

library(ggplot2)
library(ggtern)
library(data.table)
library(ggpubr)

setwd("~/Documents/code/aspire/abc")

# Plot line showing density of particles for each iteration
load(paste0(getwd(), "/priors.RDATA"))

n_distributions <- 12

for(i in 0:(n_distribution - 1)) { # n_distributions - 1 to start with 0
  load(paste0(getwd(), "/hyak/t", i, "/particles_acc_", i, ".RDATA"))
}

cols <- hue_pal()(n_distributions + 1) # Number of intermediate distributions + 1 (for prior distribution)

#### LAMBDA ####
# Create data table with density of prior distribution
prior_d_dt <- data.table(x = seq(particles_acc_0[, min(lambda)], particles_acc_0[, max(lambda)], length.out = 60000),
                         y = dbeta(x = seq(particles_acc_0[, min(lambda)], particles_acc_0[, max(lambda)], length.out = 60000), shape1 = priors[param == "lambda", alpha], shape2 = priors[param == "lambda", beta]),
                         s = "Prior")

# Create data table with weighted density of each intermediate distribution
densities_dt_list <- lapply(0:11, function(i) data.table(x = get(x = paste0("particles_acc_", i))[, density(x = lambda, weights = w)]$x,
                                                         y = get(x = paste0("particles_acc_", i))[, density(x = lambda, weights = w)]$y,
                                                         s = paste0("Intermediate distribution ", i)))

# Combine above two data tables
densities_dt <- rbindlist(l = densities_dt_list)
densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))

# Factor variable s so that colors and legend values are plotted in desired order
densities_dt[, s := factor(x = s, levels = c("Prior", sapply(0:11, function(x) paste0("Intermediate distribution ", x))))]

# Create ggplot
lambda_plot <- ggplot() +
  geom_line(data = densities_dt, aes(x = x, y = y, group = s, col = s), size = 1) +
  scale_color_manual(values = cols) +
  labs(x = "Baseline probability of HIV transmission",
       y = "Density",
       col = "") +
  theme(axis.text = element_text(size = 11, family = "Times"),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 11, family = "Times"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 11, family = "Times"))

#### COND_RR ####
# Create data table with density of prior distribution
prior_d_dt <- data.table(x = seq(particles_acc_0[, min(cond_rr)], particles_acc_0[, max(cond_rr)], length.out = 60000),
                         y = dunif(x = seq(particles_acc_0[, min(cond_rr)], particles_acc_0[, max(cond_rr)], length.out = 60000), min = priors[param == "cond_rr", min], max = priors[param == "cond_rr", max]),
                         s = "Prior")

# Create data table with weighted density of each intermediate distribution
densities_dt_list <- lapply(0:11, function(i) data.table(x = get(x = paste0("particles_acc_", i))[, density(x = cond_rr, weights = w)]$x,
                                                         y = get(x = paste0("particles_acc_", i))[, density(x = cond_rr, weights = w)]$y,
                                                         s = paste0("Intermediate distribution ", i)))

# Combine above two data tables
densities_dt <- rbindlist(l = densities_dt_list)
densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))

# Factor variable s so that colors and legend values are plotted in desired order
densities_dt[, s := factor(x = s, levels = c("Prior", sapply(0:11, function(x) paste0("Intermediate distribution ", x))))]

# Create ggplot
cond_rr_plot <- ggplot() +
  geom_line(data = densities_dt, aes(x = x, y = y, group = s, col = s), size = 1) +
  scale_color_manual(values = cols) +
  labs(x = "Relative risk of HIV-positive partner\namong women reporting condom use at last sex",
       y = "Density",
       col = "") +
  theme(axis.text = element_text(size = 11, family = "Times"),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 11, family = "Times"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 11, family = "Times"))

#### C ####
# Create data table with density of prior distribution
prior_d_dt <- data.table(x = seq(particles_acc_0[, min(c)], particles_acc_0[, max(c)], length.out = 60000),
                         y = dunif(x = seq(particles_acc_0[, min(c)], particles_acc_0[, max(c)], length.out = 60000), min = priors[param == "c", min], max = priors[param == "c", max]),
                         s = "Prior")

# Create data table with weighted density of each intermediate distribution
densities_dt_list <- lapply(0:11, function(i) data.table(x = get(x = paste0("particles_acc_", i))[, density(x = c, weights = w)]$x,
                                                         y = get(x = paste0("particles_acc_", i))[, density(x = c, weights = w)]$y,
                                                         s = paste0("Intermediate distribution ", i)))

# Combine above two data tables
densities_dt <- rbindlist(l = densities_dt_list)
densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))

# Factor variable s so that colors and legend values are plotted in desired order
densities_dt[, s := factor(x = s, levels = c("Prior", sapply(0:11, function(x) paste0("Intermediate distribution ", x))))]

# Create ggplot
c_plot <- ggplot() +
  geom_line(data = densities_dt, aes(x = x, y = y, group = s, col = s), size = 1) +
  scale_color_manual(values = cols) +
  labs(x = "Carrying capacity term of inverse logistic\ndistribution of age-dependent relative risk\nof HIV-positive partner",
       y = "Density",
       col = "") +
  theme(axis.text = element_text(size = 11, family = "Times"),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 11, family = "Times"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 11, family = "Times"))

#### S ####
# Create data table with density of prior distribution
prior_d_dt <- data.table(x = seq(particles_acc_0[, min(s)], particles_acc_0[, max(s)], length.out = 60000),
                         y = dunif(x = seq(particles_acc_0[, min(s)], particles_acc_0[, max(s)], length.out = 60000), min = priors[param == "s", min], max = priors[param == "s", max]),
                         s = "Prior")

# Create data table with weighted density of each intermediate distribution
densities_dt_list <- lapply(0:11, function(i) data.table(x = get(x = paste0("particles_acc_", i))[, density(x = s, weights = w)]$x,
                                                         y = get(x = paste0("particles_acc_", i))[, density(x = s, weights = w)]$y,
                                                         s = paste0("Intermediate distribution ", i)))

# Combine above two data tables
densities_dt <- rbindlist(l = densities_dt_list)
densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))

# Factor variable s so that colors and legend values are plotted in desired order
densities_dt[, s := factor(x = s, levels = c("Prior", sapply(0:11, function(x) paste0("Intermediate distribution ", x))))]

# Create ggplot
s_plot <- ggplot() +
  geom_line(data = densities_dt, aes(x = x, y = y, group = s, col = s), size = 1) +
  scale_color_manual(values = cols) +
  labs(x = "Scale term of inverse logistic distribution\nof age-dependent relative risk\nof HIV-positive partner",
       y = "Density",
       col = "") +
  theme(axis.text = element_text(size = 11, family = "Times"),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 11, family = "Times"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 11, family = "Times"))

#### RR_AI ####
# Create data table with density of prior distribution
prior_d_dt <- data.table(x = seq(particles_acc_0[, min(rr_ai)], particles_acc_0[, max(rr_ai)], length.out = 60000),
                         y = dunif(x = seq(particles_acc_0[, min(rr_ai)], particles_acc_0[, max(rr_ai)], length.out = 60000), min = priors[param == "rr_ai", min], max = priors[param == "rr_ai", max]),
                         s = "Prior")

# Create data table with weighted density of each intermediate distribution
densities_dt_list <- lapply(0:11, function(i) data.table(x = get(x = paste0("particles_acc_", i))[, density(x = rr_ai, weights = w)]$x,
                                                         y = get(x = paste0("particles_acc_", i))[, density(x = rr_ai, weights = w)]$y,
                                                         s = paste0("Intermediate distribution ", i)))

# Combine above two data tables
densities_dt <- rbindlist(l = densities_dt_list)
densities_dt <- rbindlist(l = list(prior_d_dt, densities_dt))

# Factor variable s so that colors and legend values are plotted in desired order
densities_dt[, s := factor(x = s, levels = c("Prior", sapply(0:11, function(x) paste0("Intermediate distribution ", x))))]

# Create ggplot
rr_ai_plot <- ggplot() +
  geom_line(data = densities_dt, aes(x = x, y = y, group = s, col = s), size = 1) +
  scale_color_manual(values = cols) +
  labs(x = "Relative risk of anal intercourse",
       y = "Density",
       col = "") +
  theme(axis.text = element_text(size = 11, family = "Times"),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 11, family = "Times"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 11, family = "Times"))

#### ALL PLOTS ####
ggarrange(lambda_plot, rr_ai_plot, c_plot, s_plot, cond_rr_plot, ncol = 2, nrow = 3, common.legend = T, legend = "right")

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Manuscript efficacy dilution/Supplement/Figures/ABC_Intermediate_Distributions.pdf", width = 8.5, height = 7)
