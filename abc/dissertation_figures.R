## Purpose: Run ASPIRE simulations with best-fitting parameter sets from t6 and final set of accepted particles.
## Author:  Kathryn Peebles
## Date:    25 November 2019

library(data.table)
library(ggplot2)
library(viridis)

# Source set-up file
source("~/Documents/code/aspire/set-up.R")

# Load particles accepted in iterations 6 to identify their best-fit parameters
load("~/Documents/code/aspire/abc/hyak/t6/particles_acc_6.RDATA")

params_vec <- c("lambda", "cond_rr", "c", "s", "rr_ai", "p_rate_rr", "base_male_hiv_incidence")
x_lab_vec <- c("Baseline probability of HIV transmission", 
               "Relative risk of HIV-positive partner\namong women reporting condom use at last sex", 
               "Carrying capacity term of inverse logistic\ndistribution of age-dependent relative risk\nof HIV-positive partner", 
               "Scale term of inverse logistic distribution\nof age-dependent relative risk\nof HIV-positive partner", 
               "Relative risk of anal intercourse",
               "Relative increase in partner change rate",
               "Reference male HIV-1 incidence")

# Identify best-fit parameter from univariate weighted densities for each parameter
for(i in 1:length(params_vec)) {
  d <- particles_acc_6[, density(x = unname(unlist(.SD[, 1])), weights = w), .SDcols = c(params_vec[i], "w")]
  
  val <- d$x[which.max(d$y)]
  
  print(paste0("Best-fitting value for ", x_lab_vec[i], " is: ", val))
}

hm <- ggtern::kde2d.weighted(x = particles_acc_6[, p_rate_rr], y = particles_acc_6[, base_male_hiv_incidence], n = 100, w = particles_acc_6[, w])

hm_max <- hm$z == max(hm$z)
hm_row <- which.max(rowSums(hm_max))
hm_col <- which.max(colSums(hm_max))
bestfit <- c(hm$x[hm_row], hm$y[hm_col])
bestfit

params_6 <- data.table(param = params_vec, value = c(0.002, 1.18, 0.26, 0.08, 6.2, 3.9, 0.077))

# Run set of 50 simulations with best-fitting values for iteration 6
results_list <- lapply(X    = 1:50, FUN = function(x) { run_sim(lambda   = params_6[param == "lambda", value],
                                                                cond_rr  = params_6[param == "cond_rr", value],
                                                                c        = params_6[param == "c", value],
                                                                s        = params_6[param == "s", value],
                                                                rr_ai    = params_6[param == "rr_ai", value],
                                                                p_rate_rr = params_6[param == "p_rate_rr", value],
                                                                base_male_hiv_incidence = params_6[param == "base_male_hiv_incidence", value],
                                                                prev_ai  = "prop_ai_obs",
                                                                prop_ai  = 0.078,
                                                                prop_adh = "prop_adh_obs",
                                                                rr_ring  = 0.3,
                                                                re_randomize = F,
                                                                i = NULL,
                                                                calibrate = F) })

# Convert to data table of incidence
inc_6 <- rbindlist(l = lapply(results_list, function(x) x$inc_sim[1:3]))
setnames(inc_6, old = "f_age_cat", new = "age_cat")
inc_6[, inc_rate := (hiv_inf/py) * 100]

save(inc_6, file = "~/Desktop/inc_6.RDATA")

rm(list = ls()[!(ls() == "inc_6")])

source("~/Documents/code/aspire/set-up.R")

# Run set of 50 simulations with best-fitting values for final iteration
results_list <- mclapply(X = 1:50, FUN = function(x) { run_sim(lambda    = params$lambda,
                                                                cond_rr  = params$cond_rr,
                                                                c        = params$c,
                                                                s        = params$s,
                                                                rr_ai    = params$rr_ai,
                                                                p_rate_rr = params$p_rate_rr,
                                                                base_male_hiv_incidence = params$base_male_hiv_incidence,
                                                                prev_ai  = "prop_ai_obs",
                                                                prop_ai  = 0.078,
                                                                prop_adh = "prop_adh_obs",
                                                                rr_ring  = 0.3,
                                                                re_randomize = F,
                                                                i = NULL,
                                                                calibrate = F) }, mc.cores = 4)

# Convert to data table of incidence
inc_12 <- rbindlist(l = lapply(results_list, function(x) x$inc_sim[1:3]))
setnames(inc_12, old = "f_age_cat", new = "age_cat")
inc_12[, inc_rate := (hiv_inf/py) * 100]

save(inc_12, file = "~/Desktop/inc_12.RDATA")

# Load data from complete calibration iterations
n_iter <- c(0:5, 7:11)

# Load tables of accepted particles
for(i in 0:12) { # n_distributions - 1 to start with 0
  load(paste0("~/Documents/code/aspire/abc/hyak/t", i, "/particles_acc_", i, ".RDATA"))
}

for(i in n_iter) {
  if(i == 0) {
    i_to_sample <- sample(x = 1:100000, size = 50, replace = F)
    
    # Create empty data table in which to store incidence
    inc_dt <- data.table(f_age_cat = NA_real_, hiv_inf = NA_real_, n = NA_real_, py = NA_real_, i = NA_real_, panel = NA_real_)
    
    for(j in i_to_sample) {
      load(paste0("~/Documents/code/aspire/abc/hyak/t0/inf_sims/inf_sim_", j, ".RData"))
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
    
    inc_dt <- inc_dt[!(age_cat %in% c("Months 15-23", "Months 24-31"))]
    
    assign(x = paste0("inc_prior"), value = inc_dt)
    rm(inc_dt)
    
    rows_to_sample <- sample(x = 1:60000, size = 50, replace = F)
    temp_particles <- eval(parse(text = paste0("particles_acc_", i)))
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
    
    inc_dt <- inc_dt[!(age_cat %in% c("Months 15-23", "Months 24-31"))]
    
    assign(x = paste0("inc_", i), value = inc_dt)
    rm(inc_dt)
  } else {
    temp_particles <- eval(parse(text = paste0("particles_acc_", i)))
    temp_particles <- temp_particles[t != 6]
    rows_to_sample <- sample(x = 1:nrow(temp_particles), size = 50, replace = F)
    
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
    
    inc_dt <- inc_dt[!(age_cat %in% c("Months 15-23", "Months 24-31"))]
    
    assign(x = paste0("inc_", i), value = inc_dt)
    rm(inc_dt)
  }
}

load("~/Documents/code/aspire/data-public/inf_obs.RData")

inf_obs <- inf_obs[1:3]

# Colors
# Dark purple: "#5F2285"
# Medium purple: "#8C30C4"
# Light purple: "#AF3AF5"
# Dark green: "#5C6C39"
# Light green: "#D5FD83"
# Alt green: "#B2D25C"

max_y_val <- max(c(inc_0[, inc_rate], inc_1[, inc_rate], inc_2[, inc_rate], inc_3[, inc_rate], inc_4[, inc_rate], inc_prior[, inc_rate]))

n_iter <- 14

cols <- viridis(25)
cols <- rev(cols[7:20])

for(i in 1:n_iter) {
  if(i == 1) {
    inc_dt <- eval(parse(text = "inc_prior"))
  } else {
    inc_dt <- eval(parse(text = paste0("inc_", i - 2)))
  }
  
  ggplot(data = inc_dt) +
    geom_jitter(data = inc_dt, aes(x = age_cat, y = inc_rate), col = cols[i], alpha = 0.5, width = 0.1) +
    geom_point(data = inf_obs, aes(x = age_cat, y = inc_rate), col = "black", size = 2.25) +
    geom_errorbar(data = inf_obs, aes(x = age_cat, ymin = inc_rate - 1.96 * se, ymax = inc_rate + 1.96 * se), width = 0.1, size = 0.75) +
    coord_cartesian(ylim = c(0, 20)) +
    labs(x = "Age", y = "Incidence rate (per 100 woman-years)", col = "") +
    theme(text = element_text(family = "sans", size = 15),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 15),
          panel.background = element_blank(),
          legend.key = element_blank(),
          strip.background = element_blank())
  
  ggsave(filename = paste0("/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/Calib_Incidence_", i, ".pdf"), width = 5, height = 5)
}

# Plot with observed data only
ggplot() +
  geom_point(data = inf_obs, aes(x = age_cat, y = inc_rate), col = "black", size = 2.25) +
  geom_errorbar(data = inf_obs, aes(x = age_cat, ymin = inc_rate - 1.96 * se, ymax = inc_rate + 1.96 * se), width = 0.1, size = 0.75) +
  coord_cartesian(ylim = c(0, 20)) +
  labs(x = "Age", y = "Incidence rate (per 100 woman-years)", col = "") +
  theme(text = element_text(family = "sans", size = 15),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 15),
        panel.background = element_blank(),
        legend.key = element_blank(),
        strip.background = element_blank())

ggsave(filename = paste0("/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/Calib_Incidence_Obs.pdf"), width = 5, height = 5)
