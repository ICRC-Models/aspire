## Purpose: ASPIRE model is currently not working as expected, but is rather producing HRs that do not match expected values when providing a given input for efficacy of the ring. Here, we provide functions to evaluate simulations in which the model has been stripped down to basics, iteratively adding complexity.
## Author:  Kathryn Peebles
## Date:    29 July 2019

# Function to run specified number of simulations
run_debug_sims <- function(n_sims, run_debug_hiv, run_debug_hiv_acts, run_debug_acts, run_debug_hiv_risk, run_debug_equal_risk, run_debug_rerandomize) {
  output_list <- lapply(X   = 1:n_sims,
                        FUN = function(x) { 
                          run_sim_secondary(rr_ai               = params$rr_ai,
                                            prev_ai             = 0,
                                            prop_ai             = 0,
                                            prop_full_adh       = 1,
                                            prop_partial_adh    = 0,
                                            prop_non_adh        = 0,
                                            rr_ring_full_adh    = 0.25,
                                            rr_ring_partial_adh = 0.65,
                                            i                   = x,
                                            output_dir          = "",
                                            debug_sim_hiv       = run_debug_hiv,
                                            debug_sim_hiv_acts  = run_debug_hiv_acts,
                                            debug_sim_acts      = run_debug_acts,
                                            debug_sim_hiv_risk  = run_debug_hiv_risk,
                                            debug_rerandomize      = run_debug_rerandomize) })
  return(output_list)
}

# Function to plot histogram of per-act RRs. Includes all infections (possible to have > 1 per woman) and acts. 
plot_hist_per_act_rr <- function(output_list, n_sims) {
  
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_total)/mean(n_acts_total)]/x$ri_dt[arm == 0, mean(n_inf_total)/mean(n_acts_total)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations"), subtitle = "Includes all infections (>1 per woman possible) and acts") +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Function to plot histogram of per-act RRs. Infections and acts are censored at the first infection per woman.
plot_hist_per_act_rr_one_inf_per_woman <- function(output_list, n_sims) {
  
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_one_inf_per_woman)/mean(n_acts_one_inf_per_woman)]/x$ri_dt[arm == 0, mean(n_inf_one_inf_per_woman)/mean(n_acts_one_inf_per_woman)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations"), subtitle = "Infections and acts censored at first infection") +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Function to plot histogram of per-act RRs. Infections are censored at the first infection per woman, but acts are not (i.e., this is the data we'd expect to get from a trial - one infection observed, but all acts reported since we don't know at which act infection occurred.)
plot_hist_per_act_rr_trial_data <- function(output_list, n_sims) {
  
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_one_inf_per_woman)/mean(n_acts_total)]/x$ri_dt[arm == 0, mean(n_inf_one_inf_per_woman)/mean(n_acts_total)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations"), subtitle = "Includes all acts and one infection per woman") +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Plot to evaluate balance by arm prior to applying reduction in risk of ring. Plot shows per-act risk ratio, calculated as the arm-specific number of infections per time step prior to reduction in risk from ring over arm-specific total acts.
plot_hist_per_act_no_ring <- function(output_list, n_sims) {
  
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_pre_ring)/mean(n_acts_total)]/x$ri_dt[arm == 0, mean(n_inf_pre_ring)/mean(n_acts_total)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, 3)) +
    labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations\nprior to application of ring"), subtitle = "Includes all infections (>1 per woman possible) and acts") +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Function to plot exposure rate ratio over time
plot_exp_rr <- function(output_list, n_sims) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$exposure_dt), idcol = "sim")
  
  min_time_val <- min(output_dt[, max(visit), by = sim]$V1)
  
  # Convert dataset from long to wide
  output_dt_wide <- dcast(data = output_dt, formula = sim + visit ~ arm, value.var = c("exposures", "pre_ring_risk"))
  
  # Create new variable for exposure rate ratio
  output_dt_wide[, `:=`(exp_rr      = exposures_1/exposures_0,
                        cum_risk_rr = pre_ring_risk_1/pre_ring_risk_0)]
  
  # Calculate mean exposure rate ratio at each time step across simulations. Include the mean only for those simulations with the maximum follow-up time observed across all sims.
  mean_dt <- output_dt_wide[visit <= min_time_val, .(exp_rr = mean(exp_rr), cum_risk_rr = mean(cum_risk_rr)), by = visit]
  
  min_y_val <- min(c(output_dt_wide[, min(exp_rr, cum_risk_rr)], 1))
  max_y_val <- output_dt_wide[, max(exp_rr, cum_risk_rr)]
  
  panel_1 <- ggplot(data = output_dt_wide) +
    geom_line(aes(x = visit, y = exp_rr, group = sim), col = "blue", alpha = 0.05) +
    geom_point(data = mean_dt, aes(x = visit, y = exp_rr), col = "blue") + 
    labs(x = "Time (months)", y = "Exposure rate ratio") +
    coord_cartesian(ylim = c(min_y_val, max_y_val)) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_blank())
  
  panel_2 <- ggplot(data = output_dt_wide) +
    geom_line(aes(x = visit, y = cum_risk_rr, group = sim), col = "blue", alpha = 0.05) +
    geom_point(data = mean_dt, aes(x = visit, y = cum_risk_rr), col = "blue") + 
    labs(x = "Time (months)", y = "Cumulative risk rate ratio") +
    coord_cartesian(ylim = c(min_y_val, max_y_val)) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_blank())
  
  plot_together <- ggarrange(panel_1,panel_2, ncol = 2, nrow = 1)
  
  return(annotate_figure(p = plot_together, top = text_grob(paste0("Exposure/risk rate ratio over time (DPV/placebo) across ", n_simulations, " simulations"), family = "Times", size = 14, face = "bold")))
  
}

# Function to plot balance of different quantities affecting HIV risk by arm
plot_balance_factors <- function(output_list, n_sims) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(n_acts)), by = arm]))
  
  panel_1 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(arm)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(arm), col = as.factor(arm)), size = 1, alpha = 0.5) +
    labs(x = "Number of sex acts\n(among those w/ HIV+ partner)", y = "Density") +
    facet_grid(arm ~ ., labeller = as_labeller(c(`0` = "Placebo", `1` = "Dapivirine"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 11),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  tab_1 <- data.table(round(sapply(0:1, function(x) output_dt[arm == x, summary(x_var)]), 0))
  tab_1[, " " := c("Min.", "1st Qu.", "Med.", "Mean", "3rd Qu.", "Max.")]
  setnames(tab_1, old = c("V1", "V2"), new = c("Placebo", "DPV"))
  setcolorder(tab_1, ncol(tab_1))
  tab_1 <- tableGrob(tab_1, theme = ttheme_default(base_size = 9), rows = NULL)
  
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(median_age)), by = arm]))
  
  panel_2 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(arm)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(arm), col = as.factor(arm)), size = 1, alpha = 0.5) +
    labs(x = "Median age\n(among those w/ HIV+ partner)", y = "") +
    facet_grid(arm ~ ., labeller = as_labeller(c(`0` = "Placebo", `1` = "Dapivirine"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 11),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  tab_2 <- data.table(round(sapply(0:1, function(x) output_dt[arm == x, summary(x_var)]), 1))
  tab_2[, " " := c("Min.", "1st Qu.", "Med.", "Mean", "3rd Qu.", "Max.")]
  setnames(tab_2, old = c("V1", "V2"), new = c("Placebo", "DPV"))
  setcolorder(tab_2, ncol(tab_2))
  tab_2 <- tableGrob(tab_2, theme = ttheme_default(base_size = 9), rows = NULL)
  
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(median_vl)), by = arm]))
  
  panel_3 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(arm)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(arm), col = as.factor(arm)), size = 1, alpha = 0.5) +
    labs(x = "Median viral load", y = "Density") +
    facet_grid(arm ~ ., labeller = as_labeller(c(`0` = "Placebo", `1` = "Dapivirine"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 11),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  tab_3 <- data.table(round(sapply(0:1, function(x) output_dt[arm == x, summary(x_var)]), 1))
  tab_3[, " " := c("Min.", "1st Qu.", "Med.", "Mean", "3rd Qu.", "Max.")]
  setnames(tab_3, old = c("V1", "V2"), new = c("Placebo", "DPV"))
  setcolorder(tab_3, ncol(tab_3))
  tab_3 <- tableGrob(tab_3, theme = ttheme_default(base_size = 9), rows = NULL)
  
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(median_prob_condom)), by = arm]))
  
  panel_4 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(arm)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(arm), col = as.factor(arm)), size = 1, alpha = 0.5) +
    labs(x = "Condom use proportion\n(among acts w/ HIV+ partner)", y = "") +
    facet_grid(arm ~ ., labeller = as_labeller(c(`0` = "Placebo", `1` = "Dapivirine"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 11),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  tab_4 <- data.table(round(sapply(0:1, function(x) output_dt[arm == x, summary(x_var)]), 3))
  tab_4[, " " := c("Min.", "1st Qu.", "Med.", "Mean", "3rd Qu.", "Max.")]
  setnames(tab_4, old = c("V1", "V2"), new = c("Placebo", "DPV"))
  setcolorder(tab_4, ncol(tab_4))
  tab_4 <- tableGrob(tab_4, theme = ttheme_default(base_size = 9), rows = NULL)
  
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(mean_sti)), by = arm]))
  
  panel_5 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(arm)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(arm), col = as.factor(arm)), size = 1, alpha = 0.5) +
    labs(x = "Mean number of STIs\n(among those w/ HIV+ partner)", y = "Density") +
    facet_grid(arm ~ ., labeller = as_labeller(c(`0` = "Placebo", `1` = "Dapivirine"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 11),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  tab_5 <- data.table(round(sapply(0:1, function(x) output_dt[arm == x, summary(x_var)]), 3))
  tab_5[, " " := c("Min.", "1st Qu.", "Med.", "Mean", "3rd Qu.", "Max.")]
  setnames(tab_5, old = c("V1", "V2"), new = c("Placebo", "DPV"))
  setcolorder(tab_5, ncol(tab_5))
  tab_5 <- tableGrob(tab_5, theme = ttheme_default(base_size = 9), rows = NULL)
  
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(prop_bv)), by = arm]))
  
  panel_6 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(arm)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(arm), col = as.factor(arm)), size = 1, alpha = 0.5) +
    labs(x = "Proportion of acts with BV\n(among those w/ HIV+ partner)", y = "") +
    facet_grid(arm ~ ., labeller = as_labeller(c(`0` = "Placebo", `1` = "Dapivirine"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 11),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  tab_6 <- data.table(round(sapply(0:1, function(x) output_dt[arm == x, summary(x_var)]), 3))
  tab_6[, " " := c("Min.", "1st Qu.", "Med.", "Mean", "3rd Qu.", "Max.")]
  setnames(tab_6, old = c("V1", "V2"), new = c("Placebo", "DPV"))
  setcolorder(tab_6, ncol(tab_6))
  tab_6 <- tableGrob(tab_6, theme = ttheme_default(base_size = 9), rows = NULL)
  
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(n_partners_hiv)), by = arm]))
  
  panel_7 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(arm)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(arm), col = as.factor(arm)), size = 1, alpha = 0.5) +
    labs(x = "Number of women with at least one\nHIV+ partner and at least one sex act", y = "Density") +
    facet_grid(arm ~ ., labeller = as_labeller(c(`0` = "Placebo", `1` = "Dapivirine"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  tab_7 <- data.table(round(sapply(0:1, function(x) output_dt[arm == x, summary(x_var)]), 0))
  tab_7[, " " := c("Min.", "1st Qu.", "Med.", "Mean", "3rd Qu.", "Max.")]
  setnames(tab_7, old = c("V1", "V2"), new = c("Placebo", "DPV"))
  setcolorder(tab_7, ncol(tab_7))
  tab_7 <- tableGrob(tab_7, theme = ttheme_default(base_size = 9), rows = NULL)
  
  plot_together <- ggarrange(panel_1, tab_1, panel_2, tab_2, panel_3, tab_3, panel_4, tab_4, panel_5, tab_5, panel_6, tab_6, panel_7, tab_7, ncol = 4, nrow = 4)
  
  return(annotate_figure(p = plot_together, top = text_grob(paste0("Comparison of risk factors in placebo and dapivirine arms across ", n_simulations, " simulations"), family = "Times", size = 14, face = "bold")))
}

# Function to calculate median per-act risk by arm
calculate_median_risk_arm <- function(output_list, n_sims) {
  
  dt <- data.table(arm = c(0, 1))
  
  # Median across simulations of mean across timesteps within simulation of median VL across people within timestep
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(median_vl)), by = arm]))
  dt[, vl := round(sapply(0:1, function(x) output_dt[arm == x, median(x_var)]), 1)]
  
  # Median age
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(median_age)), by = arm]))
  dt[, age := round(sapply(0:1, function(x) output_dt[arm == x, median(x_var)]), 1)]
  
  # Median proportion condom use
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(median_prob_condom)), by = arm]))
  dt[, prop_condom := round(sapply(0:1, function(x) output_dt[arm == x, median(x_var)]), 3)]
  
  # Median number of STIs
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(mean_sti)), by = arm]))
  dt[, n_sti := round(sapply(0:1, function(x) output_dt[arm == x, median(x_var)]), 3)]
  
  # Median proportion of BV-positive women
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(prop_bv)), by = arm]))
  dt[, prop_bv := round(sapply(0:1, function(x) output_dt[arm == x, median(x_var)]), 3)]
  
  # Number of acts
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$b_dt[, .(x_var = mean(n_acts)), by = arm]))
  dt[, n_acts := round(sapply(0:1, function(x) output_dt[arm == x, median(x_var)]), 0)]
  
  placebo_risk <- 1 - (1 - params$lambda) ^ exp(log(params$rr_vl) *    (dt[arm == 0, vl] - 4.0) +
                                                log(params$rr_age) *   (dt[arm == 0, age] - 30) +
                                                log(params$rr_condom) * dt[arm == 0, prop_condom] +
                                                log(params$rr_sti) *    dt[arm == 0, n_sti] +
                                                log(params$rr_bv) *     dt[arm == 0, prop_bv])
  
  dpv_risk <- 1 - (1 - params$lambda) ^ exp(log(params$rr_vl) *    (dt[arm == 1, vl] - 4.0) +
                                            log(params$rr_age) *   (dt[arm == 1, age] - 30) +
                                            log(params$rr_condom) * dt[arm == 1, prop_condom] +
                                            log(params$rr_sti) *    dt[arm == 1, n_sti] +
                                            log(params$rr_bv) *     dt[arm == 1, prop_bv])
  
  per_act_rr <- dpv_risk/placebo_risk
  
  cumulative_rr <- (dpv_risk * dt[arm == 1, n_acts])/(placebo_risk * dt[arm == 0, n_acts])
  
  return(list(dpv_risk = dpv_risk, placebo_risk = placebo_risk, per_act_rr = per_act_rr, cumulative_rr = cumulative_rr))
}

# Function to evaluate factors predictive of HIV-positive partner status. Compare their distribution by whether a woman has one or more HIV-positive partners.
check_assignment_pos_partner <- function(output_list, n_sims) {
  # Proportion of women with an HIV-positive partner, stratified by age (</>= 27)
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$hiv_assign_dt[, .(x_var = mean(n_part_hiv > 0)), by = f_age > 26]))
  
  panel_1 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(f_age)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(f_age), col = as.factor(f_age)), size = 1, alpha = 0.5) +
    labs(x = "Proportion of women with at least one HIV+ partner", y = "Density") +
    facet_grid(f_age ~ ., labeller = as_labeller(c(`TRUE` = "Age >= 27", `FALSE` = "Age < 27"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 10),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  # Mean reduction in probability of having an HIV-positive partner given age, stratified by having an HIV-positive partner
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$hiv_assign_dt[, .(x_var = mean(m_hiv_rr)), by = n_part_hiv > 0]))
  
  panel_2 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(n_part_hiv)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(n_part_hiv), col = as.factor(n_part_hiv)), size = 1, alpha = 0.5) +
    labs(x = "Mean reduction in probability of\nhaving an HIV-positive partner given age", y = "Density") +
    facet_grid(n_part_hiv ~ ., labeller = as_labeller(c(`TRUE` = "Has at least one\nHIV-positive partner", `FALSE` = "Has no HIV-\npositive partners"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 10),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  # Proportion of women with an HIV-positive partner, stratified by reported use of a condom in the last week at baseline
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$hiv_assign_dt[, .(x_var = mean(n_part_hiv > 0)), by = b_condom_lweek]))
  
  panel_3 <- ggplot(data = output_dt) +
    geom_histogram(aes(x = x_var, y = ..density.., fill = as.factor(b_condom_lweek)), alpha = 0.75, bins = 30) +
    geom_density(aes(x = x_var, fill = as.factor(b_condom_lweek), col = as.factor(b_condom_lweek)), size = 1, alpha = 0.5) +
    labs(x = "Proportion of women with at least one HIV+ partner", y = "Density") +
    facet_grid(b_condom_lweek ~ ., labeller = as_labeller(c(`0` = "No condom use\nin last week", `1` = "Condom use\nin last week"))) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 10),
          strip.background = element_rect(fill = NA, color = "black"),
          legend.position = "none")
  
  plot_together <- ggarrange(panel_1, panel_2, panel_3, ncol = 3, nrow = 1)
  
  return(annotate_figure(p = plot_together, top = text_grob(paste0("Distribution of factors predictive of HIV-positive partner status across ", n_simulations, " simulations"), family = "Times", size = 14, face = "bold")))
}

# Function to plot risk by arm over time. One panel shows risk calculated by cumulative risk, and the other shows new infections
plot_risk_time <- function(output_list, n_sims) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$ri_dt), idcol = "id")
  
  max_y_val <- output_dt[, max(c(risk, inf))]
  
  min_time_val <- min(output_dt[, max(time), by = id]$V1)
  
  panel_1 <- ggplot(data = output_dt) +
    geom_line(aes(x = time, y = risk, col = as.factor(arm), group = interaction(id, as.factor(arm))), alpha = 0.025) +
    geom_smooth(data = output_dt[time <= min_time_val], aes(x = time, y = risk, col = as.factor(arm)), method = "loess", formula = y ~ x, se = F) +
    geom_smooth(data = output_dt, aes(x = time, y = risk, col = as.factor(arm)), method = "loess", formula = y ~ x, se = F, size = 0.5) +
    labs(x = "Time (months)", y = "Cumulative risk", col = "Arm") +
    scale_color_discrete(labels = c("Placebo", "Dapivirine")) +
    theme(text = element_text(family = "Times")) +
    coord_cartesian(ylim = c(0, max_y_val))
  
  panel_2 <- ggplot(data = output_dt) +
    geom_line(aes(x = time, y = inf, col = as.factor(arm), group = interaction(id, as.factor(arm))), alpha = 0.025) +
    geom_smooth(data = output_dt[time <= min_time_val], aes(x = time, y = risk, col = as.factor(arm)), method = "loess", formula = y ~ x, se = F) +
    geom_smooth(data = output_dt, aes(x = time, y = risk, col = as.factor(arm)), method = "loess", formula = y ~ x, se = F, size = 0.5) +
    labs(x = "Time (months)", y = "New infections", col = "Arm") +
    scale_color_discrete(labels = c("Placebo", "Dapivirine")) +
    theme(text = element_text(family = "Times")) +
    coord_cartesian(ylim = c(0, max_y_val))
  
  plot_together <- ggarrange(panel_1, panel_2, ncol = 2, common.legend = T, legend = "right")
  
  return(annotate_figure(p = plot_together, top = text_grob(paste0("Risk over time across ", n_simulations, " simulations\nBold line is smoothed mean for simulations with the min. follow-up across all simulations"), family = "Times", size = 14, face = "bold")))
}

# Function to plot number of active partners over time
plot_n_active_partners <- function(output_list, n_sims) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$n_active_part)), idcol = "id")
  output_dt[, month := 1:.N, by = id]
  
  min_time_val <- min(output_dt[, max(month), by = id]$V1)
  
  ret_plot <- ggplot(data = output_dt) +
    geom_line(aes(x = month, y = V1, group = id), alpha = 0.025, color = "blue") +
    geom_smooth(data = output_dt[month > 1 & month <= min_time_val], aes(x = month, y = V1), method = "loess", formula = y ~ x, se = F) +
    geom_smooth(data = output_dt[month > 1], aes(x = month, y = V1), method = "loess", formula = y ~ x, se = F, linetype = "dashed") +
    geom_hline(yintercept = 3173, col = "red") +
    labs(x = "Time (months)", y = "Number of active relationships") +
    theme(text = element_text(family = "Times")) +
    coord_cartesian(ylim = c(3000, 3500))
  
  return(ret_plot)
}

# Function to plot number of active partners over time by arm
plot_n_active_partners_arm <- function(output_list, n_sims) {
  output_dt_dpv <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$n_active_part_dpv)), idcol = "id")
  output_dt_dpv[, month := 1:.N, by = id]
  output_dt_dpv[, arm := "Dapivirine"]
  
  output_dt_pbo <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$n_active_part_pbo)), idcol = "id")
  output_dt_pbo[, month := 1:.N, by = id]
  output_dt_pbo[, arm := "Placebo"]
  
  output_dt <- rbindlist(l = list(output_dt_dpv, output_dt_pbo))
  
  min_time_val <- min(output_dt[, max(month), by = id]$V1)
  
  ret_plot <- ggplot(data = output_dt) +
    geom_line(aes(x = month, y = V1, group = interaction(id, arm), color = arm), alpha = 0.025) +
    geom_smooth(data = output_dt[month > 1 & month <= min_time_val], aes(x = month, y = V1, color = arm), method = "loess", formula = y ~ x, se = F) +
    geom_smooth(data = output_dt[month > 1], aes(x = month, y = V1, color = arm), method = "loess", formula = y ~ x, se = F, linetype = "dashed") +
    geom_hline(yintercept = 3173, col = "red") +
    labs(x = "Time (months)", y = "Number of active relationships", color = "Arm") +
    theme(text = element_text(family = "Times"),
          legend.key = element_blank()) +
    coord_cartesian(ylim = c(0, 1800))
  
  return(ret_plot)
}

# Function to plot female age against male age
plot_ages <- function(output_list, n_sims) {
  
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$male_dt)), idcol = "sim_id")
  output_dt[, f_age_cat := c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49")[findInterval(x = f_age, vec = c(15, 20, 25, 30, 35, 40, 45, 50))]]
  
  ret_plot <- ggplot(data = output_dt) +
    geom_boxplot(aes(x = m_age, y = f_age)) +
    labs(x = "Male partner age", y = "Female age") +
    theme(text = element_text(family = "Times"))

  return(ret_plot)
}

# Function to plot histogram of HRs stratified by site
plot_hist_hr <- function(output_list, n_sims) {
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(hr = sapply(output_list, function(x) x$hr)), aes(x = hr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "Hazard ratio", y = "Count", title = paste0("Distribution of hazard ratio across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Function to plot histogram of HRs stratified by site, BV, and number of sex acts >/< median number of monthly sex acts
plot_hist_hr_strat <- function(output_list, n_sims) {
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(hr = sapply(output_list, function(x) x$hr_strat)), aes(x = hr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = "Hazard ratio", y = "Count", title = paste0("Distribution of hazard ratio across ", n_sims, " simulations"), subtitle = "HR stratified by site, BV status, and binary number of average monthly sex acts") +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Function to create histogram of trial duration
plot_hist_trial_duration <- function(output_list, n_sims) {
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(max_visit = sapply(output_list, function(x) x$study_dt[, max(visit)])), aes(x = max_visit), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, max(sapply(output_list, function(x) x$study_dt[, max(visit)])) + 1)) +
    labs(x = "Trial duration (months)", y = "Count", title = paste0("Distribution of trial duration across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Function to create histogram of final number of infections
plot_hist_inf <- function(output_list, n_sims) {
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(n_inf = sapply(output_list, function(x) x$study_dt[, sum(hiv)])), aes(x = n_inf), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, 200)) +
    labs(x = "Number of HIV infections", y = "Count", title = paste0("Distribution of final number of infections across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Function to create summary table of per-act relative risk values. All infections and acts are included.
create_table_per_act_rr <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), RR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$ri_dt[arm == 1, mean(n_inf_total)/mean(n_acts_total)]/x$ri_dt[arm == 0, mean(n_inf_total)/mean(n_acts_total)], 2)))))))
}

# Function to create summary table of per-act relative risk values. Infections and acts are censored at first infection per woman.
create_table_per_act_rr_one_inf_per_woman <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), RR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$ri_dt[arm == 1, mean(n_inf_one_inf_per_woman)/mean(n_acts_one_inf_per_woman)]/x$ri_dt[arm == 0, mean(n_inf_one_inf_per_woman)/mean(n_acts_one_inf_per_woman)], 2)))))))
}

# Function to create summary table of per-act relative risk values. Infections are censored at first infection per woman, but all acts are included (similar to data that would be collected in a clinical trial).
create_table_per_act_rr_trial_data <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), RR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$ri_dt[arm == 1, mean(n_inf_one_inf_per_woman)/mean(n_acts_total)]/x$ri_dt[arm == 0, mean(n_inf_one_inf_per_woman)/mean(n_acts_total)], 2)))))))
}

# Function to create summary table of per-act relative risk values prior to application of ring RR.  All infections and acts are included.
create_table_per_act_rr_no_ring <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), RR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$ri_dt[arm == 1, mean(n_inf_pre_ring)/mean(n_acts_total)]/x$ri_dt[arm == 0, mean(n_inf_pre_ring)/mean(n_acts_total)], 2)))))))
}

# Function to create summary table of hazard ratios
create_table_hr <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), HR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$hr, 2)))))))
}

# Function to create summary table of hazard ratios, stratified by site, BV status, and number of monthly sex acts >/< median numbere of sex acts
create_table_hr_strat <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), HR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$hr_strat, 2)))))))
}

# Function to create summary table of trial duration
create_table_dur <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), Duration = as.numeric(unname(summary(sapply(output_list, function(x) x$study_dt[, max(visit)]))))))
}

# Function to create summary table of final number of infections
create_table_inf <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), `N infections` = as.numeric(unname(summary(sapply(output_list, function(x) x$study_dt[, sum(hiv)]))))))
}
