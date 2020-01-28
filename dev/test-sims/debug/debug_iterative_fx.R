## Purpose: ASPIRE model is currently not working as expected, but is rather producing HRs that do not match expected values when providing a given input for efficacy of the ring. Here, we provide functions to evaluate simulations in which the model has been stripped down to basics, iteratively adding complexity.
## Author:  Kathryn Peebles
## Date:    29 July 2019

# Function to run specified number of simulations
run_debug_sims <- function(n_sims, run_debug_rerandomize, run_prev_ai, run_prop_ai, run_prop_adh, run_rr_ring) {
  output_list <- lapply(X   = 1:n_sims,
                        FUN = function(x) { 
                          run_sim_secondary(rr_ai        = params$rr_ai,
                                            prev_ai      = run_prev_ai,
                                            prop_ai      = run_prop_ai,
                                            prop_adh     = run_prop_adh,
                                            rr_ring      = run_rr_ring,
                                            i            = x,
                                            output_dir   = "",
                                            re_randomize = run_debug_rerandomize) })
  return(output_list)
}

# Function to plot histogram of per-act RRs. Includes all infections (possible to have > 1 per woman) and acts. 
plot_hist_per_act_rr <- function(output_list, n_sims, common_x_axis = F) {
  
  if(common_x_axis) {
    all_rr <- c(sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf)/mean(n_acts)]/x$ri_dt[arm == 0, mean(n_inf)/mean(n_acts)]),
                sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]/x$ri_dt[arm == 0, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]),
                sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_any_ai)/mean(n_acts_any_ai)]/x$ri_dt[arm == 0, mean(n_inf_any_ai)/mean(n_acts_any_ai)]))
    
    max_x_val <- max(all_rr)
    
    ret_plot <- ggplot() +
      geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf)/mean(n_acts)]/x$ri_dt[arm == 0, mean(n_inf)/mean(n_acts)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
      coord_cartesian(xlim = c(0, max_x_val)) +
      labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations"), subtitle = "Includes all infections (>1 per woman possible) and acts") +
      theme(text = element_text(family = "Times"),
            panel.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major = element_line(color = "grey", size = 0.2),
            panel.grid.minor = element_line(color = "grey", size = 0.1),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  } else {
    ret_plot <- ggplot() +
      geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf)/mean(n_acts)]/x$ri_dt[arm == 0, mean(n_inf)/mean(n_acts)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
      coord_cartesian(xlim = c(0, 1)) +
      labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations"), subtitle = "Includes all infections (>1 per woman possible) and acts") +
      theme(text = element_text(family = "Times"),
            panel.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            panel.grid.major = element_line(color = "grey", size = 0.2),
            panel.grid.minor = element_line(color = "grey", size = 0.1),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5))
  }
  
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

# Plot histogram of per-act relative risk of infection among women engaged exclusively in vaginal intercourse
plot_hist_per_act_rr_vi_excl <- function(output_list, n_sims) {
  
  all_rr <- c(sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf)/mean(n_acts)]/x$ri_dt[arm == 0, mean(n_inf)/mean(n_acts)]),
              sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]/x$ri_dt[arm == 0, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]),
              sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_any_ai)/mean(n_acts_any_ai)]/x$ri_dt[arm == 0, mean(n_inf_any_ai)/mean(n_acts_any_ai)]))
  
  max_x_val <- max(all_rr)
  
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]/x$ri_dt[arm == 0, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, max_x_val)) +
    labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations\namong women engaged in exclusively RVI"), subtitle = "Includes all infections (>1 per woman possible) and acts") +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Plot histogram of per-act relative risk of infection among women engaged exclusively in vaginal intercourse and with high adherence.
plot_hist_per_act_rr_vi_excl_adh <- function(output_list, n_sims, total = F) {
  
  all_rr <- c(sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf)/mean(n_acts)]/x$ri_dt[arm == 0, mean(n_inf)/mean(n_acts)]),
              sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]/x$ri_dt[arm == 0, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]),
              sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_any_ai)/mean(n_acts_any_ai)]/x$ri_dt[arm == 0, mean(n_inf_any_ai)/mean(n_acts_any_ai)]))
  
  max_x_val <- max(all_rr)
  
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_vi_excl_adh)/mean(n_acts_vi_excl_adh)]/x$ri_dt[arm == 0, mean(n_inf_vi_excl_adh)/mean(n_acts_vi_excl_adh)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, max_x_val)) +
    labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations\namong women engaged in exclusively RVI and with high adherence"), subtitle = "Includes all infections (>1 per woman possible) and acts") +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Plot histogram of per-act relative risk of infection among women engaged in both vaginal and anal intercourse
plot_hist_per_act_rr_any_ai <- function(output_list, n_sims) {
  
  all_rr <- c(sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf)/mean(n_acts)]/x$ri_dt[arm == 0, mean(n_inf)/mean(n_acts)]),
              sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]/x$ri_dt[arm == 0, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]),
              sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_any_ai)/mean(n_acts_any_ai)]/x$ri_dt[arm == 0, mean(n_inf_any_ai)/mean(n_acts_any_ai)]))
  
  max_x_val <- max(all_rr)
  
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_any_ai)/mean(n_acts_any_ai)]/x$ri_dt[arm == 0, mean(n_inf_any_ai)/mean(n_acts_any_ai)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, max_x_val)) +
    labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations\namong women engaged in any RAI"), subtitle = "Includes all infections (>1 per woman possible) and acts") +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Plot histogram of per-act relative risk of infection among women engaged in both vaginal and anal intercourse and with high adherence
plot_hist_per_act_rr_any_ai_adh <- function(output_list, n_sims) {
  
  all_rr <- c(sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf)/mean(n_acts)]/x$ri_dt[arm == 0, mean(n_inf)/mean(n_acts)]),
              sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]/x$ri_dt[arm == 0, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]),
              sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_any_ai)/mean(n_acts_any_ai)]/x$ri_dt[arm == 0, mean(n_inf_any_ai)/mean(n_acts_any_ai)]))
  
  max_x_val <- max(all_rr)
  
  ret_plot <- ggplot() +
    geom_histogram(data = data.table(rr = sapply(output_list, function(x) x$ri_dt[arm == 1, mean(n_inf_any_ai_adh)/mean(n_acts_any_ai_adh)]/x$ri_dt[arm == 0, mean(n_inf_any_ai_adh)/mean(n_acts_any_ai_adh)])), aes(x = rr), bins = 30, fill = "grey", color = "darkgrey") +
    coord_cartesian(xlim = c(0, max_x_val)) +
    labs(x = "Relative risk", y = "Count", title = paste0("Distribution of per-act relative risk across ", n_sims, " simulations\namong women engaged in any RAI and with high adherence"), subtitle = "Includes all infections (>1 per woman possible) and acts") +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(ret_plot)
}

# Function to plot the simulated and expected prevalence of RAI
plot_prev_ai <- function(output_list, n_sims, prev_ai) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$study_dt %>% group_by(site) %>% summarise(prev_ai = mean(prob_ai > 0)))), idcol = "sim")
  
  load(file = "~/Documents/code/aspire/data-public/ai_dt.RData")
  obs_dt <- prop_ai_dt[, .SD, .SDcols = c("site", prev_ai)]
  setnames(obs_dt, old = prev_ai, new = "prev_ai")
  
  # Factor site to print by country
  country_dt <- unique(output_list[[1]]$study_dt[, .(country, site)])
  output_dt <- merge(x = output_dt, y = country_dt, by = "site", all = T)
  output_dt[, site := factor(site, levels = rev(c("Blantyre", "Lilongwe", "RK Khan", "Botha's Hill", "Umkomaas", "Isipingo", "Tongaat", "Verulam", "eThekwini", "WHRI", "Emavundleni Cent", "MU-JHU", "Spilhaus", "Seke South", "Zengeza")))]
  
  obs_dt <- merge(x = obs_dt, y = country_dt, by = "site", all = T)
  
  ret_plot <- ggplot(data = output_dt) +
    geom_point(aes(x = prev_ai, y = site, col = country), alpha = 0.1) +
    geom_point(data = obs_dt, aes(x = eval(prev_ai), y = site)) +
    labs(x = "Prevalence of RAI", y = "", col = "Country", title = paste0("RAI prevalence by site across ", n_sims, " simulations"), subtitle = "Simulated data are in colored dots and expected values are in black") +
    scale_color_discrete(labels = c("Malawi", "South Africa", "Uganda", "Zimbabwe")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.key = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
  return(ret_plot)
}

# Function to plot simulated and expected prevalence of adherence - OLD
plot_prev_adh_old <- function(output_list, n_sims, prop_adh) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$adh_dt %>% group_by(site, f_age_cat %in% c("18-21", "22-26"), post_adh_intervention) %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  
  output_dt[, age_cat := ifelse(`f_age_cat %in% c(\"18-21\", \"22-26\")`, "18-26", "27-45")]
  
  load(file = "~/Documents/code/aspire/data-public/adh_dt.RData")
  obs_dt <- prop_adh_dt[, .SD, .SDcols = c("site", "age_cat", "post_adh_intervention", prop_adh)]
  setnames(obs_dt, old = prop_adh, new = "prop_adh")
  
  # Factor site to print by country
  country_dt <- unique(output_list[[1]]$study_dt[, .(country, site)])
  output_dt <- merge(x = output_dt, y = country_dt, by = "site", all = T)
  output_dt[, site := factor(site, levels = c("Blantyre", "Lilongwe", "RK Khan", "Botha's Hill", "Umkomaas", "Isipingo", "Tongaat", "Verulam", "eThekwini", "WHRI", "Emavundleni Cent", "MU-JHU", "Spilhaus", "Seke South", "Zengeza"))]
  
  obs_dt <- merge(x = obs_dt, y = country_dt, by = "site", all = T)
  
  ret_plot <- ggplot(data = output_dt) +
    geom_point(aes(x = site, y = prop_adh, col = country, shape = age_cat), position = position_dodge(width = 0.5), alpha = 0.1) +
    geom_point(data = obs_dt, aes(x = site, y = eval(prop_adh), shape = age_cat), position = position_dodge(width = 0.5)) +
    labs(x = "", y = "Proportion adherent", col = "Country", shape = "Age group", title = paste0("Proportion adherent by site across ", n_sims, " simulations"), subtitle = "Simulated data are in colored dots and expected values are in black") +
    scale_color_discrete(labels = c("Malawi", "South Africa", "Uganda", "Zimbabwe")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.key = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank()) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    facet_grid(. ~ post_adh_intervention, labeller = as_labeller(c(`0` = "Pre-adherence interventions", `1` = "Post-adherence interventions"))) +
    coord_cartesian(ylim = c(0, 1))
  
  return(ret_plot)
}

# Function to plot simulated and expected prevalence of adherence
plot_prev_adh <- function(output_list, n_sims) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% group_by(visit) %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  
  load("~/Documents/code/aspire/data-private/adh_dt.RData")
  obs_dt <- as.data.table(dt %>% group_by(visit) %>% summarise(prop_adh = mean(adh)))
  obs_dt[, visit := as.integer(gsub(pattern = "01", replacement = "", x = visit))]
  
  output_dt[, scenario := "Simulated"]
  obs_dt[, scenario := "Observed"]
  
  ggplot() +
    geom_line(data = output_dt, aes(x = visit, y = prop_adh, group = sim, color = scenario), alpha = 0.15) +
    geom_smooth(data = output_dt[visit <= 33], aes(x = visit, y = prop_adh, color = scenario), method = "loess", formula = y ~ x, se = F) +
    geom_line(data = obs_dt, aes(x = visit, y = prop_adh, color = scenario)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Visit", y = "Proportion adherent", color = "", title = paste0("Proportion adherent (marginal) across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "grey", size = 0.15),
          legend.key = element_blank())
}

# Function to plot the proportion of overall visits at which women were adherent
plot_proportion_visits_adh <- function(output_list, n_sims, target_adh) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  
  ggplot(data = output_dt) +
    geom_boxplot(aes(x = 1, y = prop_adh)) +
    geom_hline(aes(yintercept = target_adh), color = "grey", linetype = "dashed") +
    labs(x = "", y = "Proportion visits at\nwhich women were adherent") +
    coord_cartesian(ylim = c(target_adh - 0.05, target_adh + 0.05)) +
    theme(text = element_text(family = "Times"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

# Function to plot the proportion of overall visits at which women were adherent, using all hypothetical adherence values
plot_proportion_visits_adh_all <- function(output_dir, output_lists, n_sims, target_adhs) {
  
  output_dt <- data.table(expand.grid(sim = 1:n_sims, adh = target_adhs))
  
  for(i in 1:length(output_lists)) {
    load(paste0(output_dir, "/results_list_", output_lists[i], ".RData"))
    
    output_dt[adh == target_adhs[i], prop_adh := rbindlist(l = lapply(results_list, function(x) as.data.table(x$long_dt %>% summarise(prop_adh = mean(adh)))))[, prop_adh]]
  }
  
  ggplot(data = output_dt, aes()) +
    geom_boxplot(aes(x = adh, y = prop_adh, group = adh, fill = as.factor(adh))) +
    labs(x = "Target proportion adherent", y = "Simulated proportion adherent", title = paste0("Target and simulated adherence across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.grid.major = element_line(color = "grey", size = 0.1),
          panel.grid.minor.y = element_line(color = "grey", size = 0.1),
          panel.border = element_rect(color = "black", fill = NA),
          legend.position = "none")
}

# Functions to plot observed and simulated adherence by characteristics associated with adherence
plot_adh_by_characteristic <- function(output_list, n_sims) {
  
  # Ring worries
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% group_by(visit, ring_worries) %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  
  load("~/Documents/code/aspire/data-private/adh_dt.RData")
  obs_dt <- as.data.table(dt %>% group_by(visit, ring_worries) %>% summarise(prop_adh = mean(adh)))
  obs_dt[, visit := as.integer(gsub(pattern = "01", replacement = "", x = visit))]
  
  output_dt[, scenario := "Simulated"]
  obs_dt[, scenario := "Observed"]
  
  plot_ring_worries <- ggplot() +
    geom_line(data = output_dt, aes(x = visit, y = prop_adh, group = interaction(sim, ring_worries), color = ring_worries, linetype = scenario), alpha = 0.15) +
    geom_smooth(data = output_dt[visit <= 33], aes(x = visit, y = prop_adh, color = ring_worries, linetype = scenario), method = "loess", formula = y ~ x, se = F) +
    geom_line(data = obs_dt, aes(x = visit, y = prop_adh, color = ring_worries, linetype = scenario)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Visit", y = "Proportion adherent", color = "Worried about using ring", linetype = "", title = paste0("Proportion adherent by ring worries across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "grey", size = 0.15),
          legend.key = element_blank())
  
  # Age
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% group_by(visit, f_age_cat == "27-45") %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  setnames(output_dt, old = "f_age_cat == \"27-45\"", new = "age_cat")
  output_dt[, age_cat := ifelse(age_cat == "TRUE", "27-45", "18-26")]
  
  obs_dt <- as.data.table(dt %>% group_by(visit, age_cat) %>% summarise(prop_adh = mean(adh)))
  obs_dt[, visit := as.integer(gsub(pattern = "01", replacement = "", x = visit))]
  
  output_dt[, scenario := "Simulated"]
  obs_dt[, scenario := "Observed"]
  
  plot_age <- ggplot() +
    geom_line(data = output_dt, aes(x = visit, y = prop_adh, group = interaction(sim, age_cat), color = age_cat, linetype = scenario), alpha = 0.15) +
    geom_smooth(data = output_dt[visit <= 33], aes(x = visit, y = prop_adh, color = age_cat, linetype = scenario), method = "loess", formula = y ~ x, se = F) +
    geom_line(data = obs_dt, aes(x = visit, y = prop_adh, color = age_cat, linetype = scenario)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Visit", y = "Proportion adherent", color = "Age", linetype = "", title = paste0("Proportion adherent by age across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "grey", size = 0.15),
          legend.key = element_blank())
  
  # Participant aware in study
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% group_by(visit, pknow_study) %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  
  obs_dt <- as.data.table(dt %>% group_by(visit, pknow_study) %>% summarise(prop_adh = mean(adh)))
  obs_dt[, visit := as.integer(gsub(pattern = "01", replacement = "", x = visit))]
  
  output_dt[, scenario := "Simulated"]
  obs_dt[, scenario := "Observed"]
  
  plot_pknow <- ggplot() +
    geom_line(data = output_dt, aes(x = visit, y = prop_adh, group = interaction(sim, pknow_study), color = pknow_study, linetype = scenario), alpha = 0.15) +
    geom_smooth(data = output_dt[visit <= 33], aes(x = visit, y = prop_adh, color = pknow_study, linetype = scenario), method = "loess", formula = y ~ x, se = F) +
    geom_line(data = obs_dt, aes(x = visit, y = prop_adh, color = pknow_study, linetype = scenario)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Visit", y = "Proportion adherent", color = "Partner aware of study participation", linetype = "", title = paste0("Proportion adherent by partner awareness of\nstudy participation across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "grey", size = 0.15),
          legend.key = element_blank())
  
  # Same primary sex partner in past three months
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% group_by(visit, same_psp) %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  
  obs_dt <- as.data.table(dt %>% group_by(visit, same_psp) %>% summarise(prop_adh = mean(adh)))
  obs_dt[, visit := as.integer(gsub(pattern = "01", replacement = "", x = visit))]
  
  output_dt[, scenario := "Simulated"]
  obs_dt[, scenario := "Observed"]
  
  plot_same_psp <- ggplot() +
    geom_line(data = output_dt, aes(x = visit, y = prop_adh, group = interaction(sim, same_psp), color = same_psp, linetype = scenario), alpha = 0.15) +
    geom_smooth(data = output_dt[visit <= 33], aes(x = visit, y = prop_adh, color = same_psp, linetype = scenario), method = "loess", formula = y ~ x, se = F) +
    geom_line(data = obs_dt, aes(x = visit, y = prop_adh, color = same_psp, linetype = scenario)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Visit", y = "Proportion adherent", color = "Same primary partner\npast three months", linetype = "", title = paste0("Proportion adherent by same primary partner\nin prior 3 months across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "grey", size = 0.15),
          legend.key = element_blank())
  
  # Use of OCP at enrollment
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% group_by(visit, bl_ocp) %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  output_dt[, bl_ocp := ifelse(bl_ocp == 1, "OCP", "Other")]
  
  obs_dt <- as.data.table(dt %>% group_by(visit, fp_method) %>% summarise(prop_adh = mean(adh)))
  setnames(obs_dt, old = "fp_method", new = "bl_ocp")
  obs_dt[, visit := as.integer(gsub(pattern = "01", replacement = "", x = visit))]
  
  output_dt[, scenario := "Simulated"]
  obs_dt[, scenario := "Observed"]
  
  plot_bl_ocp <- ggplot() +
    geom_line(data = output_dt, aes(x = visit, y = prop_adh, group = interaction(sim, bl_ocp), color = bl_ocp, linetype = scenario), alpha = 0.15) +
    geom_smooth(data = output_dt[visit <= 33], aes(x = visit, y = prop_adh, color = bl_ocp, linetype = scenario), method = "loess", formula = y ~ x, se = F) +
    geom_line(data = obs_dt, aes(x = visit, y = prop_adh, color = bl_ocp, linetype = scenario)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Visit", y = "Proportion adherent", color = "Baseline contraceptive method", linetype = "", title = paste0("Proportion adherent by baseline contraceptive method\nacross ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "grey", size = 0.15),
          legend.key = element_blank())
  
  # Site
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% group_by(visit, site) %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  
  obs_dt <- as.data.table(dt %>% group_by(visit, site) %>% summarise(prop_adh = mean(adh)))
  obs_dt[, visit := as.integer(gsub(pattern = "01", replacement = "", x = visit))]
  
  output_dt[, scenario := "Simulated"]
  obs_dt[, scenario := "Observed"]
  
  plot_site <- ggplot() +
    geom_line(data = output_dt, aes(x = visit, y = prop_adh, group = interaction(sim, site), color = site, linetype = scenario), alpha = 0.05) +
    geom_smooth(data = output_dt[visit <= 33], aes(x = visit, y = prop_adh, color = site, linetype = scenario), method = "loess", formula = y ~ x, se = F) +
    geom_line(data = obs_dt, aes(x = visit, y = prop_adh, color = site, linetype = scenario)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Visit", y = "Proportion adherent", color = "Site", linetype = "", title = paste0("Proportion adherent by site across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "grey", size = 0.15),
          legend.key = element_blank(),
          strip.background = element_blank(),
          legend.text = element_text(size = 10)) +
    facet_grid(. ~ scenario)
  
  # Calendar time > 12 months since study activation
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% group_by(visit, post_adh_intervention) %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  output_dt[, post_adh_intervention := ifelse(post_adh_intervention == 1, "Yes", "No")]
  
  obs_dt <- as.data.table(dt %>% group_by(visit, post_adh_intervention) %>% summarise(prop_adh = mean(adh)))
  obs_dt[, post_adh_intervention := ifelse(post_adh_intervention == 1, "Yes", "No")]
  obs_dt[, visit := as.integer(gsub(pattern = "01", replacement = "", x = visit))]
  
  output_dt[, scenario := "Simulated"]
  obs_dt[, scenario := "Observed"]
  
  plot_cal_time <- ggplot() +
    geom_line(data = output_dt, aes(x = visit, y = prop_adh, group = interaction(sim, post_adh_intervention), color = post_adh_intervention, linetype = scenario), alpha = 0.15) +
    geom_smooth(data = output_dt[visit <= 33], aes(x = visit, y = prop_adh, color = post_adh_intervention, linetype = scenario), method = "loess", formula = y ~ x, se = F) +
    geom_line(data = obs_dt, aes(x = visit, y = prop_adh, color = post_adh_intervention, linetype = scenario)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Visit", y = "Proportion adherent", color = "Calendar time > 12 months\nsince study activation", linetype = "", title = paste0("Proportion adherent by calendar time across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "grey", size = 0.15),
          legend.key = element_blank())
  
  # Time on study > 3 months
  output_dt <- rbindlist(l = lapply(output_list, function(x) as.data.table(x$long_dt %>% group_by(visit, time_study_gt_three_mos) %>% summarise(prop_adh = mean(adh)))), idcol = "sim")
  output_dt[, time_study_gt_three_mos := ifelse(time_study_gt_three_mos == 1, "Yes", "No")]
  
  obs_dt <- as.data.table(dt %>% group_by(visit, time_study_gt_three_mos) %>% summarise(prop_adh = mean(adh)))
  obs_dt[, time_study_gt_three_mos := ifelse(time_study_gt_three_mos == 1, "Yes", "No")]
  obs_dt[, visit := as.integer(gsub(pattern = "01", replacement = "", x = visit))]
  
  output_dt[, scenario := "Simulated"]
  obs_dt[, scenario := "Observed"]
  
  plot_time_study <- ggplot() +
    geom_line(data = output_dt, aes(x = visit, y = prop_adh, group = interaction(sim, time_study_gt_three_mos), color = time_study_gt_three_mos), alpha = 0.15) +
  geom_smooth(data = output_dt[visit <= 33], aes(x = visit, y = prop_adh, color = time_study_gt_three_mos), method = "loess", formula = y ~ x, se = F) +
    geom_point(data = obs_dt, aes(x = visit, y = prop_adh, color = time_study_gt_three_mos)) +
    coord_cartesian(ylim = c(0.5, 1)) +
    labs(x = "Visit", y = "Proportion adherent", color = "On study > 3 months\n(obs. data in points)", title = paste0("Proportion adherent by time on study across ", n_sims, " simulations")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid = element_line(color = "grey", size = 0.15),
          legend.key = element_blank())
  
  plot_list <- list(plot_age, plot_site, plot_pknow, plot_same_psp, plot_bl_ocp, plot_ring_worries, plot_cal_time, plot_time_study)
  names(plot_list) <- c("plot_age", "plot_site", "plot_pknow", "plot_same_psp", "plot_bl_ocp", "plot_ring_worries", "plot_cal_time", "plot_time_study")
  
  return(plot_list)
}

# Function to plot exposure rate ratio over time
plot_exp_rr <- function(output_list, n_sims) {
  output_dt_exp  <- rbindlist(l = lapply(output_list, function(x) x$exposure_dt), idcol = "sim")
  output_dt_risk <- rbindlist(l = lapply(output_list, function(x) x$cumulative_risk_dt), idcol = "sim")
  
  output_dt <- merge(x = output_dt_exp, y = output_dt_risk, by = c("sim", "visit", "arm"), all = T)
  
  min_time_val <- min(output_dt[, max(visit), by = sim]$V1)
  
  # Convert dataset from long to wide
  output_dt_wide <- dcast(data = output_dt, formula = sim + visit ~ arm, value.var = c("exposures", "mean_cum_risk"))
  
  # Create new variable for exposure rate ratio
  output_dt_wide[, `:=`(exp_rr      = exposures_1/exposures_0,
                        cum_risk_rr = mean_cum_risk_1/mean_cum_risk_0)]
  
  # Calculate mean exposure rate ratio at each time step across simulations. Include the mean only for those simulations with the maximum follow-up time observed across all sims.
  mean_dt <- output_dt_wide[visit <= min_time_val, .(exp_rr = mean(exp_rr, na.rm = T), cum_risk_rr = mean(cum_risk_rr, na.rm = T)), by = visit]
  
  # min_y_val <- min(c(output_dt_wide[, min(exp_rr, cum_risk_rr, na.rm = T)], 1))
  # max_y_val <- min(7, output_dt_wide[!is.infinite(exp_rr), max(exp_rr, cum_risk_rr, na.rm = T)])
  
  panel_1 <- ggplot(data = output_dt_wide[!is.na(exp_rr)]) +
    geom_line(aes(x = visit, y = exp_rr, group = sim), col = "blue", alpha = 0.05) +
    geom_point(data = mean_dt, aes(x = visit, y = exp_rr), col = "blue") + 
    labs(x = "Time (months)", y = "Exposure rate ratio") +
    coord_cartesian(ylim = c(0.5, 3)) +
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
    coord_cartesian(ylim = c(0.5, 3)) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          strip.background = element_blank())
  
  plot_together <- ggarrange(panel_1, panel_2, ncol = 2, nrow = 1)
  
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
    geom_line(aes(x = month, y = V1, group = id), alpha = 0.5, color = "blue") +
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

# Function to plot follow-up time by country
plot_fup_time_country <- function(output_list, n_sims) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$study_dt[, .(fup_years = sum(visit * 30)/365.25), by = country]), idcol = "sim")
  
  # Create data table for observed follow-up time (calculated from variable fu_yrs hiv.RData)
  dt_obs <- data.table(country = c("mal", "sa", "uga", "zim"), fup_years = c(330.8329, 2425.8493, 456.6849, 1066.8274))
  
  ret_plot <- ggplot(data = output_dt) +
    geom_point(aes(x = country, y = fup_years, col = country), alpha = 0.05) +
    geom_point(data = dt_obs, aes(x = country, y = fup_years)) +
    labs(x = "", y = "Years of follow-up", title = paste0("Follow-up time by country across ", n_sims, " simulations"), subtitle = "Simulated data in colors, observed data in black") +
    scale_x_discrete(labels = c("mal" = "Malawi", "sa" = "South Africa", "uga" = "Uganda", "zim" = "Zimbabwe")) +
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.grid.major = element_line(color = "grey", size = 0.2),
          panel.grid.minor = element_line(color = "grey", size = 0.1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none") +
    guides(colour = guide_legend(override.aes = list(alpha = 1)))
  
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
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), RR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$ri_dt[arm == 1, mean(n_inf)/mean(n_acts)]/x$ri_dt[arm == 0, mean(n_inf)/mean(n_acts)], 2)))))))
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

create_table_per_act_rr_vi_excl <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), RR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$ri_dt[arm == 1, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)]/x$ri_dt[arm == 0, mean(n_inf_vi_excl)/mean(n_acts_vi_excl)], 2)))))))
}

create_table_per_act_rr_vi_excl_adh <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), RR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$ri_dt[arm == 1, mean(n_inf_vi_excl_adh)/mean(n_acts_vi_excl_adh)]/x$ri_dt[arm == 0, mean(n_inf_vi_excl_adh)/mean(n_acts_vi_excl_adh)], 2)))))))
}

create_table_per_act_rr_any_ai <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), RR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$ri_dt[arm == 1, mean(n_inf_any_ai)/mean(n_acts_any_ai)]/x$ri_dt[arm == 0, mean(n_inf_any_ai)/mean(n_acts_any_ai)], 2)))))))
}

create_table_per_act_rr_any_ai_adh <- function(output_list) {
  return(data.table(`Summary measure` = c("Minimum", "1st Quartile", "Median", "Mean", "3rd Quartile", "Maximum"), RR = as.numeric(unname(summary(sapply(output_list, function(x) round(x$ri_dt[arm == 1, mean(n_inf_any_ai_adh)/mean(n_acts_any_ai_adh)]/x$ri_dt[arm == 0, mean(n_inf_any_ai_adh)/mean(n_acts_any_ai_adh)], 2)))))))
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

# Function to create a barplot of how each factor contributes to efficacy dilution. Specify user-input per-act RR; dataset that contains model-output per-act RR (including imbalance in arms) and HR with heterogeneity in risk; dataset that contains proportion efficacy dilution attributable to RAI and non-adherence and simulated HR.
plot_bar_efficacy_dilution <- function(sims_directory_best_eff, sims_directory_rai_nonadh, results_best_possible_eff, results_rai_nonadh, rr_ring, n_sims) {
  load(paste0("~/Documents/code/aspire/dev/test-sims/debug/", sims_directory_best_eff, "/", results_best_possible_eff, ".RData"))
  output_best_possible_eff <- copy(results_list)
  
  load(paste0("~/Documents/code/aspire/dev/test-sims/debug/", sims_directory_rai_nonadh, "/", results_rai_nonadh, ".RData"))
  output_rai_nonadh <- copy(results_list)
  
  # Estimate median per-act RR from best-possible simulations (i.e., perfect adherence, no RAI)
  med_per_act_rr <- median(sapply(output_best_possible_eff, function(x) x$ri_dt[arm == 1, mean(n_inf)/mean(n_acts)]/x$ri_dt[arm == 0, mean(n_inf)/mean(n_acts)]))
  
  # Estimate median HR from best-possible simulations (i.e., perfect adherence, no RAI)
  med_hr_best_possible <- median(sapply(output_best_possible_eff, function(x) x$hr))
  
  # Estimate median HR from simulations with non-adherence and RAI
  med_hr_trial <- median(sapply(output_rai_nonadh, function(x) x$hr))
  
  # Estimate the median proportion of efficacy dilution attributable to RAI and non-adherence
  med_dilution_nonadh <- median(sapply(output_rai_nonadh, function(x) x$prop_eff_dilution_non_adh))
  
  # Create table to store the proportion of efficacy dilution (from a perfect product) attributable to each factor
  dt <- data.table("Source of efficacy dilution" = c("Ring failure to provide protection", "Imbalance in arms", "Heterogeneity in risk", "Non-adherence", "RAI", ""), value = c(rr_ring, med_per_act_rr - rr_ring, med_hr_best_possible - med_per_act_rr, (med_hr_trial - med_hr_best_possible) * med_dilution_nonadh, (med_hr_trial - med_hr_best_possible) * (1 - med_dilution_nonadh), 1 - med_hr_trial))
  
  # Factor source of efficacy dilution to plot in desired order
  dt[, `Source of efficacy dilution` := factor(`Source of efficacy dilution`, levels = rev(c("Ring failure to provide protection", "Imbalance in arms", "Heterogeneity in risk", "Non-adherence", "RAI", "")))]
  
  col_func <- colorRampPalette(colors = c("#B2D35C", "#DBEABC", "#B0A2C5", "#844594"))
  
  cols <- c("white", col_func(nrow(dt) - 1))
  
  panel_1 <- ggplot(data = dt, aes(x = 1, y = value, fill = `Source of efficacy dilution`)) + 
    geom_bar(stat = "identity") + 
    geom_hline(yintercept = 0.73, linetype = "dashed") +
    scale_fill_manual(values = cols) + 
    labs(x = "", y = "Hazard ratio", fill = "", title = paste0("Median"))  + 
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(), 
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12)) +
    coord_flip()
  
  # Include plot of interquartile range of simulations with non-adherence and RAI
  sims_dt <- data.table(hr_trial = sapply(output_rai_nonadh, function(x) x$hr), dilution_nonadh = sapply(output_rai_nonadh, function(x) x$prop_eff_dilution_non_adh))
  
  sims_dt <- sims_dt[hr_trial >= quantile(x = hr_trial, probs = 0.25) & hr_trial <= quantile(x = hr_trial, probs = 0.75)]
  
  # Create table to store the proportion of efficacy dilution (from a perfect product) attributable to each factor
  dt <- data.table("Ring failure to provide protection" = rep(rr_ring, nrow(sims_dt)), "Imbalance in arms" = rep(med_per_act_rr - rr_ring, nrow(sims_dt)), "Heterogeneity in risk" = rep(med_hr_best_possible - med_per_act_rr, nrow(sims_dt)), "Non-adherence" = sims_dt[, (hr_trial - med_hr_best_possible) * dilution_nonadh], "RAI" = sims_dt[, (hr_trial - med_hr_best_possible) * (1 - dilution_nonadh)], " " = sims_dt[, 1 - hr_trial], sim = 1:nrow(sims_dt))
  
  melt_dt <- melt(dt, id.vars = "sim", measure.vars = names(dt)[!(names(dt) %in% "sim")])
  
  # Factor source of efficacy dilution to plot in desired order
  melt_dt[, variable := factor(variable, levels = rev(c("Ring failure to provide protection", "Imbalance in arms", "Heterogeneity in risk", "Non-adherence", "RAI", " ")))]
  
  panel_2 <- ggplot(data = melt_dt, aes(x = sim, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.73, linetype = "dashed") +
    scale_fill_manual(values = cols) + 
    labs(x = "", y = "Hazard ratio", fill = "", title = paste0("Individual simulations with non-adherence and RAI\n(incl. interquartile range of the HR)")) + 
    theme(text = element_text(family = "Times"),
          panel.background = element_blank(), 
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12)) +
    coord_flip()
  
  plot_together <- ggarrange(panel_1, panel_2, nrow = 2)
  
  return(annotate_figure(p = plot_together, top = text_grob(paste0("Attribution of each source of efficacy dilution across ", n_sims, " simulations"), family = "Times", size = 14, face = "bold")))
}

# Function to plot incidence by age
plot_inc_by_age <- function(output_list) {
  # Calculate incidence by age group in the placebo arm
  inc_dt_list <- lapply(output_list, function(x) as.data.table(x$study_dt[arm == 0] %>% group_by(f_age_cat) %>% summarise(inc_rate = (sum(hiv)/sum(as.numeric(visit * 30)/365)) * 100)))
  
  lapply(inc_dt_list, function(x) setnames(x, old = "f_age_cat", new = "age_cat"))
  
  # Combine all items in inc_list to one data.table
  inc_sim_dt <- rbindlist(l = inc_dt_list)
  
  # Specify whether data is from a simulation or observed
  inc_sim_dt[, data_source := "sim"]
  
  load("~/Documents/code/aspire/data-public/inf_obs.RData")
  
  inc_obs <- inf_obs[1:3, .(age_cat, hiv_inf, n, py, inc_rate)]
  inc_obs[, data_source := "obs"]
  
  # Combine all incidence data
  inc_dt <- rbindlist(l = list(inc_sim_dt, inc_obs[, .(age_cat, inc_rate, data_source)]))
  
  # Plot simulated and observed age-specific incidence
  ret_plot <- ggplot(data = inc_dt) +
    geom_jitter(data = inc_dt[data_source == "sim"], aes(x = age_cat, y = inc_rate, col = data_source), width = 0.1, alpha = 0.5) +
    geom_point(data = inc_dt[data_source == "obs"], aes(x = age_cat, y = inc_rate, col = data_source), size = 2.25) +
    geom_errorbar(data = inf_obs[1:3], aes(x = age_cat, ymin = inc_rate - 1.96 * se, ymax = inc_rate + 1.96 * se), width = 0.1, size = 0.75) +
    scale_color_manual(labels = c("Observed incidence", "Simulated incidence"), values = c("black", "#AF3AF5")) +
    coord_cartesian(ylim = c(0, 9)) +
    labs(x = "Age group", y = "Incidence rate (per 100 woman-years)", col = "") +
    theme(text = element_text(family = "Times"),
          panel.border = element_rect(color = "black", fill = NA),
          panel.background = element_blank(),
          panel.grid = element_line(size = 0.05, colour = "black"),
          legend.key = element_blank())
  
  return(ret_plot)
}

# Function to plot incidence in months 15-23 and 24-31
plot_inc_time_period <- function(output_list) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$long_dt[visit %in% 15:31 & arm == 0]), idcol = "sim")
  
  # Create time period variable
  output_dt[, time_period := ifelse(visit %in% 15:23, "15-23 months", "24-31 months")]
  
  # Keep last row by id, sim, and time period
  output_dt <- output_dt[, .SD[.N], by = c("id", "sim", "time_period")]
  
  inc_sim_dt_1523 <- as.data.table(output_dt[time_period == "15-23 months"] %>% group_by(sim) %>% summarise(inc_rate = sum(hiv)/sum(((visit - 14) * 30)/365) * 100))
  inc_sim_dt_2431 <- as.data.table(output_dt[time_period == "24-31 months"] %>% group_by(sim) %>% summarise(inc_rate = sum(hiv)/sum(((visit - 23) * 30)/365) * 100))
  
  load("~/Documents/code/aspire/data-public/inf_obs.RData")
  
  inc_dt <- rbindlist(l = list(inf_obs[4:5, .(time_period = c("15-23 months", "24-31 months"), inc_rate, data_source = "obs")], inc_sim_dt_1523[, .(time_period = "15-23 months", inc_rate, data_source = "sim")], inc_sim_dt_2431[, .(time_period = "24-31 months", inc_rate, data_source = "sim")]))
  
  inf_obs[, time_period := c(NA, NA, NA, "15-23 months", "24-31 months")]
  
  # Plot simulated and observed age-specific incidence
  ret_plot <- ggplot(data = inc_dt) +
    geom_jitter(data = inc_dt[data_source == "sim"], aes(x = time_period, y = inc_rate, col = data_source), width = 0.1, alpha = 0.5) +
    geom_point(data = inc_dt[data_source == "obs"], aes(x = time_period, y = inc_rate, col = data_source), size = 2.25) +
    geom_errorbar(data = inf_obs[4:5], aes(x = time_period, ymin = inc_rate - 1.96 * se, ymax = inc_rate + 1.96 * se), width = 0.1, size = 0.75) +
    scale_color_manual(labels = c("Observed incidence", "Simulated incidence"), values = c("black", "#AF3AF5")) +
    coord_cartesian(ylim = c(0, 9)) +
    labs(x = "Time period", y = "Incidence rate (per 100 woman-years)", col = "") +
    theme(text = element_text(family = "Times"),
          panel.border = element_rect(color = "black", fill = NA),
          panel.background = element_blank(),
          panel.grid = element_line(size = 0.05, colour = "black"),
          legend.key = element_blank())
  
  return(ret_plot)
}

# Function to plot incidence over time
plot_inc_time <- function(output_list) {
  output_dt <- rbindlist(l = lapply(output_list, function(x) x$long_dt[arm == 0, .(id, visit, hiv)]), idcol = "sim")
  
  min_visit_val <- min(output_dt[, max(visit), by = "sim"]$V1)
  
  inc_dt <- data.table(visit = 1:min_visit_val)
  
  for(i in 1:min_visit_val) {
    temp_dt <- output_dt[visit %in% 1:i]
    
    # Keep only last observation for each participant
    temp_dt <- temp_dt[, .SD[.N], by = c("sim", "id")]
    
    inc_dt[visit == i, inc_rate := temp_dt[, sum(hiv)/(sum(visit * 30)/365) * 100]]
  }
  
  max_y_val <- inc_dt[, max(inc_rate)]
  
  # Load observed data
  load("~/Documents/code/aspire/data-public/inc_obs_time.RData")
  
  ggplot(data = inc_dt) +
    geom_line(aes(x = visit, y = inc_rate), col = "#AF3AF5") +
    geom_point(data = inc_obs_time_dt, aes(x = visit, y = est)) +
    geom_errorbar(data = inc_obs_time_dt, aes(x = visit, ymin = lower, ymax = upper), width = 0.4) +
    coord_cartesian(ylim = c(0, max_y_val)) +
    labs(x = "Visit", y = "Incidence (per 100 py)") +
    theme(text = element_text(family = "Times"),
          panel.border = element_rect(color = "black", fill = NA),
          panel.background = element_blank(),
          panel.grid = element_line(size = 0.05, colour = "black"),
          legend.key = element_blank())
}
