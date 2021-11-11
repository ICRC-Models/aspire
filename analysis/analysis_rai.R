## Purpose: Analyze data and create figures for RAI brief report
## Author:  Kathryn Peebles
## Date:    31 October 2019

# Attach packages
library(data.table)
library(ggplot2)
library(dplyr)
library(scales)
library(viridis)
library(colorspace)

# Specify input directory of data objects
input_dir <- "~/Documents/code/aspire/analysis/data"
fig_dir   <- "~/Documents/Education/UW Epi Program/ASPIRE/Manuscript - RAI math model/Figures"

# ASPIRE supplement colors
# Dark purple: "#5F2285"
# Medium purple: "#8C30C4"
# Light purple: "#AF3AF5"
# Dark green: "#5C6C39"
# Light green: "#D5FD83"
# Alt green: "#B2D25C"

# MTN colors
# Light green:  #E0EAC0
# Dark green:   #BAD47B
# Dark purple:  #9C6EA9
# Light purple: #BB9FC7

aspire_green  <- "#B2D35C"
aspire_light_green <- "#DBEABC"
aspire_purple <- "#B0A2C5"
aspire_dark_purple <- "#844594"

#### PREVALENCE OF RAI AND NON-ADHERENCE ####
# Load data from simulations in which the proportion adherent and prevalence of RAI were systematically varied
load(paste0(input_dir, "/sim_results_adh_ai_dt.RData"))

# Function to calculate geometric mean
geo_mean <- function(x) { return((prod(x, na.rm = T) ^ (1/length(x)))) }

# Calculate the geometric mean hazard ratio for each combination of simulated efficacy, prevalence of adherence, and prevalence of RAI
hr_dt <- as.data.table(dt %>% group_by(prop_adh, prev_ai, rr_ring) %>% summarise(geo_mean_hr = geo_mean(hr)))
hr_dt[, prop_non_adh := 1 - prop_adh]

hr_loess <- loess(geo_mean_hr ~ prev_ai * prop_non_adh * rr_ring, data = hr_dt, degree = 2, span = 2.5)
hr_fit   <- as.data.table(expand.grid(list(prev_ai      = seq(0, 30, 0.5),
                                           prop_non_adh = seq(0, 0.50, 0.005),
                                           rr_ring      = seq(0.25, 0.35, 0.05))))
hr_fit[, hr := as.numeric(predict(hr_loess, newdata = hr_fit))]

# Refactor efficacy values to plot in desired order
hr_fit[, rr_ring := factor(rr_ring, levels = c(0.25, 0.30, 0.35))]

col_func <- colorRampPalette(colors = c(lighten("#B6A0C4", amount = 0.5), "#9670A5"))
cols <- col_func(n = 15)
purple_cols <- cols[seq(1, 12, length.out = 5)]

col_func <- colorRampPalette(colors = c(lighten("#B6A0C4", amount = 0.5), "#DFEAC0"))
# cols <- col_func(n = dt[, length(unique(rr_ai))])

light_purple <- lighten("#B6A0C4", amount = 0.5)

hr_fit[, prev_ai := prev_ai/100]

# Figure in color
ggplot(data = hr_fit) +
  geom_raster(aes(x = prev_ai, y = prop_non_adh, fill = hr)) +
  scale_fill_gradientn(colors = c("#A683B4", "#A683B4", light_purple, "#A683B4", "#A683B4"), values = rescale(x = c(hr_fit[, min(hr)], 0.73, hr_fit[, max(hr)]))) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Prevalence of receptive anal intercourse", y = "Prevalence of non-adherence", fill = "Simulated HR") +
  facet_grid(. ~ rr_ring, labeller = as_labeller(c("0.25" = "75% efficacy\nper vaginal exposure", "0.3" = "70% efficacy\nper vaginal exposure", "0.35" = "65% efficacy\nper vaginal exposure"))) +
  theme(text = element_text(family = "Times"),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggsave(filename = paste0(fig_dir, "/figure_hr_adh_rai_eff_color.pdf"), width = 8, height = 4)

# Figure in color for dissertation
hr_fit[, rr_ring := factor(rr_ring, levels = c(0.35, 0.3, 0.25))]

ggplot(data = hr_fit) +
  geom_raster(aes(x = prev_ai, y = prop_non_adh, fill = hr)) +
  scale_fill_gradientn(colors = c("#A683B4", "#A683B4", light_purple, "#A683B4", "#A683B4"), values = rescale(x = c(hr_fit[, min(hr)], 0.73, hr_fit[, max(hr)]))) +
  labs(x = "Prevalence of anal intercourse", y = "Prevalence of non-adherence", fill = "Simulated HR") +
  facet_grid(. ~ rr_ring, labeller = as_labeller(c("0.25" = "Efficacy = 75%", "0.3" = "Efficacy = 70%", "0.35" = "Efficacy = 65%"))) +
  theme(text = element_text(family = "sans"),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/Figure_HR_Adh_AI_Efficacy.pdf", width = 8, height = 4)

# Figure in greyscale
ggplot(data = hr_fit) +
  geom_raster(aes(x = prev_ai, y = prop_non_adh, fill = hr)) +
  scale_fill_gradientn(colors = c("gray27", "gray45", "gray94", "gray45", "gray27"), values = rescale(x = c(hr_fit[, min(hr)], 0.73, hr_fit[, max(hr)]))) +
  labs(x = "Prevalence of receptive anal intercourse", y = "Prevalence of non-adherence", fill = "Simulated HR") +
  facet_grid(. ~ rr_ring, labeller = as_labeller(c("0.25" = "Efficacy = 75%", "0.3" = "Efficacy = 70%", "0.35" = "Efficacy = 65%"))) +
  theme(text = element_text(family = "Times"),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggsave(filename = paste0(fig_dir, "/figure_hr_adh_rai_eff_grayscale.pdf"), width = 8, height = 4)

# Calculate the range of non-adherence and prevalence of RAI in simulations consistent with trial results. Define consistent as +/- 1% in loess smooths
hr_fit[hr >= 0.73 * 0.99 & hr <= 0.73 * 1.01, .(min_prop_non_adh = min(prop_non_adh), max_prop_non_adh = max(prop_non_adh), range = max(prop_non_adh) - min(prop_non_adh)), by = rr_ring] # Seems that the range of possible values increases with reduced efficacy. Think about why that might be - with lower efficacy, there's a wider range of possible proportion adherent.
hr_fit[hr >= 0.73 * 0.99 & hr <= 0.73 * 1.01, .(min_prev_ai = min(prev_ai), max_prev_ai =  max(prev_ai), range = max(prev_ai) - min(prev_ai)), by = rr_ring]

#### PER-EXPOSURE EFFECT ####
# Load data from simulations in which, among women engaging in RAI, the proportion of acts that are RAI systematically vary
load(paste0(input_dir, "/sim_results_prop_ai_dt.RData"))

dt_median_iqr <- as.data.table(dt %>% group_by(prop_ai, rr_ring, re_randomize) %>% summarise(median_per_exp_rr = median(per_exp_rr_ai_adh), lb = quantile(per_exp_rr_ai_adh, prob = 0.25), ub = quantile(per_exp_rr_ai_adh, prob = 0.75)))

dt_median_iqr[, eff_med_iqr := paste0(round(1 - median_per_exp_rr, 2) * 100, "% (IQR: ", round(1 - ub, 2) * 100, ", ", round(1 - lb, 2) * 100, ")")]

print(paste0("In simulations with 70% per-vaginal exposure efficacy and stratified re-randomization, the per-exposure effect of the ring among women engaged in RAI for 6.3% of her total acts was ", dt_median_iqr[prop_ai == 0.063 & rr_ring == 0.30 & re_randomize == "TRUE", eff_med_iqr], ", a ", dt_median_iqr[prop_ai == 0.063 & rr_ring == 0.30 & re_randomize == "TRUE", round(((1 - rr_ring) - (1 - median_per_exp_rr))/(1 - rr_ring), 2) * 100], "% reduction in effectiveness."))

print(paste0("The average per-exposure effect of the ring was highest among women engaged in RAI for 5% of their acts (", dt_median_iqr[prop_ai == 0.05 & rr_ring == 0.30 & re_randomize == "TRUE", eff_med_iqr], " relative to 70% per-vaginal exposure) and lowest for women for whom 30% of acts were RAI (", dt_median_iqr[prop_ai == 0.30 & rr_ring == 0.30 & re_randomize == "TRUE", eff_med_iqr], " relative to 70% per-vaginal exposure)"))

print(paste0("At the median frequency of RAI of one in every 16 acts, the ring reduced HIV-1 risk by approximately ", dt_median_iqr[prop_ai == 0.063 & re_randomize == "TRUE", .(paste0(min(round(1 - median_per_exp_rr, 2) * 100), "-", max(round(1 - median_per_exp_rr, 2) * 100)))], "% across variation in per-vaginal exposure efficacy of 65-75%, respectively"))

print(paste0("this risk reduction diminished with increasing RAI frequency to ", dt_median_iqr[prop_ai == 0.3 & re_randomize == "TRUE", .(paste0(min(round(1 - median_per_exp_rr, 2) * 100), "-", max(round(1 - median_per_exp_rr, 2) * 100)))], "% among women for whom one in every three acts is RAI"))

# Factor efficacy to plot from low to high
dt[, rr_ring := factor(rr_ring, levels = c(0.35, 0.30, 0.25))]

# Factor randomization so that balanced arms plot first
dt[, re_randomize := factor(re_randomize, levels = c("TRUE", "FALSE"))]

dt[, efficacy := 1 - per_exp_rr_ai_adh]

# Color figure
ggplot(data = dt[prop_ai %in% c(seq(0.05, 0.30, 0.05), 0.063)]) +
  geom_boxplot(aes(x = as.factor(prop_ai), y = efficacy, fill = re_randomize), outlier.shape = NA, coef = 0, show.legend = F) +
  geom_point(size = -1, aes(x = as.factor(prop_ai), y = efficacy, fill = re_randomize)) +
  geom_hline(aes(yintercept = 1 - as.numeric(as.character(rr_ring))), linetype = "dashed") +
  scale_fill_manual(labels = c("Re-balanced randomization", "As-randomized"), values = c(lighten("#B6A0C4", amount = 0.5), "#9670A5")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Proportion of acts that are receptive anal intercourse", y = "Overall HIV-1 risk reduction", fill = "") +
  scale_x_discrete(breaks = sort(c(seq(0, 0.30, 0.05), 0.063)), labels = c("0", "0.05", "0.063", "0.1", "0.15", "0.2", "0.25", "0.3")) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(. ~ as.factor(rr_ring), labeller = as_labeller(c("0.25" = "75% efficacy\nper vaginal exposure", "0.3" = "70% efficacy\nper vaginal exposure", "0.35" = "65% efficacy\nper vaginal exposure"))) +
  theme(text = element_text(family = "Times"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "grey", size = 0.1),
        panel.grid.minor.y = element_line(color = "grey", size = 0.05),
        legend.key = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))

ggsave(filename = paste0(fig_dir, "/per_act_rr_color.pdf"), width = 8, height = 4)

# Grayscale figure
ggplot(data = dt) +
  geom_boxplot(aes(x = as.factor(prop_ai), y = per_exp_rr_ai_adh, fill = re_randomize), outlier.shape = NA, coef = 0, show.legend = F) +
  geom_point(size = -1, aes(x = as.factor(prop_ai), y = per_exp_rr_ai_adh, fill = re_randomize)) +
  scale_fill_manual(labels = c("Baseline imbalance in arms", "Balanced arms"), values = c("gray55", "gray80")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Proportion of acts that are receptive anal intercourse", y = "Per-exposure ring effect", fill = "") +
  scale_x_discrete(breaks = sort(c(seq(0, 0.30, 0.05), 0.063)), labels = c("0", "0.05", "0.063", "0.1", "0.15", "0.2", "0.25", "0.3")) +
  facet_grid(. ~ as.factor(rr_ring), labeller = as_labeller(c("0.25" = "75% efficacy", "0.3" = "70% efficacy", "0.35" = "65% efficacy"))) +
  theme(text = element_text(family = "Times"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "grey", size = 0.1),
        panel.grid.minor = element_line(color = "grey", size = 0.05),
        legend.key = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))

ggsave(filename = paste0(fig_dir, "/per_act_rr_grayscale.pdf"), width = 8, height = 5)

# Dissertation figure #1: 6.3% proportion AI, 70% efficacy
dt_median_iqr[, `:=`(efficacy    = 1 - median_per_exp_rr,
                     efficacy_lb = 1 - lb,
                     efficacy_ub = 1 - ub)]

# Factor rr_ring to plot from low to high
dt_median_iqr[, rr_ring := factor(rr_ring, levels = c(0.35, 0.30, 0.25))]

temp_dt <-copy(dt_median_iqr)
temp_dt <- temp_dt[!(prop_ai == 0.063 & rr_ring == 0.30), `:=`(efficacy = NA, efficacy_lb = NA, efficacy_ub = NA)]

ggplot(data = temp_dt[re_randomize == "TRUE" & prop_ai %in% c(seq(0.05, 0.30, 0.05), 0.063)]) +
  geom_bar(aes(x = as.factor(prop_ai), y = efficacy), fill = "#B6A0C4", stat = "identity", color = "black") +
  geom_errorbar(aes(x = as.factor(prop_ai), ymin = efficacy_lb, ymax = efficacy_ub), width = 0.3) +
  geom_hline(aes(yintercept = 1 - as.numeric(as.character(rr_ring))), linetype = "dashed") +
  geom_hline(aes(yintercept = 0), linetype = "solid") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Proportion of acts that are anal intercourse", y = "Average per-exposure ring effect") +
  scale_x_discrete(breaks = sort(c(seq(0, 0.30, 0.05), 0.063)), labels = c("0", "0.05", "0.063", "0.1", "0.15", "0.2", "0.25", "0.3")) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(. ~ as.factor(rr_ring), labeller = as_labeller(c("0.25" = "75% efficacy\nper vaginal exposure", "0.3" = "70% efficacy\nper vaginal exposure", "0.35" = "65% efficacy\nper vaginal exposure"))) +
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "grey", size = 0.1),
        panel.grid.minor.y = element_line(color = "grey", size = 0.05),
        legend.key = element_blank(),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none")

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/Per_Act_RR_1.pdf", width = 9, height = 4)

# Dissertation figure #2: with stratified re-randomization only, all 70% efficacy
temp_dt <-copy(dt_median_iqr)
temp_dt <- temp_dt[!(rr_ring == 0.30), `:=`(efficacy = NA, efficacy_lb = NA, efficacy_ub = NA)]

ggplot(data = temp_dt[re_randomize == "TRUE" & prop_ai %in% c(seq(0.05, 0.30, 0.05), 0.063)]) +
  geom_bar(aes(x = as.factor(prop_ai), y = efficacy), fill = "#B6A0C4", stat = "identity", color = "black") +
  geom_errorbar(aes(x = as.factor(prop_ai), ymin = efficacy_lb, ymax = efficacy_ub), width = 0.3) +
  geom_hline(aes(yintercept = 1 - as.numeric(as.character(rr_ring))), linetype = "dashed") +
  geom_hline(aes(yintercept = 0), linetype = "solid") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Proportion of acts that are anal intercourse", y = "Average per-exposure ring effect") +
  scale_x_discrete(breaks = sort(c(seq(0, 0.30, 0.05), 0.063)), labels = c("0", "0.05", "0.063", "0.1", "0.15", "0.2", "0.25", "0.3")) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(. ~ as.factor(rr_ring), labeller = as_labeller(c("0.25" = "75% efficacy\nper vaginal exposure", "0.3" = "70% efficacy\nper vaginal exposure", "0.35" = "65% efficacy\nper vaginal exposure"))) +
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "grey", size = 0.1),
        panel.grid.minor.y = element_line(color = "grey", size = 0.05),
        legend.key = element_blank(),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none")

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/Per_Act_RR_2.pdf", width = 9, height = 4)

# Dissertation figure #3: with stratified re-randomization only, all
ggplot(data = dt_median_iqr[re_randomize == "TRUE" & prop_ai %in% c(seq(0.05, 0.30, 0.05), 0.063)]) +
  geom_bar(aes(x = as.factor(prop_ai), y = efficacy), fill = "#B6A0C4", stat = "identity", color = "black") +
  geom_errorbar(aes(x = as.factor(prop_ai), ymin = efficacy_lb, ymax = efficacy_ub), width = 0.3) +
  geom_hline(aes(yintercept = 1 - as.numeric(as.character(rr_ring))), linetype = "dashed") +
  geom_hline(aes(yintercept = 0), linetype = "solid") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Proportion of acts that are anal intercourse", y = "Average per-exposure ring effect") +
  scale_x_discrete(breaks = sort(c(seq(0, 0.30, 0.05), 0.063)), labels = c("0", "0.05", "0.063", "0.1", "0.15", "0.2", "0.25", "0.3")) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(. ~ as.factor(rr_ring), labeller = as_labeller(c("0.25" = "75% efficacy\nper vaginal exposure", "0.3" = "70% efficacy\nper vaginal exposure", "0.35" = "65% efficacy\nper vaginal exposure"))) +
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "grey", size = 0.1),
        panel.grid.minor.y = element_line(color = "grey", size = 0.05),
        legend.key = element_blank(),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none")

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/Per_Act_RR_3.pdf", width = 9, height = 4)

# Calculate mean improvement in efficacy with re-balanced arms
summ_dt <- dcast(data = dt_median_iqr, formula = prop_ai + rr_ring ~ re_randomize, value.var = c("median_per_exp_rr", "lb", "ub"))
summ_dt[, diff_rerandomize := median_per_exp_rr_FALSE - median_per_exp_rr_TRUE]
diff_rerandomize_summary <- summ_dt[, summary(diff_rerandomize)]
print(paste0("Stratified rerandomization resulted in an absolute difference in the estimated per-exposure ring effect of ", unname(round(diff_rerandomize_summary["Median"], 2)) * 100, "% (IQR: ", unname(round(diff_rerandomize_summary["1st Qu."], 2)) * 100, ", ", unname(round(diff_rerandomize_summary["3rd Qu."], 2)) * 100, ")"))

# Create figure showing the proportion of efficacy retained at varying proportions of acts that are RAI
dt[, prop_eff := round((1 - per_exp_rr_ai_adh)/(1 - rr_ring), 3)]

prop_eff_dt <- as.data.table(dt %>% group_by(prop_ai, re_randomize) %>% summarise(prop_eff_med = median(prop_eff), prop_eff_lb = quantile(prop_eff, prob = 0.25), prop_eff_ub = quantile(prop_eff, prob = 0.75)))

ggplot(data = prop_eff_dt) +
  geom_point(aes(x = prop_ai, y = prop_eff_med)) +
  geom_errorbar(aes(x = prop_ai, ymin = prop_eff_lb, ymax = prop_eff_ub), width = 0.01) +
  labs(x = "Proportion of acts that are receptive anal intercourse", y = "Proportion efficacy retained") +
  scale_x_continuous(breaks = seq(0.05, 0.30, 0.05), labels = formatC(seq(0.05, 0.30, 0.05), digits = 2, format = "f")) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(. ~ re_randomize, labeller = as_labeller(c("FALSE" = "Baseline imbalance in arms", "TRUE" = "Balanced arms"))) +
  theme(text = element_text(family = "Times", size = 13),
        strip.background = element_blank(),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 11),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "grey", size = 0.1),
        panel.grid.minor = element_line(color = "grey", size = 0.05),
        legend.key = element_blank())

ggsave(filename = paste0(fig_dir, "/prop_retained_efficacy.pdf"), width = 8, height = 4)

# Create character variable for median (IQR)
prop_eff_dt[, med_iqr := paste0(round(prop_eff_med, 2) * 100, "% (IQR: ", round(prop_eff_lb, 2) * 100, ", ", round(prop_eff_ub, 2) * 100, ")")]

print(paste0("This translates to ", dt[prop_ai == 0.05 & re_randomize == "TRUE", .(round(c(median(prop_eff), quantile(prop_eff, probs = c(0.25, 0.75))), 2))], " to ", dt[prop_ai == 0.30 & re_randomize == "TRUE", .(round(c(median(prop_eff), quantile(prop_eff, probs = c(0.25, 0.75))), 2))], " of the per-exposure effect of the ring among women with high adherence engaging in only RVI"))

print(paste0("Among women for women 7.8% of all vaginal acts are RAI, the ring retains ", prop_eff_dt[prop_ai == 0.078 & re_randomize == "TRUE", med_iqr], " of the full efficacy conferred to women engaged in only vaginal intercourse."))

print(paste0("The proportion of efficacy retained was highest for women engaged in RAI for 5% of acts ", prop_eff_dt[prop_ai == 0.05 & re_randomize == "TRUE", med_iqr], " and declined to ", prop_eff_dt[prop_ai == 0.30 & re_randomize == "TRUE", med_iqr], " retained among women who engaged in RAI for 30% of their acts."))

#### POPULATION PER-EXPOSURE EFFECT ####
# Subset dt_median_iqr to represent defined ranges, keeping only re-randomized simulations
dt_median_iqr_sub <- dt_median_iqr[prop_ai %in% c(0.01, 0.027, 0.063, 0.167, 0.30) & re_randomize == "TRUE"]

dt_median_iqr_sub_wide <- dcast(data = dt_median_iqr_sub, formula = prop_ai ~ rr_ring, value.var = "eff_med_iqr")

# Create data table with proportion of population engaging in RAI/VI only and among those engaged in RAI, the proportion of their acts that are RAI.
bar_ymin <- c(0, 82.0, 82.0 + 0.25 * 18, 82.0 + 2 * (0.25 * 18), 82.0 + 3 * (0.25 * 18), 82.0 + 3 * (0.25 * 18) + 0.124 * 18)
bar_ymax <- c(82.0, 82.0 + 0.25 * 18, 82.0 + 2 * (0.25 * 18), 82.0 + 3 * (0.25 * 18), 82.0 + 3 * (0.25 * 18) + 0.124 * 18, 100.0)

dt_figure <- data.table(ymin = c(0.0, 81.0, bar_ymin),
                        ymax = c(81.0, 100.0, bar_ymax),
                        act_type  = c("vi", "ai", "vi", "ai", "ai", "ai", "ai", "ai"),
                        prop_type = c("vi", "ai", "vi_1", "ai_1", "ai_2", "ai_3", "ai_4", "ai_5"))

# Factor prop_type
dt_figure[, prop_type := factor(prop_type, levels = c("vi", "ai", "vi_1", "ai_1", "ai_2", "ai_3", "ai_4", "ai_5"))]
dt_figure[, act_type := factor(act_type, levels = c("vi", "ai"))]

col_func <- colorRampPalette(colors = c(lighten("#B6A0C4", amount = 0.5), "#9670A5"))
cols <- col_func(n = 15)
purple_cols <- cols[seq(1, 12, length.out = 5)]

ggplot(data = dt_figure) +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 100), fill = "white") +
  geom_rect(aes(xmin = 2, xmax = 3, ymin = ymin, ymax = ymax, fill = act_type)) +
  geom_rect(aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = prop_type)) +
  coord_polar(theta = "y", start = 2.71) +
  labs(fill = "Average per-exposure\nring effect") +
  scale_fill_manual(values = c(lighten(col = "#9670A5", amount = 0), purple_cols, lighten(col = "#BED387", amount = 0), lighten(col = "#DFEAC0", amount = 0)), labels = c("Any receptive anal intercourse", "63-67%", "53-63%", "36-53%", "26-36%", "<26%", "Only receptive vaginal intercourse", "70%")) +
  theme(text = element_text(family = "Times"),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 11))

ggsave(filename = paste0(fig_dir, "/population_per_exp.pdf"), width = 6, height = 4)

# Dissertation figure #1: Inner circle only
ggplot(data = dt_figure[1:2, .SD, .SDcols = 1:3]) +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 100), fill = "white") +
  geom_rect(aes(xmin = 2, xmax = 3, ymin = ymin, ymax = ymax, fill = act_type)) +
  geom_rect(aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax), fill = "white") +
  coord_polar(theta = "y") +
  labs(fill = "Average per-exposure\nring effect") +
  scale_fill_manual(values = c("#BED387", "#9670A5")) +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

ggsave("/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/Donut_Chart_Prop_AI_1.pdf", width = 4, height = 4)

# Dissertation figure #2: Inner and outer circles
ggplot(data = dt_figure) +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 100), fill = "white") +
  geom_rect(aes(xmin = 2, xmax = 3, ymin = ymin, ymax = ymax, fill = act_type)) +
  geom_rect(aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = prop_type)) +
  coord_polar(theta = "y") +
  labs(fill = "Average per-exposure\nring effect") +
  scale_fill_manual(values = c(lighten(col = "#9670A5", amount = 0), purple_cols, lighten(col = "#BED387", amount = 0), lighten(col = "#DFEAC0", amount = 0))) +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

ggsave("/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/Donut_Chart_Prop_AI_2.pdf", width = 4, height = 4)

# Create separate figures for inner and outer circles to get separate legends for manually created figure
ggplot(data = dt_figure[1:2, .SD, .SDcols = 1:3]) +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 100), fill = "white") +
  geom_rect(aes(xmin = 2, xmax = 3, ymin = ymin, ymax = ymax, fill = act_type)) +
  geom_rect(aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax), fill = "white") +
  coord_polar(theta = "y") +
  labs(fill = "Average per-exposure ring effect") +
  scale_fill_manual(values = c("#BED387", "#9670A5"), labels = c("Only receptive vaginal intercourse", "Any receptive anal intercourse")) +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(family = "Times"),
        legend.text = element_text(size = 11))

ggsave(filename = paste0(fig_dir, "/population_per_exp_legend_1.pdf"), width = 6, height = 4)

ggplot(data = dt_figure[3:8, .SD, .SDcols = c(1, 2, 4)]) +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 100), fill = "white") +
  geom_rect(aes(xmin = 2, xmax = 3, ymin = ymin, ymax = ymax), fill = "white") +
  geom_rect(aes(xmin = 3, xmax = 4, ymin = ymin, ymax = ymax, fill = prop_type)) +
  coord_polar(theta = "y") +
  labs(fill = "Average per-exposure ring effect") +
  scale_fill_manual(values = c(lighten(col = "#DFEAC0", amount = 0), purple_cols), labels = c("70%", "63-67%", "53-63%", "36-53%", "26-36%", "<26%")) +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(family = "Times"),
        legend.text = element_text(size = 11))

ggsave(filename = paste0(fig_dir, "/population_per_exp_legend_2.pdf"), width = 6, height = 4)

#### SENSITIVITY ANALYSIS: PER-EXPOSURE EFFECT ####
# Load data from simulations in which, among women engaging in RAI, the proportion of acts that are RAI systematically vary
load(paste0(input_dir, "/sim_results_prop_ai_dt.RData"))

dt_median_iqr <- as.data.table(dt %>% group_by(prop_ai, rr_ring, re_randomize) %>% summarise(median_per_exp_rr = median(per_exp_rr_ai_adh), lb = quantile(per_exp_rr_ai_adh, prob = 0.25), ub = quantile(per_exp_rr_ai_adh, prob = 0.75)))

dt_median_iqr_primary <- copy(dt_median_iqr)
dt_primary <- copy(dt)

# Load data from sensitivity analyses
load(paste0(input_dir, "/sim_results_prop_ai_sens_dt.RData"))

dt_median_iqr <- as.data.table(dt %>% group_by(prop_ai, rr_ring, re_randomize, rr_ai_round) %>% summarise(median_per_exp_rr = median(per_exp_rr_ai_adh), lb = quantile(per_exp_rr_ai_adh, prob = 0.25), ub = quantile(per_exp_rr_ai_adh, prob = 0.75)))

dt_median_iqr_sens <- copy(dt_median_iqr)
dt_sens <- copy(dt)

dt_median_iqr_sens[, eff_med_iqr := paste0(round(1 - median_per_exp_rr, 2) * 100, "% (IQR: ", round(1 - ub, 2) * 100, ", ", round(1 - lb, 2) * 100, ")")]

print(paste0("In sensitivity analyses with 70% per-vaginal exposure efficacy and stratified re-randomization, the per-exposure effect of the ring among women engaged in RAI for 6.3% of her total acts was ", dt_median_iqr_sens[prop_ai == 0.063 & rr_ring == 0.30 & re_randomize == "TRUE" & rr_ai_round == 16, eff_med_iqr], ", a ", dt_median_iqr_sens[prop_ai == 0.063 & rr_ring == 0.30 & re_randomize == "TRUE" & rr_ai_round == 16, round(((1 - rr_ring) - (1 - median_per_exp_rr))/(1 - rr_ring), 2) * 100], "% reduction in effectiveness."))

print(paste0("In sensitivity analyses, the average per-exposure effect of the ring was highest among women engaged in RAI for 5% of their acts (", dt_median_iqr_sens[prop_ai == 0.05 & rr_ring == 0.30 & re_randomize == "TRUE" & rr_ai_round == 16, eff_med_iqr], " relative to 70% per-vaginal exposure) and lowest for women for whom 30% of acts were RAI (", dt_median_iqr_sens[prop_ai == 0.30 & rr_ring == 0.30 & re_randomize == "TRUE" & rr_ai_round == 16, eff_med_iqr], " relative to 70% per-vaginal exposure)"))

# Concatenate data for plotting
dt_primary[, rr_ai := "6.9"]

dt_rr_ai <- data.table(rr_ai = c(5.4, 5.8, 8.2, 10.2, 11.4, 11.7, 12.9, 13.9, 15.1, 15.8), rr_ai_round = c(5, 6, 8, 10:16))

dt_sens <- merge(x = dt_sens, y = dt_rr_ai, by = "rr_ai_round")
dt_sens <- dt_sens[, .(sim, prop_ai, rr_ring, re_randomize, per_exp_rr_ai_adh, rr_ai)]

dt <- rbindlist(l = list(dt_primary, dt_sens))

# Limit data to ring efficacy of 70% and prop_ai of 5-30%
dt <- dt[rr_ring == 0.3 & prop_ai %in% c(seq(0.05, 0.30, 0.05), 0.063) & re_randomize == "TRUE"]

# Factor rr_ai to plot in desired order
dt[, rr_ai := factor(rr_ai, levels = c("5.4", "5.8", "6.9", "8.2", "10.2", "11.4", "11.7", "12.9", "13.9", "15.1", "15.8"))]

dt[, efficacy := 1 - per_exp_rr_ai_adh]

col_func <- colorRampPalette(colors = c(lighten("#B6A0C4", amount = 0.2), "#DFEAC0"))
cols <- col_func(n = dt[, length(unique(rr_ai))])

# Color figure
ggplot(data = dt) +
  geom_boxplot(aes(x = as.factor(prop_ai), y = efficacy, fill = rr_ai), outlier.shape = NA, coef = 0, show.legend = F) +
  geom_point(size = -1, aes(x = as.factor(prop_ai), y = efficacy, fill = rr_ai)) +
  scale_fill_manual(values = cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Proportion of acts that are receptive anal intercourse", y = "Overall HIV-1 risk reduction", fill = expression(RR[AI])) +
  scale_x_discrete(breaks = sort(c(seq(0, 0.30, 0.05), 0.063)), labels = c("0", "0.05", "0.063", "0.10", "0.15", "0.20", "0.25", "0.30")) +
  scale_y_continuous(labels = scales::percent) +
  theme(text = element_text(family = "Times"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "grey", size = 0.1),
        panel.grid.minor = element_line(color = "grey", size = 0.05),
        legend.key = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 9.5),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right") +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))

ggsave(filename = paste0(fig_dir, "/per_act_rr_color_sens.pdf"), width = 8, height = 4)

# Dissertation figure alternative: bar chart
dt_median_iqr_figure <- as.data.table(dt %>% group_by(prop_ai, rr_ai) %>% summarise(median_eff = median(efficacy), lb = quantile(efficacy, prob = 0.25), ub = quantile(efficacy, prob = 0.75)))

ggplot(data = dt_median_iqr_figure) +
  geom_bar(aes(x = as.factor(prop_ai), y = median_eff, fill = rr_ai), stat = "identity", position = position_dodge(width = 0.9), color = "black", size = 0.3) +
  geom_errorbar(aes(x = as.factor(prop_ai), ymin = lb, ymax = ub, group = rr_ai), position = position_dodge(width = 0.9), width = 0.4, size = 0.2) +
  geom_hline(aes(yintercept = 0.70), linetype = "dashed", size = 0.2) +
  scale_fill_manual(values = cols) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Proportion of acts that are receptive anal intercourse", y = "Average per-exposure ring effect", fill = expression(RR[AI])) +
  scale_x_discrete(breaks = sort(c(seq(0, 0.30, 0.05), 0.063)), labels = c("0", "0.05", "0.063", "0.10", "0.15", "0.20", "0.25", "0.30")) +
  scale_y_continuous(labels = scales::percent) +
  theme(text = element_text(family = "sans"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = "grey", size = 0.1),
        panel.grid.minor.y = element_line(color = "grey", size = 0.05),
        legend.key = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 9.5),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right")

ggsave(filename = "/Volumes/GoogleDrive/My Drive/UW Epi Program/Dissertation/Dissertation oral defense/Figures/Per_Act_RR_Sens_Bar_Chart.pdf", width = 8, height = 4)

# Grayscale figure
ggplot(data = dt) +
  geom_boxplot(aes(x = as.factor(prop_ai), y = per_exp_rr_ai_adh, fill = analysis), outlier.shape = NA, coef = 0, show.legend = F) +
  geom_point(size = -1, aes(x = as.factor(prop_ai), y = per_exp_rr_ai_adh, fill = analysis)) +
  scale_fill_manual(values = c("gray55", "gray80")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Proportion of acts that are receptive anal intercourse", y = "Per-exposure ring effect", fill = "RR of RAI") +
  scale_x_discrete(breaks = sort(c(seq(0, 0.30, 0.05), 0.063)), labels = c("0", "0.05", "0.063", "0.10", "0.15", "0.20", "0.25", "0.30")) +
  theme(text = element_text(family = "Times"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "grey", size = 0.1),
        panel.grid.minor = element_line(color = "grey", size = 0.05),
        legend.key = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 9.5),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = "right") +
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5)))

ggsave(filename = paste0(fig_dir, "/per_act_rr_grayscale_sens.pdf"), width = 6, height = 3)

# Calculate mean reduction in efficacy with higher RR of RAI
summ_dt <- merge(x = dt_median_iqr_sens[rr_ring == 0.3, .(prop_ai, re_randomize, median_per_exp_rr)], y = dt_median_iqr_primary[rr_ring == 0.3, .(prop_ai, re_randomize, median_per_exp_rr)], by = c("prop_ai", "re_randomize"), all.x = T, all.y = F)
summ_dt[, percent_diff_rr_rai := (median_per_exp_rr.x - median_per_exp_rr.y)/median_per_exp_rr.y]
diff_rerandomize_summary <- summ_dt[, summary(percent_diff_rr_rai)]
print(paste0("Assuming a higher relative risk of RAI value of 14.3 resulted in an absolute difference in the estimated per-exposure ring effect of ", unname(round(diff_rerandomize_summary["Median"], 2)) * 100, "% (IQR: ", unname(round(diff_rerandomize_summary["1st Qu."], 2)) * 100, ", ", unname(round(diff_rerandomize_summary["3rd Qu."], 2)) * 100, ")"))

