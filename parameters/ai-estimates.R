## Purpose: Calculate proportion of women reporting any anal sex by site.
## Author:  Kathryn Peebles
## Date:    14 October 2019

# Attach packages
library(data.table)
library(dplyr)

# Load data
acasi_dt <- as.data.table(read.csv(file = "/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Data/ACASI.csv", na.strings = ""))
load("~/Documents/code/aspire/data-private/hiv.RData")
hiv_dt <- as.data.table(hiv)

# Format ptid
acasi_dt[, ptid := gsub(pattern = "-", replacement = "", x = ptid)]
acasi_dt[, ptid := gsub(pattern = " ", replacement = "", x = ptid)]
acasi_dt[, ptid := as.numeric(ptid)]

hiv_dt[, ptid := as.numeric(ptid)]

# Keep data from surveys 1 and 2
acasi_dt <- acasi_dt[survey == 1 | survey == 2]

# Keep observations from women confirmed to be HIV-negative at enrollment and with at least one follow-up HIV test
acasi_dt <- acasi_dt[ptid %in% hiv_dt[, unique(ptid)], .(ptid, survey, QANAL)]

# Create numeric number of anal acts variable
acasi_dt[, n_anal_acts := ifelse(QANAL == "Skipped by participant  ", NA, as.numeric(as.character(QANAL)))]

# Convert data from long to wide
acasi_dt_wide <- dcast(data = acasi_dt, formula = ptid ~ survey, value.var = "n_anal_acts")

# Merge with hiv_dt
dt <- merge(x = hiv_dt[, .(ptid, fu_days, site)], y = acasi_dt_wide, by = "ptid", all = T)

# Create variable for any anal sex at either enrollment or month 3. Any participant who reports "no" at either M0 or M3 but has missing data for the other month hasa missing value, since they could have engaged in AI at the missing month. Those who report "no" at baseline and have less than three months of follow-up have a value of no AI (i.e., no AI reported during the study).
dt[, anal_any_m0_m3 := ifelse(`1` > 0 | `2` > 0, 1,
                              ifelse((`1` == 0 & `2` == 0) | (`1` == 0 & fu_days < 90), 0, NA_real_))]
dt[`1` == 0 & is.na(`2`) & fu_days < 90, anal_any_m0_m3 := 0]

# Check coding
dt[, table(`1`, `2`, anal_any_m0_m3, useNA = 'always')]

# Format site variable to be concordant with formatting in f_dt
dt[, site := tstrsplit(x = site, split = ": ", keep = 2)]

# Estimate RAI prevalence by site
prop_ai_dt <- as.data.table(dt %>% group_by(site) %>% summarise(prop_ai_obs = mean(anal_any_m0_m3, na.rm = T)))

# For baseline simulations, use observed prevalence of adherence. In sensitivity analyses, use values of 50% to 100% adherence. Modify site prevalence in proportion to its observed prevalence in order to maintain the same relative adherence by site.
mean_ai_obs <- prop_ai_dt[, mean(prop_ai_obs)] # Value of 21% is slightly different than 19% estimated from dt[, mean(anal_any_m0_m3 == 1, na.rm = T)] due to variation in number of people by site. Weighted average should be equal. However, because the difference is so small, will use an unweighted average. In simulations, we output the simulated proportion of women engaging in AI, so we can center median values of interest at the corresponding median prevalence of AI on the x-y plane.

# Increase or decrease adherence proportionally to observed adherence. For high levels of adherence (90% and 95%), this results in greater than 100% adherence.
prop_ai_dt[, `:=`(`0`   = 0,
                  `5`   = prop_ai_obs * (1 - (mean_ai_obs - 0.05)/mean_ai_obs),
                  `10`  = prop_ai_obs * (1 - (mean_ai_obs - 0.10)/mean_ai_obs),
                  `15`  = prop_ai_obs * (1 - (mean_ai_obs - 0.15)/mean_ai_obs),
                  `20`  = prop_ai_obs * (1 - (mean_ai_obs - 0.20)/mean_ai_obs),
                  `25`  = prop_ai_obs * (1 - (mean_ai_obs - 0.25)/mean_ai_obs),
                  `30`  = prop_ai_obs * (1 - (mean_ai_obs - 0.30)/mean_ai_obs),
                  `100` = 1)]

save(prop_ai_dt, file = "~/Documents/code/aspire/data-public/ai_dt.RData")

