## Purpose: Calculate proportion of women by site with plasma dapivirine > 95 pg/mL. Use average over time.
## Author:  Kathryn Peebles
## Date:    14 October 2019

# Attach packages
library(data.table)
library(dplyr)

# Load data
adh_dt <- as.data.table(read.csv(file = "/Volumes/GoogleDrive/My Drive/UW Epi Program/IS - EB/Analysis/plasma_dpv.csv"))
load("~/Documents/code/aspire/data-private/hiv.RData")
hiv_dt <- as.data.table(hiv)

# Format ptid
adh_dt[, ptid := as.numeric(gsub(pattern = "-", replacement = "", x = ptid))]
hiv_dt[, ptid := as.numeric(ptid)]

# Create dichotomous variable indicating adherent or not
adh_dt[, adh := ifelse(conc > 95, 1, 0)]

# Merge relevant variables on id
dt <- merge(x = adh_dt[, .(ptid, adh, visit)], y = hiv_dt[, .(ptid, site, arm)], by = "ptid", all = T)

# Remove observations from women not confirmed to be HIV-negative at enrollment and with at least one follow-up HIV test (i.e., limit to cohort of 2,614 women in mITT analyses). Also remove observations from women in the placebo arm.
dt <- dt[ptid %in% hiv_dt[, unique(ptid)] & arm == "Dapivirine"]

# There are 15 women with no adherence measures (have values of NA for adh and visit in adh_dt). Remove them from dt.
dt <- dt[!is.na(adh)]

# Format site variable to be concordant with formatting in f_dt
dt[, site := tstrsplit(x = site, split = ": ", keep = 2)]

prop_adh_dt <- as.data.table(dt %>% group_by(site) %>% summarise(prop_adh_obs = mean(adh)))

# For baseline simulations, use observed prevalence of adherence. In sensitivity analyses, use values of 50% to 100% adherence. Modify site prevalence in proportion to its observed prevalence in order to maintain the same relative adherence by site.
mean_adh_obs <- prop_adh_dt[, mean(prop_adh_obs)] # Value of 81% is slightly different than 82% estimated from dt[, mean(adh)] due to variation in person-time by site. Weighted average should be equal. However, because the difference is so small, will use an unweighted average.

# Increase or decrease adherence proportionally to observed adherence. For high levels of adherence (90% and 95%), this results in greater than 100% adherence.
prop_adh_dt[, `:=`(`50` = prop_adh_obs * (1 - (mean_adh_obs - 0.50)/mean_adh_obs),
                   `55` = prop_adh_obs * (1 - (mean_adh_obs - 0.55)/mean_adh_obs),
                   `60` = prop_adh_obs * (1 - (mean_adh_obs - 0.60)/mean_adh_obs),
                   `65` = prop_adh_obs * (1 - (mean_adh_obs - 0.65)/mean_adh_obs),
                   `70` = prop_adh_obs * (1 - (mean_adh_obs - 0.70)/mean_adh_obs),
                   `75` = prop_adh_obs * (1 - (mean_adh_obs - 0.75)/mean_adh_obs),
                   `80` = prop_adh_obs * (1 - (mean_adh_obs - 0.80)/mean_adh_obs),
                   `85` = prop_adh_obs * (1 - (mean_adh_obs - 0.85)/mean_adh_obs),
                   `90` = prop_adh_obs * (1 - (mean_adh_obs - 0.90)/mean_adh_obs),
                   `95` = prop_adh_obs * (1 - (mean_adh_obs - 0.95)/mean_adh_obs),
                   `100` = 1)]

# Adjust estimates of adherence greater than 100%.
prop_adh_dt[`90` > 1, `90` := 1]
mean_adh_90 <- prop_adh_dt[, mean(`90`)]
prop_adh_dt[`90` < 1, `90` := `90` * (1 + ((0.90 - mean_adh_90)/0.90))]

prop_adh_dt[`95` > 1, `95` := 1]
mean_adh_95 <- prop_adh_dt[, mean(`95`)]
prop_adh_dt[`95` < 1, `95` := `95` * (1 + ((0.96 - mean_adh_95)/0.96))] # Since we're increasing the mean only among the subset of women with adherence less than 100%, the relative proportion increase needs to be adjusted upward to result in overall mean adherence of 95%.
prop_adh_dt[`95` > 1, `95` := 1]

save(prop_adh_dt, file = "~/Documents/code/aspire/data-public/adh_dt.RData")
