## Purpose: Calculate proportion of women by site with plasma dapivirine > 95 pg/mL. Use average over time.
## Author:  Kathryn Peebles
## Date:    14 October 2019

# Attach packages
library(data.table)
library(dplyr)
library(gee)

# Load data
adh_dt <- as.data.table(read.csv(file = "/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Data/plasma_dpv.csv"))
load("~/Documents/code/aspire/data-private/hiv.RData")
hiv_dt <- as.data.table(hiv)

# Format ptid
adh_dt[, ptid := as.numeric(gsub(pattern = "-", replacement = "", x = ptid))]
hiv_dt[, ptid := as.numeric(ptid)]

# Create dichotomous variable indicating adherent or not
adh_dt[, adh := ifelse(conc > 95, 1, 0)]

# Create dichotomous variable indicating if sample was obtained after Sept. 1, 2013 (adherence interventions were implemented on 8/1/2013 and visits were monthly, so we wouldn't expect an effect of the adherence interventions to show up until at least 30 days later).
adh_dt[, collection_dt := as.Date(x = colldt, format = "%d-%b-%y")]
adh_dt[, post_adh_intervention := ifelse(collection_dt >= as.Date("2013-09-01", format = "%Y-%m-%d"), 1, 0)]

# Merge relevant variables on id
dt <- merge(x = adh_dt[, .(ptid, adh, visit, post_adh_intervention)], y = hiv_dt[, .(ptid, site, arm, age)], by = "ptid", all = T)

# Remove observations from women not confirmed to be HIV-negative at enrollment and with at least one follow-up HIV test (i.e., limit to cohort of 2,614 women in mITT analyses). Also remove observations from women in the placebo arm.
dt <- dt[ptid %in% hiv_dt[, unique(ptid)] & arm == "Dapivirine"]

# There are 15 women with no adherence measures (have values of NA for adh and visit in adh_dt). Remove them from dt.
dt <- dt[!is.na(adh)]

# Format site variable to be concordant with formatting in f_dt
dt[, site := tstrsplit(x = site, split = ": ", keep = 2)]

# Create age category variable
dt[, age_cat := ifelse(age <= 26, "18-26", "27-45")]

model_1 <- gee(formula   = adh ~ post_adh_intervention + site + age_cat,
               family    = poisson,
               na.action = na.omit,
               id        = ptid,
               corstr    = "exchangeable",
               data      = dt[site != "Blantyre" & site != "Lilongwe"])

rel_incr_adh <- unname(exp(model_1$coefficients["post_adh_intervention"]))

prop_adh_dt <- as.data.table(dt %>% group_by(site, age_cat) %>% summarise(prop_adh_obs = mean(adh)))

# Estimate the probability of adherence following adherence interventions as pre-adherence interventions adherence * 1.12. Both Blantyre and Lilongwe began enrollment only after adherence interventions were in place and so have the same adherence pre- and post-adherence interventions.
prop_adh_dt[site != "Blantyre" & site != "Lilongwe", adh_post_intervention := prop_adh_obs * rel_incr_adh]
prop_adh_dt[site == "Blantyre" | site == "Lilongwe", adh_post_intervention := prop_adh_obs]

prop_adh_dt <- melt(prop_adh_dt, id.vars = c("site", "age_cat"), measure.vars = c("prop_adh_obs", "adh_post_intervention"))

prop_adh_dt[, post_adh_intervention := ifelse(variable == "prop_adh_obs", 0, 1)]
setnames(prop_adh_dt, old = "value", new = "prop_adh_obs")

prop_adh_dt <- prop_adh_dt[, .(site, age_cat, post_adh_intervention, prop_adh_obs)]

# # In the adherence assignment function, should merge values for post adherence interventions with the period pre-intervention in order to not have missing or dropped observations (these observations are all pre-enrollment and are dropped during censoring).
# mal_dt <- prop_adh_dt[site == "Blantyre" | site == "Lilongwe"]
# mal_dt[, post_adh_intervention := 0]
# 
# prop_adh_dt <- rbindlist(l = list(prop_adh_dt, mal_dt))

setkey(prop_adh_dt, site, age_cat, post_adh_intervention)

prop_adh_dt[prop_adh_obs > 1, prop_adh_obs := 1]

# For baseline simulations, use observed prevalence of adherence. In sensitivity analyses, use values of 50% to 100% adherence. Modify site prevalence in proportion to its observed prevalence in order to maintain the same relative adherence by site.
mean_adh_obs <- prop_adh_dt[, mean(prop_adh_obs)] # Value of 85% is higher than 82% estimated from dt[, mean(adh)] due to variation in person-time by site and possibly due to modeling adherence in the second year. Average weighted by person-time should be about equal. However, because the difference is relatively small, will use an unweighted average to re-distribute adherence.

# Increase or decrease adherence proportionally to observed adherence. For high levels of adherence (80% to 95%), this results in greater than 100% adherence.
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
prop_adh_dt[`85` > 1, `85` := 1]
mean_adh_85 <- prop_adh_dt[, mean(`85`)]
prop_adh_dt[`85` < 1, `85` := `85` * (1 + ((0.85 - mean_adh_85)/0.85))]
prop_adh_dt[`85` > 1, `85` := 1]

prop_adh_dt[`90` > 1, `90` := 1]
mean_adh_90 <- prop_adh_dt[, mean(`90`)]
prop_adh_dt[`90` < 1, `90` := `90` * (1 + ((0.9025 - mean_adh_90)/0.9025))] # Since we're increasing the mean only among the subset of women with adherence less than 100%, the relative proportion increase needs to be adjusted upward to result in overall mean adherence of 90%.
prop_adh_dt[`90` > 1, `90` := 1]

prop_adh_dt[`95` > 1, `95` := 1]
mean_adh_95 <- prop_adh_dt[, mean(`95`)]
prop_adh_dt[`95` < 1, `95` := `95` * (1 + ((0.97 - mean_adh_95)/0.97))] # Since we're increasing the mean only among the subset of women with adherence less than 100%, the relative proportion increase needs to be adjusted upward to result in overall mean adherence of 95%.
prop_adh_dt[`95` > 1, `95` := 1]

save(prop_adh_dt, file = "~/Documents/code/aspire/data-public/adh_dt.RData")
