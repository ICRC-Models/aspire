## Purpose: Create data table of study participants with fixed characteristics.
## Date:    9 March 2018
## Author:  Kathryn Peebles

# Data from Elizabeth Brown, email dated 2018-02-08. Note that there are 2,614 observations here, while Baeten, 2016 reports 2,629 women in the trial. 15 women did not report a partner in the previous 3 months, and so are excluded from model analyses.

library(data.table)

setwd("~/Documents/code/aspire/data-private")

source("~/Documents/code/aspire/fx/assign_m_hiv_rr.R")
load(paste0(getwd(), "/hiv.RData"))
dt <- as.data.table(hiv)

## 13 women are missing values for unprotected sex in last week. Impute values for them.
coef <- glm(unplast7base ~ basehivsp + baseSTD + noalc + age + married, data = hiv, family = binomial(link = "logit"))$coef
dt[is.na(unplast7base), log_odds_unp := coef[1] + 
     coef[2] * (basehivsp == "HIV negative") + 
     coef[3] * (basehivsp == "HIV positive") +
     coef[4] * (basehivsp == "participant does not know")  +
     coef[5] * baseSTD +
     coef[6] * noalc +
     coef[7] * age +
     coef[8] * married]
dt[is.na(unplast7base), prob_unp := exp(log_odds_unp)/(1 + exp(log_odds_unp))]
dt[is.na(unplast7base) & prob_unp < 0.50, unplast7base := TRUE]

f_dt <- dt[, .(site, country, arm, fu_days, nVI7base, age, basenumsp, basesyp, basetrr, basect, basegon, basebv, basehivsp, ponarv, married, unplast7base, enrolldt)]

f_dt[, id := 1:nrow(f_dt)]

f_dt[country == "Malawi",       country := "mal"]
f_dt[country == "South Africa", country := "sa"]
f_dt[country == "Uganda",       country := "uga"]
f_dt[country == "Zimbabwe",     country := "zim"]

f_dt[, site := tstrsplit(x = site, split = ": ", keep = 2)]

f_dt[, arm := ifelse(test = arm == "Placebo",
                     yes  = 0,
                     no   = 1)]

f_dt[, days_pre_adh_int := as.numeric(as.Date("2013-08-01") - enrolldt)]
f_dt[, days_pre_adh_int := ifelse(days_pre_adh_int < 0, 0, days_pre_adh_int)]
f_dt[, enrolldt := NULL]

f_dt[, f_age_cat := c("18-21", "22-26", "27-45")[findInterval(x = f_dt$age, vec = c(18, 22, 27))]] # age categories (18-21, 22-26, 27-45)

f_dt[, bv := as.numeric(basebv)]
f_dt[, n_sti := sum(basesyp, basetrr, basect, basegon), by = id]

f_dt[, condom_lweek := !unplast7base]

f_dt[, `:=`(basebv = NULL, basesyp = NULL, basetrr = NULL, basect = NULL, basegon = NULL, unplast7base = NULL)]

setkey(x = f_dt, country, age)

setnames(x = f_dt, old = c("age", "nVI7base", "basenumsp", "ponarv"), new = c("f_age", "n_acts", "n_part", "m_arv"))

# Set max number of baseline partners to 6 and minimum to 1.
f_dt[n_part > 6,  n_part := 6]
f_dt[n_part == 0, n_part := 1]

# Set max number of concurrent partners in simulation to n_part + 1. 
f_dt[, max_part := n_part + 1]

# If participant is married, set max_part and n_part equal to 1. In partner_change function, those with max_part == 1 will be uneligible for partner acquisition and dissolution
f_dt[married == T, `:=`(max_part = 1, n_part = 1)]

# Set max number of weekly acts to 35.
f_dt[n_acts > 35, n_acts := 35]

# Model runs in 30-day increments. Scale weekly number of acts to monthly number of acts.
f_dt[, n_acts := n_acts * 30/7]

# Set min number of monthly acts to 1.
f_dt[n_acts == 0, n_acts := 1]

# Number of acts is total acts. Convert to per-partner number of monthly acts. This assumes that acts are evenly distributed across partners, that per-partner number of acts at baseline is applicable to all follow-up points, and that woman-specific number of acts per partner does not vary across a woman's relationships.
f_dt[, pp_acts := n_acts/n_part]
f_dt[, n_acts := NULL]

# HIV status
f_dt[basehivsp == "" | basehivsp == "participant does not know", m_hiv := "unknown"]
f_dt[basehivsp == "HIV positive", m_hiv := "positive"]
f_dt[basehivsp == "HIV negative", m_hiv := "negative"]
f_dt[, basehivsp := NULL]

# Assign relative probability of having an HIV-positive male partner
f_dt[, m_hiv_rr := sapply(f_age, function(x) { assign_m_hiv_rr(age = x, l = 0.95) })]

setcolorder(x = f_dt, neworder = c("id", "country", "site", "f_age", "f_age_cat", "n_part", "max_part", "pp_acts", "married", "m_hiv_rr", "m_hiv", "m_arv", "n_sti", "bv", "condom_lweek", "arm", "fu_days", "days_pre_adh_int"))

save(f_dt, file = paste0(getwd(), "/f_dt.RData"))
