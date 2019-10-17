## Purpose: Create data table of study participants with fixed characteristics.
## Date:    9 March 2018
## Author:  Kathryn Peebles

# Attach packages
library(data.table)
library(dplyr)
library(lubridate)
library(zoo)

setwd("~/Documents/code/aspire/data-private")

# Load datasets
load(paste0(getwd(), "/hiv.RData")) # baseline variables
load(paste0(getwd(), "/ba.RData"))  # follow-up CRF
load(paste0(getwd(), "/bba.RData")) # baseline CRF

# Merge baseline ACASI data with other baseline variables.
hiv <- as.data.table(hiv)
hiv[, ptid := as.numeric(ptid)]

bba <- as.data.table(bba)
bba[, ptid := as.numeric(gsub(pattern = "-", replacement = "", x = ptid))]

dt <- merge(x = hiv, y = bba[, .(ptid, BBA3mvsn, BBAlvsc, numsp)], by = "ptid", all.x = T)

# Create baseline variables to be used in prediction logit models for missing data.
dt[, `:=`(b_noalc        = ifelse(noalc, 1, 0),
          b_married      = ifelse(married, 1, 0),
          b_sti          = ifelse(baseSTD, 1, 0),
          m_arv          = ifelse(ponarv, 1, 0),
          b_n_part       = ifelse(numsp == 0, 1, ifelse(numsp > 6, 6, numsp)))] # Set max number of baseline partners to 6 and minimum to 1 for use in partner_change function.

# Create empty follow-up dataset to fill in monthly visits. Allow up to 150 visits.
visits_dt <- as.data.table(expand.grid(ptid = unique(dt$ptid), visit = 0:150))
setkey(visits_dt, ptid, visit)

# Merge date of enrollment, censor status, censor date, and termination reason
visits_dt <- merge(x = dt[, .(ptid, enrolldt, censor, censordt = as.Date(censordt, format = "%d%B%Y"), termrsn)], y = visits_dt, by = "ptid", all = T)
setkey(visits_dt, ptid, visit)

# Add variable for model calendar time. Calendar time starts from the date of first study enrollment and proceeds in monthly intervals (= 365.25 days per year/12 months per year).
visits_dt[, first_enroll_dt := min(enrolldt)]
visits_dt[visit == 0, model_calendar_time := first_enroll_dt]
visits_dt[visit != 0, model_calendar_time := first_enroll_dt + ((365.25/12) * visit)]

# Create variable to indicate if model calendar time is on or after date of enrollment (by each woman)
visits_dt[, post_enrollment := ifelse(model_calendar_time >= enrolldt, 1, 0)]

# Create variable to indicate if model calendar time is after the date of early study termination (by woman). For women who ended the study per schedule, this indicator is 0 for all visits. For women who seroconverted, this indicator is 0 for all visits (women may have been more likely to terminate study participation following seroconversion, but were on study until seroconversion. We will treat this as an indication of willingness to continue participation in the absence of seroconversion).
visits_dt[, early_term_censor := ifelse(termrsn != "Sched. exit visit" & censor == 0, as.numeric(model_calendar_time > censordt), 0)]

# Check coding
visits_dt[censor == 1, all(early_term_censor == 0)] # Should be true
visits_dt[termrsn == "Sched. exit visit", all(early_term_censor == 0)] # Should be true
# View(visits_dt[termrsn != "Sched. exit visit" & censor == 0, .(ptid, visit, enrolldt, termrsn, censor, censordt, model_calendar_time, early_term_censor, post_enrollment)])

# There are three participants for whom the coding above assigns them 0 follow-up time. Manually adjust to ensure a minimum of 30 days of follow-up.
# View(visits_dt[, all(early_term_censor == post_enrollment), by = ptid][V1 == T])
visits_dt[ptid == 306300207 & visit == 11, post_enrollment := 1]
visits_dt[ptid == 312402776 & visit == 17, post_enrollment := 1]
visits_dt[ptid == 318303968 & visit == 21, post_enrollment := 1]

# Merge data table with all monthly visit dates with baseline variables. Carry forward values of baseline variables.
dt <- merge(x = dt[, .SD, .SDcols = which(names(dt) != "visit")], y = visits_dt[, .(ptid, visit, post_enrollment, early_term_censor)], by = "ptid", all.y = T)

dt[country == "Malawi",       country := "mal"]
dt[country == "South Africa", country := "sa"]
dt[country == "Uganda",       country := "uga"]
dt[country == "Zimbabwe",     country := "zim"]

dt[, site := tstrsplit(x = site, split = ": ", keep = 2)]

dt[, arm := ifelse(arm == "Placebo", 0, 1)]

# Create age variables
dt[, f_age := age]
dt[, f_age_cat := c("18-21", "22-26", "27-45")[findInterval(x = dt$age, vec = c(18, 22, 27))]] # age categories (18-21, 22-26, 27-45)

dt[, b_bv := as.numeric(basebv)]
dt[, b_n_sti := sum(basesyp, basetrr, basect, basegon), by = c("ptid", "visit")]

dt[, b_condom_lweek := as.numeric(!unplast7base)]

# Create country-specific condom use proportion as the median proportion of women reporting condom use at last vaginal sex by country
dt_cond <- dt[visit == 0, .(country, BBAlvsc)]
dt_cond[, condom := ifelse(BBAlvsc == "NA", NA,
                           ifelse(BBAlvsc == "neither", 0, 1))]

dt_cond <- as.data.table(dt_cond %>% group_by(country) %>% summarise(n_condom = sum(condom, na.rm = T), n_country = sum(!is.na(condom))))
dt_cond[, prop_condom := round(n_condom/n_country, 3)]

dt <- merge(x = dt, y = dt_cond[, .(country, prop_condom)], by = "country", all.x = T)

# Set max number of concurrent partners in simulation to b_n_part + 1. 
dt[, max_part := b_n_part + 1]

# If participant is married, set max_part and b_n_part equal to 1. In partner_change function, those with max_part == 1 will be ineligible for partner acquisition and dissolution
dt[married == T, `:=`(max_part = 1, b_n_part = 1)]

# Create monthly average acts variables for vaginal sex acts. Women reported number of vaginal sex acts at baseline in the previous 3 months, so divide by 3 to estimate monthly average. This assumes that number of baseline vaginal acts is representative of number of vaginal acts throughout study participation.
dt[, n_vi_monthly_avg := round(BBA3mvsn/3, 2)]

# Set min number of monthly acts to 1/3 (i.e., women are assumed to engage in a minimum of 1 sex act every 3 months).
dt[n_vi_monthly_avg == 0, n_vi_monthly_avg := round(1/3, 2)]

# Number of acts is total acts. Convert to per-partner number of monthly acts. This assumes that acts are evenly distributed across partners, that per-partner number of acts at baseline is applicable to all follow-up points, and that woman-specific number of acts per partner does not vary across a woman's relationships.
dt[, pp_vi_acts_avg := round(n_vi_monthly_avg/b_n_part, 2)]

# Male partner HIV status
dt[basehivsp == "" | basehivsp == "participant does not know", m_hiv := "unknown"]
dt[basehivsp == "HIV positive", m_hiv := "positive"]
dt[basehivsp == "HIV negative", m_hiv := "negative"]

ids <- data.table(ptid = unique(dt$ptid), id = 1:length(unique(dt$ptid)))
dt <- merge(x = dt, y = ids, by = "ptid", all.x = T)

f_dt <- dt[, .(id, visit, post_enrollment, early_term_censor, site, country, arm, b_sti, b_noalc, f_age, m_arv, b_married, b_n_part, f_age_cat, b_bv, b_n_sti, b_condom_lweek, max_part, pp_vi_acts_avg, m_hiv, prop_condom)]

setcolorder(x = f_dt, neworder = c("id", "visit", "post_enrollment", "early_term_censor", "country", "site", "f_age", "f_age_cat", "b_n_part", "max_part", "pp_vi_acts_avg", "b_married", "b_noalc", "m_hiv", "m_arv", "b_n_sti", "b_sti", "b_bv", "b_condom_lweek", "prop_condom", "arm"))

save(f_dt, file = paste0(getwd(), "/f_dt.RData"))
