## Purpose: Create data table of study participants with fixed characteristics.
## Date:    9 March 2018
## Author:  Kathryn Peebles

# Data from Elizabeth Brown, email dated 2018-02-08. Note that there are 2,614 observations here, while Baeten, 2016 reports 2,629 women in the trial.

library(data.table)
library(lubridate)
library(zoo)

setwd("~/Documents/code/aspire/data-private")

# Load datasets
load(paste0(getwd(), "/hiv.RData")) # baseline variables
load(paste0(getwd(), "/ba.RData"))  # follow-up CRF
load(paste0(getwd(), "/bba.RData")) # baseline CRF
acasi <- as.data.table(read.csv(file = paste0(getwd(), "/ACASI.csv"))) # ACASI

# Merge baseline ACASI data with other baseline variables.
hiv <- as.data.table(hiv)
hiv[, ptid := as.numeric(ptid)]

bba <- as.data.table(bba)
bba[, ptid := as.numeric(gsub(pattern = "-", replacement = "", x = ptid))]

acasi <- as.data.table(acasi)
acasi[, ptid := gsub(pattern = "-", replacement = "", x = ptid)]
acasi[, ptid := as.numeric(gsub(pattern = " ", replacement = "", x = ptid))]

dt <- merge(x = hiv,  y = acasi[survey == 1, .(ptid, QEXCH)], by = "ptid", all.x = T)
dt <- merge(x = dt,   y = bba[, .(ptid, BBA3mvsn, numsp)])

# Create baseline variables to be used in prediction logit models for missing data and for adherence. For adherence model, values are from predictive model from Elizabeth Brown, but I've used coding from predictive models from Betsy (see document aspiresummPlasmaOnlyFINAL.pdf, page 2, for coding notes). Betsy's report lists coding as, e.g., if var = 1, new var = 1, 0 otherwise, so I assume missing values are coded as 0 below.
dt[, `:=`(b_edu          = ifelse(edu == "Secondary or above", 1, 0),
          b_noalc        = ifelse(noalc, 1, 0),
          b_married      = ifelse(married, 1, 0),
          b_npart_eq_two = ifelse(npart == "(1,2]", 1, 0),
          b_npart_gt_two = ifelse(npart == "(2,100]", 1, 0),
          b_sti          = ifelse(baseSTD, 1, 0),
          b_trr          = ifelse(basetrr, 1, 0),
          b_ct           = ifelse(basect, 1, 0),
          b_gon          = ifelse(basegon, 1, 0),
          b_depo         = ifelse(bdepo, 1, 0),
          b_neten        = ifelse(bneten, 1, 0),
          b_circ         = ifelse(circ, 1, 0),
          m_arv          = ifelse(ponarv, 1, 0),
          b_enroll_Apr1  = ifelse(enrolldt > as.Date("2013-04-01"), 1, 0),
          b_pknow        = ifelse(pknow, 1, ifelse(basenumsp == 0, 1, 0)),
          b_unpsex       = ifelse(is.na(unplast7base), NA_real_, ifelse(unplast7base, 1, 0)),
          b_trans_sex    = ifelse(QEXCH == "                        " | QEXCH == "Skipped by participant  ", NA_real_, ifelse(QEXCH == "No                      ", 0, 1)),
          b_n_part       = ifelse(numsp == 0, 1, ifelse(numsp > 6, 6, numsp)))] # Set max number of baseline partners to 6 and minimum to 1 for use in partner_change function.

# Create empty follow-up dataset to fill in monthly visits. 
visits_dt <- as.data.table(expand.grid(ptid = unique(dt$ptid), visit = 0:34))
setkey(visits_dt, ptid, visit)

# Bring in time on-study to estimate the number of visits and the interval between visits
dt[, visit := 0]
visits_dt <- merge(x = visits_dt, y = dt[, .(ptid, visit, enrolldt, censordt, fu_days)], by = c("ptid", "visit"), all = T)
visits_dt[, n_visits := round(fu_days/30)]
visits_dt[, n_visits := ifelse(n_visits == 0, 1, n_visits)]
visits_dt[, monthly_visit_interval := fu_days/n_visits]

# Fill in visit dates for all visits beyond enrollment visit
visits_dt[, monthly_visit_interval := na.locf(monthly_visit_interval), by = ptid]
visits_dt[, enrolldt := na.locf(enrolldt), by = ptid]
visits_dt[visit == 0, visit_date := enrolldt]
visits_dt[visit != 0, visit_date := enrolldt + (monthly_visit_interval * visit)]

visits_dt[, censor_dt := na.locf(as.Date(censordt, format = "%d%b%Y")), by = ptid]
visits_dt[, n_visits := na.locf(n_visits), by = ptid]

# Check that each participant's last visit date and her censor date are aligned
visits_dt[which(censor_dt == visit_date), any(visit != n_visits)] # Should be F
visits_dt[, fu_days := na.locf(fu_days), by = ptid]
all(visits_dt[, any((visit_date - enrolldt) == fu_days), by = ptid]$V1) # Should be T

# Create indicator for on_study
visits_dt[, on_study := as.numeric(visit_date <= censor_dt)]

# Merge data table with all monthly visit dates with baseline variables. Carry forward values of baseline variables.
dt <- merge(x = dt[, .SD, .SDcols = which(names(dt) != "visit")], y = visits_dt[, .(ptid, visit, visit_date, on_study)], by = "ptid", all.y = T)

# Merge ACASI data for baseline and Q1. (Add visit variable)
acasi[, visit := ifelse(survey == 1, 0, ifelse(survey == 2, 3, NA_real_))]
dt <- merge(x = dt, y = acasi[visit %in% c(0, 3), .(ptid, visit, QANAL)], by = c("ptid", "visit"), all.x = T)

# Create clean anal acts variable
dt[, anal_n_acts := ifelse(QANAL == "Skipped by participant  ", NA_real_, as.numeric(as.character(QANAL)))]

dt[country == "Malawi",       country := "mal"]
dt[country == "South Africa", country := "sa"]
dt[country == "Uganda",       country := "uga"]
dt[country == "Zimbabwe",     country := "zim"]

dt[, site := tstrsplit(x = site, split = ": ", keep = 2)]

dt[, arm := ifelse(arm == "Placebo", 0, 1)]

# Create indicator variables for enroll date after April 1, 2013 and visit date of previous visit being after August 1, 2013 (adherence procedures changed on these dates, with subsequent improvements in adherence).
dt[, enrolldt_after_Apr12013 := as.numeric(enrolldt > as.Date("2013-04-01", format = "%Y-%m-%d"))]
dt[, prev_visit_date_after_Aug12013 := as.numeric(shift(visit_date) > as.Date("2013-08-01", format = "%Y-%m-%d")), by = ptid]

dt[, b_f_age_cat := c("18-21", "22-26", "27-45")[findInterval(x = dt$age, vec = c(18, 22, 27))]] # age categories (18-21, 22-26, 27-45)

dt[, b_bv := as.numeric(basebv)]
dt[, b_n_sti := sum(basesyp, basetrr, basect, basegon), by = c("ptid", "visit")]

dt[, b_condom_lweek := as.numeric(!unplast7base)]

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

# Fill in ages for visits subsequent to baseline.
setnames(x = dt, old = "age", new = "b_f_age")
dt[, f_age := b_f_age + (30/365) * visit]

ids <- data.table(ptid = unique(dt$ptid), id = 1:length(unique(dt$ptid)))
dt <- merge(x = dt, y = ids, by = "ptid", all.x = T)

f_dt <- dt[, .(id, visit, enrolldt, visit_date, on_study, site, country, arm, b_sti, b_trr, b_ct, b_gon, b_circ, b_noalc, b_f_age, f_age, b_pknow, m_arv, b_depo, b_neten, b_edu, b_married, b_npart_eq_two, b_npart_gt_two, b_trans_sex, b_n_part, anal_n_acts, enrolldt_after_Apr12013, prev_visit_date_after_Aug12013, b_f_age_cat, b_bv, b_n_sti, b_condom_lweek, max_part, pp_vi_acts_avg, m_hiv)]

setcolorder(x = f_dt, neworder = c("id", "visit", "visit_date", "on_study", "enrolldt", "country", "site", "b_f_age", "b_f_age_cat", "f_age", "b_n_part", "max_part", "b_npart_eq_two", "b_npart_gt_two", "pp_vi_acts_avg", "anal_n_acts", "b_trans_sex", "b_married", "b_edu", "b_noalc", "b_pknow", "m_hiv", "m_arv", "b_n_sti", "b_sti", "b_trr", "b_ct", "b_gon", "b_bv", "b_circ", "b_depo", "b_neten", "b_condom_lweek", "arm", "enrolldt_after_Apr12013", "prev_visit_date_after_Aug12013"))

save(f_dt, file = paste0(getwd(), "/f_dt.RData"))
