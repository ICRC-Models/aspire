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
preg_dt <- as.data.table(read.csv(file = paste0(getwd(), "/preg_data.csv")))
load("/Volumes/GoogleDrive/My Drive/UW Epi Program/ASPIRE/Data/dt_fp_baseline.RData")

# Identify women who became pregnant. Format ptid and pregnancy results
preg_dt[, ptid := gsub(pattern = "-", replacement = "", x = ptid)]
preg_dt[, ptid := as.numeric(gsub(pattern = "'", replacement = "", x = ptid))]
preg_dt[, pregnant := ifelse(MLRhcg == "'               '" | MLRhcg == "'Not Available  '", NA_real_,
                             ifelse(MLRhcg == "'Negative       '", 0, 1))]

preg_dt[, any_preg := as.numeric(any(pregnant == 1, na.rm = T)), by = ptid]

preg_dt <- unique(preg_dt[, .(ptid, any_preg)])

# Merge baseline ACASI data with other baseline variables.
hiv <- as.data.table(hiv)
hiv[, ptid := as.numeric(ptid)]

bba <- as.data.table(bba)
bba[, ptid := as.numeric(gsub(pattern = "-", replacement = "", x = ptid))]

# Recode variables in bba to be used in adherence prediction models: same primary partner in last three months (yes or no primary partner vs. no), primary partner aware of study participation (yes vs. no/not sure/no primary partner), participant worried about having ring inside every day for at least a year (not at all vs. somewhat/very worried)
bba[, same_psp := ifelse(BBApsp == "N" | BBApsp == "No" | BBA3mpsp == "Yes", "Yes or no primary partner", "No")]
bba[, table(BBApsp, BBA3mpsp, same_psp, useNA = "always")]

bba[, pknow_study := ifelse(BBApsp == "N" | BBApsp == "No" | BBAstdsp == "No" | BBAstdsp == "Not Sure", "No/Not sure/No primary partner", "Yes")]
bba[, table(BBApsp, BBAstdsp, pknow_study, useNA = "always")]

bba[, ring_worries := ifelse(BBAvgrwr == "" | BBAvgrwr == "not at all worried", "Not at all worried/missing", "Somewhat or very worried")] # Two missing responses. I've grouped them with "not at all worried" to avoid missingness. Sufficiently small number that this shouldn't matter, particularly as this is a prediction, rather than inferential, model.
bba[, table(ring_worries, BBAvgrwr, useNA = 'always')]

# Create variable for OCP vs. other methods
dt_fp_baseline[, bl_ocp := ifelse(bl_method == "COC", 1, 0)]

dt <- merge(x = hiv, y = bba[, .(ptid, BBA3mvsn, BBAlvsc, numsp, same_psp, pknow_study, ring_worries)], by = "ptid", all.x = T)
dt <- merge(x = dt,  y = preg_dt, by = "ptid", all.x = T)
dt <- merge(x = dt, y = dt_fp_baseline[, .(ptid, bl_ocp)], all.x = T)

# Create baseline variables to be used in prediction logit models for missing data.
dt[, `:=`(b_noalc        = ifelse(noalc, 1, 0),
          b_married      = ifelse(married, 1, 0),
          b_sti          = ifelse(baseSTD, 1, 0),
          m_arv          = ifelse(ponarv, 1, 0),
          b_n_part       = ifelse(numsp == 0, 1, ifelse(numsp > 6, 6, numsp)))] # Set max number of baseline partners to 6 and minimum to 1 for use in partner_change function.

# Create empty follow-up dataset to fill in monthly visits. Allow up to 60 visits (approximately twice the number of maximum observed visits).
visits_dt <- as.data.table(expand.grid(ptid = unique(dt$ptid), visit = 0:60))
setkey(visits_dt, ptid, visit)

# Merge date of enrollment, censor status, censor date, and termination reason
visits_dt <- merge(x = dt[, .(ptid, enrolldt, censor, censordt = as.Date(censordt, format = "%d%B%Y"), termrsn, termdays)], y = visits_dt, by = "ptid", all = T)
setkey(visits_dt, ptid, visit)

# Add variable for model calendar time. Calendar time starts from the date of first study enrollment and proceeds in monthly intervals (= 365.25 days per year/12 months per year).
visits_dt[, first_enroll_dt := min(enrolldt)]
visits_dt[visit == 0, model_calendar_time := first_enroll_dt]
visits_dt[visit != 0, model_calendar_time := first_enroll_dt + ((365.25/12) * visit)]

# Create variable to indicate if model calendar time is on or after date of enrollment (by each woman)
visits_dt[, post_enrollment := ifelse(model_calendar_time >= enrolldt, 1, 0)]

# Create variable to indicate if model calendar time is on or after Sept. 1, 2013 (adherence interventions were implemented on 8/1/2013 and visits were monthly, so we wouldn't expect an effect of the adherence interventions to show up until at least 30 days later).
visits_dt[, post_adh_intervention := ifelse(model_calendar_time >= as.Date("2013-09-01", format = "%Y-%m-%d"), 1, 0)]

# Create variable to indicate if the participant visit is after three months of follow-up
visits_dt[, time_study_gt_three_mos := ifelse(cumsum(post_enrollment) >= 4, 1, 0), by = ptid]

# Create variable to indicate if model calendar time is after the date of early study termination (by woman). For women who ended the study per schedule, this indicator is 0 for all visits. For women who seroconverted, this indicator is 0 for all visits (women may have been more likely to terminate study participation following seroconversion, but were on study until seroconversion. We will treat this as an indication of willingness to continue participation in the absence of seroconversion).
visits_dt[, early_term_censor := ifelse(termrsn != "Sched. exit visit" & censor == 0, as.numeric(model_calendar_time > censordt), 0)]

# Check coding
visits_dt[censor == 1, all(early_term_censor == 0)] # Should be true
visits_dt[termrsn == "Sched. exit visit", all(early_term_censor == 0)] # Should be true
# View(visits_dt[termrsn != "Sched. exit visit" & censor == 0, .(ptid, visit, enrolldt, termrsn, censor, censordt, model_calendar_time, early_term_censor, post_enrollment)])

# Create variable to indicate observed visits on study (obs_on_study)
visits_dt[, visits_obs := round(termdays/30)]
visits_dt[, visits_post_enrollment := cumsum(post_enrollment), by = ptid]
visits_dt[, obs_on_study := as.numeric(visits_post_enrollment != 0 & visits_post_enrollment <= visits_obs)]

# There are three participants for whom the coding above assigns them 0 follow-up time. Manually adjust to ensure a minimum of 30 days of follow-up.
# View(visits_dt[, all(early_term_censor == post_enrollment), by = ptid][V1 == T])
visits_dt[ptid == 306300207 & visit == 11, post_enrollment := 1]
visits_dt[ptid == 312402776 & visit == 17, post_enrollment := 1]
visits_dt[ptid == 318303968 & visit == 21, post_enrollment := 1]

# Merge data table with all monthly visit dates with baseline variables. Carry forward values of baseline variables.
dt <- merge(x = dt[, .SD, .SDcols = which(names(dt) != "visit")], y = visits_dt[, .(ptid, visit, post_enrollment, early_term_censor, post_adh_intervention, time_study_gt_three_mos, obs_on_study)], by = "ptid", all.y = T)

# Calculate adherence probability
params_dt <- as.data.table(read_excel("~/Documents/code/aspire/parameters/parameters_adh.xlsx", range = "A1:B49", col_names = T))
params <- lapply(params_dt[, name], function(x) { x = params_dt[name == x, value] })
names(params) <- params_dt[, name]

dt[, log_odds_adh := 
         params$adh_coef_Intercept +
         params$adh_coef_post_adh_intervention * post_adh_intervention +
         params$adh_coef_time_study_gt_three_mos * time_study_gt_three_mos +
         params$`adh_coef_age_cat27-45` * (age >= 27) +
         params$adh_coef_fp_methodOther * (1 - bl_ocp) +
         params$adh_coef_pknow_studyYes * (pknow_study == "Yes") +
         params$`adh_coef_ring_worriesSomewhat or very worried` * (ring_worries == "Somewhat or very worried") +
         params$`adh_coef_same_pspYes or no primary partner` * (same_psp == "Yes or no primary partner") +
         params$`adh_coef_siteBotha's Hill` * (site == "South Africa: Botha's Hill") +
         params$`adh_coef_siteEmavundleni Cent` * (site == "South Africa: Emavundleni Cent") +
         params$adh_coef_siteeThekwini * (site == "South Africa: eThekwini") +
         params$adh_coef_siteIsipingo * (site == "South Africa: Isipingo") +
         params$adh_coef_siteLilongwe * (site == "Malawi: Lilongwe") +
         params$`adh_coef_siteMU-JHU` * (site == "Uganda: MU-JHU") +
         params$`adh_coef_siteRK Khan` * (site == "South Africa: RK Khan") +
         params$`adh_coef_siteSeke South` * (site == "Zimbabwe: Seke South") +
         params$adh_coef_siteSpilhaus * (site == "Zimbabwe: Spilhaus") +
         params$adh_coef_siteTongaat * (site == "South Africa: Tongaat") +
         params$adh_coef_siteUmkomaas * (site == "South Africa: Umkomaas") +
         params$adh_coef_siteVerulam * (site == "South Africa: Verulam") +
         params$adh_coef_siteWHRI * (site == "South Africa: WHRI") +
         params$adh_coef_siteZengeza * (site == "Zimbabwe: Zengeza")]

dt[, prob_adh := exp(log_odds_adh)/(1 + exp(log_odds_adh))]

# Each woman has up to three adherence values, with an increase in the probability of adherence after three months of study participation and following implementation of adherence interventions. Identify these adherence periods as 1, 2, or 3 for carry-forward values in model simulations.
adh_periods_dt <- dt[, .(prob_adh = sort(unique(prob_adh))), by = ptid]
adh_periods_dt[, incr_prob_adh := c(0, diff(prob_adh)), by = ptid]
adh_periods_dt[, adh_period := 1:.N, by = ptid]

dt <- merge(x = dt, y = adh_periods_dt, by = c("ptid", "prob_adh"), all.x = T)

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

# Create country-specific proportion of women reporting condom use at last sex. Use values for quarter one (numbers calculated by Elizabeth Brown).
dt_cond <- data.table(country = c("mal", "sa", "uga", "zim"), prop_condom = c(0.41, 0.52, 0.29, 0.45))

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

f_dt <- dt[, .(id, visit, post_enrollment, early_term_censor, obs_on_study, post_adh_intervention, time_study_gt_three_mos, prob_adh, incr_prob_adh, adh_period, ring_worries, pknow_study, same_psp, bl_ocp, country, site, arm, any_preg, f_age, f_age_cat, b_n_part, max_part, pp_vi_acts_avg, b_married, b_noalc, m_hiv, m_arv, b_n_sti, b_sti, b_bv, b_condom_lweek, prop_condom)]

save(f_dt, file = paste0(getwd(), "/f_dt.RData"))
